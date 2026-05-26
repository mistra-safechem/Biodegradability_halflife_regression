"""Logging utilities using Python's standard logging module.

Provides a factory function that returns a named Logger configured with both a
rotating file handler (persistent record) and an optional stream handler
(console echo).  All callers should obtain a logger via get_logger() rather
than instantiating handlers directly.
"""

import logging
import sys
from datetime import datetime
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Public constants
# ---------------------------------------------------------------------------
LOG_FORMAT = "[%(asctime)s] %(levelname)-8s %(name)s — %(message)s"
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
SECTION_CHAR = "="
SECTION_WIDTH = 70


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _make_section_str(title: str, char: str = SECTION_CHAR, width: int = SECTION_WIDTH) -> str:
    border = char * width
    return f"\n{border}\n  {title}\n{border}"


class _SectionLogger(logging.Logger):
    """Logger subclass with a .section() convenience method."""

    def section(self, title: str, char: str = SECTION_CHAR) -> None:
        """Log a clearly delimited section header at INFO level."""
        self.info(_make_section_str(title, char=char))


# Register the subclass so getLogger returns it for new loggers
logging.setLoggerClass(_SectionLogger)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_logger(
    name: str,
    log_dir: Optional[str | Path] = None,
    log_filename: Optional[str] = None,
    level: int = logging.DEBUG,
    console: bool = True,
    max_bytes: int = 10 * 1024 * 1024,  # 10 MB
    backup_count: int = 3,
) -> _SectionLogger:
    """Return a configured Logger instance.

    A new logger is created only if one with *name* does not already exist;
    subsequent calls with the same name return the cached instance unchanged.

    Parameters
    ----------
    name:
        Logger name — use a descriptive string such as ``"soil_hsbd"`` so that
        log messages are easy to filter.
    log_dir:
        Directory for the log file.  Created (including parents) if absent.
        When ``None`` no file handler is added.
    log_filename:
        Filename inside *log_dir*.  Defaults to ``"{name}.log"``.
    level:
        Minimum severity level captured by the logger (default: DEBUG).
    console:
        When ``True`` (default) a StreamHandler writing to stdout is attached.
        Set to ``False`` for batch / headless runs where stdout is not useful.
    max_bytes:
        Maximum size of the log file before rotation (default: 10 MB).
    backup_count:
        Number of rotated backup files to keep (default: 3).

    Returns
    -------
    _SectionLogger
        A Logger with an additional ``.section(title)`` method.

    Examples
    --------
    >>> log = get_logger("water_run", log_dir="logs/water_hsbd/20260421")
    >>> log.info("Preprocessing complete — %d samples", n)
    >>> log.section("Model Training")
    >>> log.warning("AD coverage below 80%%")
    """
    logger: _SectionLogger = logging.getLogger(name)  # type: ignore[assignment]

    # Avoid adding duplicate handlers on repeated calls
    if logger.handlers:
        return logger

    logger.setLevel(level)
    formatter = logging.Formatter(fmt=LOG_FORMAT, datefmt=DATE_FORMAT)

    # --- File handler ---
    if log_dir is not None:
        log_path = Path(log_dir)
        log_path.mkdir(parents=True, exist_ok=True)
        filename = log_filename or f"{name}.log"
        file_handler = RotatingFileHandler(
            log_path / filename,
            maxBytes=max_bytes,
            backupCount=backup_count,
            encoding="utf-8",
        )
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # --- Console handler ---
    if console:
        stream_handler = logging.StreamHandler(sys.stdout)
        stream_handler.setLevel(level)
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    return logger
