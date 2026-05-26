"""Logging utilities for appending output to logfiles.

LEGACY: copied to src/legacy for backwardscompatibility, see also the readme.md in this folder.
"""

from datetime import datetime
from pathlib import Path
from typing import Optional


def log_to_file(message: str, logfile: str, include_timestamp: bool = True) -> None:
    """
    Append a message to a logfile.

    Parameters
    ----------
    message : str
        The message to append to the logfile
    logfile : str
        Path to the logfile (will be created if it doesn't exist)
    include_timestamp : bool, optional
        Whether to prepend a timestamp to the message, by default True

    Examples
    --------
    >>> log_to_file("Starting model training", "run.log")
    >>> log_to_file("Training complete", "run.log", include_timestamp=True)
    """
    logfile_path = Path(logfile)

    # Create parent directories if they don't exist
    logfile_path.parent.mkdir(parents=True, exist_ok=True)

    # Format message with optional timestamp
    if include_timestamp:
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        formatted_message = f"[{timestamp}] {message}"
    else:
        formatted_message = message

    # Append to logfile
    with open(logfile_path, "a", encoding="utf-8") as f:
        f.write(formatted_message + "\n")


def log_section(title: str, logfile: str, char: str = "=") -> None:
    """
    Log a section header to a logfile.

    Parameters
    ----------
    title : str
        The section title
    logfile : str
        Path to the logfile
    char : str, optional
        Character to use for the border, by default "="

    Examples
    --------
    >>> log_section("Data Processing", "run.log")
    >>> log_section("Results", "run.log", char="-")
    """
    separator = char * 60
    log_to_file(separator, logfile, include_timestamp=False)
    log_to_file(f"  {title}", logfile, include_timestamp=True)
    log_to_file(separator, logfile, include_timestamp=False)
