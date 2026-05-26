"""Generate Markdown summaries from combined SVR log files.
Can be followed up by further summary with 4_summary.py

Script created by AI coding agent and may require manual adjustments to parsing logic if log format changes.
"""

from __future__ import annotations

import argparse
import ast
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Sequence


@dataclass
class DatasetResult:
    """Parsed metrics for a single dataset block."""

    key: str
    compartment: str
    dataset: str
    run_id: str | None

    # Full-feature model (initial train_model() run)
    cv_r2: float
    cv_r2_std: float | None
    full_test_metrics: dict  # {'R2':..., 'MAE':..., 'RMSE':..., ...}

    # Reduced-feature final model
    reduced_cv_r2: float | None
    reduced_cv_r2_std: float | None
    reduced_test_rmse_log: float | None  # RMSE in log10 units
    reduced_test_mae_log: float | None  # MAE in log10 units
    reduced_test_rmse_days: float | None
    reduced_test_mae_days: float | None

    # Feature reduction
    winner_strategy: str | None
    n_features_full: int | None
    n_features_reduced: int | None

    # AD / learning curve
    interpretation: str | None
    h_star: float | None
    train_in_ad: float | None
    test_in_ad: float | None
    bbox_coverage: float | None
    train_r2: float | None
    val_r2: float | None
    butina_clusters: int | None
    butina_singletons: int | None

    # Back-compat alias
    @property
    def test_metrics(self) -> dict:
        return self.full_test_metrics

    @property
    def name(self) -> str:
        return f"{self.compartment.title()} {self.dataset.upper()}".strip()

    # Full-feature held-out
    @property
    def test_r2(self) -> float:
        return float(self.full_test_metrics.get("R2", 0.0))

    @property
    def test_rmse(self) -> float:
        return float(self.full_test_metrics.get("RMSE", 0.0))

    @property
    def test_mae(self) -> float:
        return float(self.full_test_metrics.get("MAE", 0.0))

    # Reduced final model convenience
    @property
    def reduced_test_rmse(self) -> float:
        return self.reduced_test_rmse_log or 0.0

    @property
    def reduced_test_mae(self) -> float:
        return self.reduced_test_mae_log or 0.0


def parse_log_file(path: Path) -> List[DatasetResult]:
    results: List[DatasetResult] = []
    header: str | None = None
    body: List[str] = []

    for line in path.read_text().splitlines():
        stripped = line.rstrip("\n")
        if stripped.startswith("############################################################"):
            continue
        if stripped.startswith("# "):
            if header and body:
                results.append(_parse_block(header, body))
            header = stripped
            body = []
            continue
        if header:
            body.append(stripped)

    if header and body:
        results.append(_parse_block(header, body))

    return results


def _parse_block(header: str, lines: Sequence[str]) -> DatasetResult:
    trimmed = [ln.strip() for ln in lines if ln.strip()]
    header_value = header.lstrip("# ")
    block_key = header_value.split("/")[0]
    run_id = header_value.split("/")[1] if "/" in header_value else None
    if "_" in block_key:
        compartment, dataset = block_key.split("_", 1)
    else:
        compartment, dataset = block_key, ""

    # --- Full-feature section (before strategy comparison) ---
    full_lines = _lines_before_section(trimmed, "Feature Reduction — Strategy Comparison")
    cv_line = next((ln for ln in full_lines if ln.startswith("R² (CV):")), None)
    cv_r2, cv_std = _parse_value_with_uncertainty(cv_line)
    full_test_metrics = _extract_test_metrics(full_lines)

    # --- Reduced model section ---
    red_lines = _lines_in_section(trimmed, "Feature Reduction — Final Model")
    reduced_cv_r2, reduced_cv_r2_std = _extract_reduced_cv(red_lines)
    reduced_test_rmse_log = _extract_float_inline(red_lines, "RMSE (log10):")
    reduced_test_mae_log = _extract_float_inline(red_lines, "MAE (log10):")
    reduced_test_rmse_days = _extract_float_inline(red_lines, "RMSE (days):")
    reduced_test_mae_days = _extract_float_inline(red_lines, "MAE (days):")

    # --- Feature reduction strategy section ---
    strat_lines = _lines_in_section(trimmed, "Feature Reduction — Strategy Comparison")
    winner_strategy = _extract_winner(strat_lines)
    n_features_full = _extract_baseline_n(strat_lines)
    n_features_reduced = _extract_winner_n(strat_lines)

    # --- AD / learning curve (from reduced section + AD section) ---
    interpretation = _extract_after_prefix(trimmed, "Interpretation:")
    h_star = _extract_float_after(trimmed, "h* threshold")
    train_in_ad = _extract_percentage(trimmed, "Train inside AD")
    test_in_ad = _extract_percentage(trimmed, "Test  inside AD")
    if test_in_ad is None:
        test_in_ad = _extract_percentage(trimmed, "Test inside AD")
    bbox_coverage = _extract_percentage(trimmed, "Test coverage inside training bounding box")
    train_r2 = _extract_value(trimmed, "Training R²:")
    val_r2 = _extract_value(trimmed, "Validation R²:")
    butina_clusters, butina_singletons = _extract_butina(trimmed)

    return DatasetResult(
        key=block_key,
        compartment=compartment,
        dataset=dataset,
        run_id=run_id,
        cv_r2=cv_r2 or 0.0,
        cv_r2_std=cv_std,
        full_test_metrics=full_test_metrics,
        reduced_cv_r2=reduced_cv_r2,
        reduced_cv_r2_std=reduced_cv_r2_std,
        reduced_test_rmse_log=reduced_test_rmse_log,
        reduced_test_mae_log=reduced_test_mae_log,
        reduced_test_rmse_days=reduced_test_rmse_days,
        reduced_test_mae_days=reduced_test_mae_days,
        winner_strategy=winner_strategy,
        n_features_full=n_features_full,
        n_features_reduced=n_features_reduced,
        interpretation=interpretation,
        h_star=h_star,
        train_in_ad=train_in_ad,
        test_in_ad=test_in_ad,
        bbox_coverage=bbox_coverage,
        train_r2=train_r2,
        val_r2=val_r2,
        butina_clusters=butina_clusters,
        butina_singletons=butina_singletons,
    )


# ---------------------------------------------------------------------------
# Section slicing helpers
# ---------------------------------------------------------------------------


def _lines_before_section(lines: Sequence[str], section_marker: str) -> List[str]:
    """Return lines before the first line containing section_marker."""
    out = []
    for ln in lines:
        if section_marker in ln:
            break
        out.append(ln)
    return out


def _lines_in_section(lines: Sequence[str], section_marker: str) -> List[str]:
    """Return lines from section_marker header until next === separator or EOF.
    Skips the === separator line that immediately follows the marker header.
    """
    in_section = False
    skip_next_sep = False
    out = []
    for ln in lines:
        if section_marker in ln:
            in_section = True
            skip_next_sep = True
            continue
        if in_section:
            if ln.startswith("====") or ln.startswith("####"):
                if skip_next_sep:
                    skip_next_sep = False
                    continue  # skip opening separator
                break  # closing separator → end of section
            skip_next_sep = False
            out.append(ln)
    return out


# ---------------------------------------------------------------------------
# Existing parse helpers (unchanged)
# ---------------------------------------------------------------------------


def _parse_value_with_uncertainty(
    line: str | None,
) -> tuple[float | None, float | None]:
    if not line:
        return None, None
    try:
        value_part = line.split(":", 1)[1].strip()
        value_str, rest = value_part.split(" ", 1)
        value = float(value_str)
        std = None
        if "+/-" in rest:
            std = float(rest.split("+/-", 1)[1].strip(" )"))
        return value, std
    except Exception:
        return None, None


def _extract_test_metrics(lines: Sequence[str]) -> dict:
    for idx, line in enumerate(lines):
        if line.startswith("{'R2'") or line.startswith('{"R2"'):
            try:
                return ast.literal_eval(line)
            except Exception:
                continue
        if "Held-out test set performance" in line and idx + 1 < len(lines):
            metrics_line = lines[idx + 1]
            try:
                return ast.literal_eval(metrics_line)
            except Exception:
                continue
    return {}


def _extract_after_prefix(lines: Sequence[str], prefix: str) -> str | None:
    for line in lines:
        if line.startswith(prefix):
            return line.split(":", 1)[1].strip()
    return None


def _extract_float_after(lines: Sequence[str], prefix: str) -> float | None:
    for line in lines:
        if line.startswith(prefix):
            try:
                return float(line.split(":", 1)[1].strip())
            except ValueError:
                continue
    return None


def _extract_percentage(lines: Sequence[str], prefix: str) -> float | None:
    for line in lines:
        normalized = line.replace("  ", " ")
        if normalized.startswith(prefix):
            after = normalized.split(":", 1)[1].strip()
            after = after.replace("%", "")
            try:
                return float(after)
            except ValueError:
                continue
    return None


def _extract_value(lines: Sequence[str], prefix: str) -> float | None:
    for line in lines:
        if line.startswith(prefix):
            remainder = line.split(":", 1)[1].strip()
            number = remainder.split(" ")[0]
            try:
                return float(number)
            except ValueError:
                continue
    return None


# ---------------------------------------------------------------------------
# New parse helpers
# ---------------------------------------------------------------------------


def _extract_reduced_cv(lines: Sequence[str]) -> tuple[float | None, float | None]:
    """Parse 'R² (CV): 0.7903 (+/- 0.0556)' from reduced model section."""
    cv_line = next((ln for ln in lines if ln.startswith("R² (CV):")), None)
    return _parse_value_with_uncertainty(cv_line)


def _extract_float_inline(lines: Sequence[str], key: str) -> float | None:
    """Extract float from a line containing 'key value' possibly mid-line.
    e.g. '  RMSE (log10): 0.2612  MAE (log10): 0.1854'
    """
    for line in lines:
        if key in line:
            after = line.split(key, 1)[1].strip()
            token = after.split()[0]
            try:
                return float(token)
            except ValueError:
                continue
    return None


def _extract_winner(lines: Sequence[str]) -> str | None:
    """Parse 'Winner: A  top-SVR (importance>0)  (N=131, ...)' → strategy label."""
    for line in lines:
        if line.startswith("Winner:"):
            raw = line.split("Winner:", 1)[1].strip()
            # Strip trailing parenthetical (N=...) if present
            raw = re.sub(r"\s*\(N=\d+.*\)$", "", raw).strip()
            return raw
    return None


def _extract_baseline_n(lines: Sequence[str]) -> int | None:
    """Parse 'Baseline (full feature set, 188 features)' → 188."""
    for line in lines:
        m = re.search(r"Baseline \(full feature set,\s*(\d+)\s*features\)", line)
        if m:
            return int(m.group(1))
    return None


def _extract_winner_n(lines: Sequence[str]) -> int | None:
    """Parse N from 'Winner: ... (N=131, ...)' line."""
    for line in lines:
        if line.startswith("Winner:"):
            m = re.search(r"N=(\d+)", line)
            if m:
                return int(m.group(1))
    return None


def _extract_butina(lines: Sequence[str]) -> tuple[int | None, int | None]:
    """Parse 'Butina clusters: 29, singletons: 28 (93.3%)'."""
    for line in lines:
        m = re.search(r"Butina clusters:\s*(\d+),\s*singletons:\s*(\d+)", line)
        if m:
            return int(m.group(1)), int(m.group(2))
    return None, None


# ---------------------------------------------------------------------------
# Markdown builders
# ---------------------------------------------------------------------------


def _run_date(run_id: str | None) -> str:
    if run_id and len(run_id) >= 8:
        d = run_id[:8]
        return f"{d[:4]}-{d[4:6]}-{d[6:8]}"
    return "unknown date"


def build_markdown(results: Sequence[DatasetResult]) -> str:
    dates = ", ".join(sorted({_run_date(r.run_id) for r in results}))
    header = f"# Log review generated from SVR experiments ({dates}) recorded in logs/combined_logs.txt.\n"
    body_lines = [header, ""]
    for dataset in results:
        body_lines.append(_build_paragraph(dataset))
        body_lines.append("")
    body_lines.append(f"**Recommendations**: {build_recommendations(results)}")
    body_lines.append("")
    warns = build_warnings(results)
    if warns:
        body_lines.append("## Warnings")
        for w in warns:
            body_lines.append(f"- WARNING: {w}")
        body_lines.append("")
    body_lines.extend(build_table(results))
    body_lines.append("")
    return "\n".join(body_lines).strip() + "\n"


def _build_paragraph(result: DatasetResult) -> str:
    # Full-feature held-out
    full_r2 = _approx(result.test_r2)
    full_rmse = _format_with_commas(result.test_rmse)
    comparison = _comparison_phrase(result)

    # Reduced model
    feat_clause = _feat_clause(result)
    red_cv = (
        f"{result.reduced_cv_r2:.3f} ±{result.reduced_cv_r2_std:.4f}"
        if result.reduced_cv_r2 is not None and result.reduced_cv_r2_std is not None
        else (f"{result.reduced_cv_r2:.3f}" if result.reduced_cv_r2 is not None else "n/a")
    )
    red_rmse_days = _format_with_commas(result.reduced_test_rmse_days) if result.reduced_test_rmse_days is not None else "n/a"
    red_mae_days = _format_with_commas(result.reduced_test_mae_days) if result.reduced_test_mae_days is not None else "n/a"

    # AD / learning
    learning = _learning_phrase(result.interpretation)
    h_star = f"{result.h_star:.2f}" if result.h_star is not None else "n/a"
    test_ad = f"{result.test_in_ad:.1f}" if result.test_in_ad is not None else "n/a"
    ad_clause = _ad_phrase(result.test_in_ad)

    # Butina note
    butina_note = ""
    if result.butina_clusters is not None:
        butina_note = f" Butina clustering: {result.butina_clusters} clusters, {result.butina_singletons} singletons."

    # Val R² warning
    val_note = ""
    if result.val_r2 is not None and result.val_r2 < -1:
        val_note = f" Learning curve validation R² = {result.val_r2:.2f} — severe overfitting signal on small training subsets."

    return (
        f"**{result.name}**: Full-feature held-out R² ≈ {full_r2}, RMSE ≈ {full_rmse} days, "
        f"{comparison} "
        f"{feat_clause} "
        f"Reduced model: CV R² = {red_cv}, RMSE ≈ {red_rmse_days} days, MAE ≈ {red_mae_days} days. "
        f"Learning curves {learning}{val_note} "
        f"Williams plot: h* ≈ {h_star}, {test_ad}% of test set inside AD — {ad_clause}.{butina_note}"
    )


def _feat_clause(result: DatasetResult) -> str:
    full = result.n_features_full
    red = result.n_features_reduced
    strat = result.winner_strategy or "n/a"
    if full and red:
        return f"Feature reduction ({strat}) cut {full} → {red} features."
    return f"Feature reduction strategy: {strat}."


def _comparison_phrase(result: DatasetResult) -> str:
    delta = result.test_r2 - result.cv_r2
    cv_str = f"{result.cv_r2:.3f}" if result.cv_r2 else "n/a"
    if delta > 0.1:
        return f"notably stronger than 5-fold CV mean of {cv_str}, suggesting favorable split."
    if delta < -0.1:
        return f"below 5-fold CV mean of {cv_str}, indicating domain shift or overfitting."
    return f"closely aligned with 5-fold CV mean of {cv_str}."


def _learning_phrase(interpretation: str | None) -> str:
    if not interpretation:
        return "follow the recorded data trend."
    interp = interpretation.rstrip(".")
    lowered = interp.lower()
    if "still improving" in lowered:
        return "continue to improve — more data should still yield gains."
    if "plateau" in lowered:
        return "are plateauing — diminishing returns from more data."
    if "overfitting" in lowered:
        return "highlight overfitting risk — stronger regularization warranted."
    return f"indicate {interp}."


def _ad_phrase(test_in_ad: float | None) -> str:
    if test_in_ad is None:
        return "AD coverage not reported"
    if test_in_ad >= 70:
        return "most predictions trustworthy"
    if test_in_ad >= 40:
        return "sizable fraction of predictions remain extrapolations"
    if test_in_ad >= 15:
        return "many predictions fall outside reliable AD"
    return "vast majority of predictions fall outside reliable AD"


def build_warnings(results: Sequence[DatasetResult]) -> List[str]:
    warns: List[str] = []
    for r in results:
        if r.test_r2 < 0.2:
            warns.append(f"{r.name} — very low full-feature test R² ({r.test_r2:.3f}); model generalises poorly.")
        if r.test_in_ad is not None and r.test_in_ad < 50:
            warns.append(f"{r.name} — only {r.test_in_ad:.1f}% of test set inside AD; many predictions are extrapolations.")
        if r.val_r2 is not None and r.val_r2 < -1:
            warns.append(
                f"{r.name} — learning curve validation R² = {r.val_r2:.2f}; severe overfitting on small training subsets."
            )
        if r.reduced_cv_r2 is not None and (r.reduced_cv_r2 - r.test_r2) > 0.15:
            warns.append(
                f"{r.name} — reduced-model CV R² ({r.reduced_cv_r2:.3f}) vs "
                f"full-feature held-out R² ({r.test_r2:.3f}) gap > 0.15; "
                "test set may be out of distribution."
            )
    return warns


def build_table(results: Sequence[DatasetResult]) -> List[str]:
    lines: List[str] = []

    # --- Table 1: Full-feature model ---
    lines += [
        "## Full-Feature Model",
        "Compartment | Dataset | Test R² | Test RMSE (days) | Test MAE (days) | CV R² (5-fold)",
        "--- | --- | --- | --- | --- | ---",
    ]
    for r in results:
        cv_str = f"{r.cv_r2:.3f}"
        if r.cv_r2_std is not None:
            cv_str += f" ±{r.cv_r2_std:.4f}"
        lines.append(
            " | ".join(
                [
                    r.compartment.title(),
                    r.dataset.upper(),
                    f"{r.test_r2:.3f}",
                    _format_with_commas(r.test_rmse),
                    _format_with_commas(r.test_mae),
                    cv_str,
                ]
            )
        )

    lines.append("")

    # --- Table 2: Reduced-feature final model ---
    lines += [
        "## Reduced-Feature Final Model",
        "Compartment | Dataset | CV R² (5-fold) | CV R² std | RMSE log10 | MAE log10 | RMSE (days) | MAE (days) | Winner Strategy | Features (full→reduced) | h* | Test in AD | Train in AD | BB Coverage",
        "--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | ---",
    ]
    for r in results:
        red_cv = f"{r.reduced_cv_r2:.3f}" if r.reduced_cv_r2 is not None else "n/a"
        red_cv_std = f"±{r.reduced_cv_r2_std:.4f}" if r.reduced_cv_r2_std is not None else "n/a"
        red_rmse_log = f"{r.reduced_test_rmse_log:.4f}" if r.reduced_test_rmse_log is not None else "n/a"
        red_mae_log = f"{r.reduced_test_mae_log:.4f}" if r.reduced_test_mae_log is not None else "n/a"
        red_rmse_days = _format_with_commas(r.reduced_test_rmse_days) if r.reduced_test_rmse_days is not None else "n/a"
        red_mae_days = _format_with_commas(r.reduced_test_mae_days) if r.reduced_test_mae_days is not None else "n/a"
        winner = r.winner_strategy or "n/a"
        feat_range = f"{r.n_features_full}→{r.n_features_reduced}" if r.n_features_full and r.n_features_reduced else "n/a"
        h_star = f"{r.h_star:.2f}" if r.h_star is not None else "n/a"
        test_ad = f"{r.test_in_ad:.1f}%" if r.test_in_ad is not None else "n/a"
        train_ad = f"{r.train_in_ad:.1f}%" if r.train_in_ad is not None else "n/a"
        bb = f"{r.bbox_coverage:.1f}%" if r.bbox_coverage is not None else "n/a"

        lines.append(
            " | ".join(
                [
                    r.compartment.title(),
                    r.dataset.upper(),
                    red_cv,
                    red_cv_std,
                    red_rmse_log,
                    red_mae_log,
                    red_rmse_days,
                    red_mae_days,
                    winner,
                    feat_range,
                    h_star,
                    test_ad,
                    train_ad,
                    bb,
                ]
            )
        )

    return lines


def build_recommendations(results: Sequence[DatasetResult]) -> str:
    low_perf = [r.name for r in results if r.test_r2 < 0.3]
    low_ad = [r.name for r in results if r.test_in_ad is not None and r.test_in_ad < 50]
    big_gap = [r.name for r in results if r.reduced_cv_r2 is not None and (r.reduced_cv_r2 - r.test_r2) > 0.15]
    bad_val = [r.name for r in results if r.val_r2 is not None and r.val_r2 < -1]

    sentences: List[str] = []
    if low_perf:
        sentences.append(f"Gather more representative data for {_comma_list(low_perf)} — generalization weakest there.")
    if low_ad:
        sentences.append(
            f"AD coverage below 50% for {_comma_list(low_ad)} — rebalance train/test split or broaden descriptor coverage."
        )
    if big_gap:
        sentences.append(
            f"CV vs held-out gap for {_comma_list(big_gap)} — consider stronger regularization or alternative kernels."
        )
    if bad_val:
        sentences.append(
            f"Severe overfitting on learning curve for {_comma_list(bad_val)} — "
            f"collect more training data or increase regularization (C)."
        )
    if not sentences:
        sentences.append("No immediate concerns detected; continue monitoring future runs.")
    return " ".join(sentences)


def _comma_list(items: Sequence[str]) -> str:
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + f" and {items[-1]}"


def _approx(value: float, decimals: int = 2) -> str:
    formatted = f"{value:.{decimals}f}"
    return formatted.rstrip("0").rstrip(".")


def _format_with_commas(value: float) -> str:
    return f"{value:,.2f}".rstrip("0").rstrip(".")


def main(argv: Sequence[str] | None = None) -> int:
    # if no argv provided, argparse will default to sys.argv, simplest for automation
    log_dir = Path(__file__).resolve().parent
    default_log = log_dir / "combined_logs.txt"

    parser = argparse.ArgumentParser(description="Summarize SVR combined logs into Markdown")
    parser.add_argument("--log-path", type=Path, default=default_log, help="Path to combined_logs.txt")
    parser.add_argument(
        "--output-path",
        type=Path,
        default=None,
        help="Destination markdown file (default: log_review_{timestamp}.md)",
    )
    parser.add_argument(
        "--print",
        dest="print_output",
        action="store_true",
        help="Print markdown to stdout",
    )

    args = parser.parse_args(argv)

    results = parse_log_file(args.log_path)
    if not results:
        raise SystemExit(f"No dataset blocks parsed from {args.log_path}")

    # Derive timestamp from earliest run_id in the parsed results
    run_ids = sorted(r.run_id for r in results if r.run_id)
    ts = run_ids[0] if run_ids else "unknown"
    output_path: Path = args.output_path or (log_dir / f"log_review_{ts}.md")

    markdown = build_markdown(results)
    output_path.write_text(markdown)
    if args.print_output:
        print(markdown)
    else:
        print(f"Log review written to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
