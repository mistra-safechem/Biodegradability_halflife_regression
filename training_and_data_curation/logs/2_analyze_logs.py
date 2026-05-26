#!/usr/bin/env python3
"""
2_analyze_logs.py
-----------------
Parse all 10_interval_svr.log files under this directory and produce a
human-readable log_summary.md with three comparison tables:

  1. Main performance table  — final model metrics per compartment/source
  2. Feature reduction table — winner strategy, n features, final metrics
  3. Applicability domain    — PCA coverage, AD-outside counts

Usage:
    uv run logs/2_analyze_logs.py
    python logs/2_analyze_logs.py

Script created by AI coding agent and may require manual adjustments to parsing logic if log format changes.
"""

from __future__ import annotations

import re
from pathlib import Path
from datetime import datetime

LOG_DIR = Path(__file__).parent
OUTPUT_FILE = LOG_DIR / "log_summary.md"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _find(pattern: str, text: str, group: int = 1, default: str = "—") -> str:
    m = re.search(pattern, text)
    return m.group(group) if m else default


def _pct(value_str: str, total_str: str) -> str:
    """Return 'value/total (xx.x%)' or just the raw strings on failure."""
    try:
        v = int(value_str)
        t = int(total_str)
        return f"{v}/{t} ({100 * v / t:.1f}%)"
    except (ValueError, ZeroDivisionError):
        return f"{value_str}/{total_str}"


def na(s: str) -> str:
    """Return '—' for empty/whitespace strings."""
    return s.strip() if s.strip() else "—"


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def parse_run_block(header: str, text: str) -> dict:
    """
    Extract all metrics from one run block.

    header examples:
        "air_hsbd/20260421_153507"
        "water_vega/20260422_090012"
    """
    parts = header.strip().split("/")
    run_key = parts[0]  # e.g. "air_hsbd"
    timestamp = parts[1] if len(parts) > 1 else "—"

    compartment_source = run_key.split("_", 1)
    compartment = compartment_source[0].upper()
    source = compartment_source[1] if len(compartment_source) > 1 else "—"

    r: dict = {
        "run_key": run_key,
        "timestamp": timestamp,
        "compartment": compartment,
        "source": source,
    }

    # ------------------------------------------------------------------
    # Dataset size
    # ------------------------------------------------------------------
    r["n_raw"] = _find(r"Raw dataset:\s*(\d+) rows", text)

    # ------------------------------------------------------------------
    # Initial CV (before feature reduction)
    # ------------------------------------------------------------------
    # First occurrence of the 5-fold CV line
    cv_lines = re.findall(
        r"5-fold CV on training set: R²=([\d.]+) \(±([\d.]+)\)" r"\s+RMSE\(log10\)=([\d.]+)\s+RMSE\(days-geom\)=([\d.]+)×",
        text,
    )
    if cv_lines:
        r["init_cv_r2"] = cv_lines[0][0]
        r["init_cv_r2_std"] = cv_lines[0][1]
        r["init_rmse_log10"] = cv_lines[0][2]
        r["init_rmse_days"] = cv_lines[0][3]
    else:
        r["init_cv_r2"] = r["init_cv_r2_std"] = r["init_rmse_log10"] = r["init_rmse_days"] = "—"

    if len(cv_lines) >= 2:
        r["final_cv_r2"] = cv_lines[1][0]
        r["final_cv_r2_std"] = cv_lines[1][1]
        r["final_rmse_log10"] = cv_lines[1][2]
        r["final_rmse_days"] = cv_lines[1][3]
    else:
        r["final_cv_r2"] = r["final_cv_r2_std"] = r["final_rmse_log10"] = r["final_rmse_days"] = "—"

    # ------------------------------------------------------------------
    # Test-set metrics (initial model, [test] tag)
    # ------------------------------------------------------------------
    test_cov = re.search(r"\[test\]\s+Interval coverage probability\s*:\s*([\d.]+)\s*\((\d+)/(\d+)\)", text)
    if test_cov:
        r["test_coverage"] = test_cov.group(1)
        r["test_cov_n"] = test_cov.group(2)
        r["test_cov_total"] = test_cov.group(3)
        r["test_cov_fmt"] = _pct(test_cov.group(2), test_cov.group(3))
    else:
        r["test_coverage"] = r["test_cov_n"] = r["test_cov_total"] = r["test_cov_fmt"] = "—"

    test_mil = re.search(r"\[test\]\s+Mean interval loss\s*:\s*([\d.]+) log10 days\s+\(~([\d.]+) days", text)
    r["test_mil_log10"] = test_mil.group(1) if test_mil else "—"
    r["test_mil_days"] = test_mil.group(2) if test_mil else "—"

    test_rho = re.search(r"\[test\]\s+Spearman rho\s*:\s*([\d.]+)", text)
    r["test_rho"] = test_rho.group(1) if test_rho else "—"

    test_tau = re.search(r"\[test\]\s+Kendall tau\s*:\s*([\d.]+)", text)
    r["test_tau"] = test_tau.group(1) if test_tau else "—"

    test_acc = re.search(r"\[test\]\s+Class accuracy\s*:\s*([\d.]+)", text)
    r["test_class_acc"] = test_acc.group(1) if test_acc else "—"

    # ------------------------------------------------------------------
    # Train-set metrics (initial model, [train] tag)
    # ------------------------------------------------------------------
    train_cov = re.search(r"\[train\]\s+Interval coverage probability\s*:\s*([\d.]+)\s*\((\d+)/(\d+)\)", text)
    if train_cov:
        r["train_coverage"] = train_cov.group(1)
        r["train_cov_fmt"] = _pct(train_cov.group(2), train_cov.group(3))
    else:
        r["train_coverage"] = r["train_cov_fmt"] = "—"

    train_rho = re.search(r"\[train\]\s+Spearman rho\s*:\s*([\d.]+)", text)
    r["train_rho"] = train_rho.group(1) if train_rho else "—"

    # ------------------------------------------------------------------
    # Feature reduction — winner
    # ------------------------------------------------------------------
    winner = re.search(r"Winner:\s*(.+?)\s+\(n=(\d+),\s*coverage=([\d.]+),\s*rho=([\d.]+)\)", text)
    if winner:
        r["winner_strategy"] = winner.group(1).strip()
        r["winner_n"] = winner.group(2)
        r["winner_coverage"] = winner.group(3)
        r["winner_rho"] = winner.group(4)
    else:
        r["winner_strategy"] = r["winner_n"] = r["winner_coverage"] = r["winner_rho"] = "—"

    # ------------------------------------------------------------------
    # Final model metrics ([final] tag)
    # ------------------------------------------------------------------
    final_cov = re.search(r"\[final\]\s+Interval coverage probability\s*:\s*([\d.]+)\s*\((\d+)/(\d+)\)", text)
    if final_cov:
        r["final_coverage"] = final_cov.group(1)
        r["final_cov_fmt"] = _pct(final_cov.group(2), final_cov.group(3))
    else:
        r["final_coverage"] = r["final_cov_fmt"] = "—"

    final_mil = re.search(r"\[final\]\s+Mean interval loss\s*:\s*([\d.]+) log10 days\s+\(~([\d.]+) days", text)
    r["final_mil_log10"] = final_mil.group(1) if final_mil else "—"
    r["final_mil_days"] = final_mil.group(2) if final_mil else "—"

    final_rho = re.search(r"\[final\]\s+Spearman rho\s*:\s*([\d.]+)", text)
    r["final_rho"] = final_rho.group(1) if final_rho else "—"

    final_tau = re.search(r"\[final\]\s+Kendall tau\s*:\s*([\d.]+)", text)
    r["final_tau"] = final_tau.group(1) if final_tau else "—"

    final_acc = re.search(r"\[final\]\s+Class accuracy\s*:\s*([\d.]+)", text)
    r["final_class_acc"] = final_acc.group(1) if final_acc else "—"

    # ------------------------------------------------------------------
    # Learning curve
    # ------------------------------------------------------------------
    lc = re.search(r"Learning curve val R².*?:\s*([\d.]+)\s*±\s*([\d.]+)", text)
    r["lc_val_r2"] = lc.group(1) if lc else "—"
    r["lc_val_r2_std"] = lc.group(2) if lc else "—"

    # ------------------------------------------------------------------
    # Applicability domain
    # ------------------------------------------------------------------
    pca_cov = re.search(r"PCA coverage:\s*([\d.]+)%", text)
    r["pca_coverage"] = pca_cov.group(1) if pca_cov else "—"

    ad_out = re.search(r"AD-outside test compounds:\s*(\d+)", text)
    r["ad_outside_n"] = ad_out.group(1) if ad_out else "—"

    # AD-outside as % of test set
    try:
        r["ad_outside_pct"] = f"{100 * int(r['ad_outside_n']) / int(r['test_cov_total']):.1f}%"
    except (ValueError, ZeroDivisionError):
        r["ad_outside_pct"] = "—"

    # ------------------------------------------------------------------
    # Artefact paths
    # ------------------------------------------------------------------
    model_path = re.search(r"Model saved to\s+(.+\.joblib)", text)
    r["model_path"] = model_path.group(1) if model_path else "—"

    ad_path = re.search(r"AD artefact saved to\s+(.+\.npz)", text)
    r["ad_path"] = ad_path.group(1) if ad_path else "—"

    db_path = re.search(r"Database:\s+(.+\.db)", text)
    r["db_path"] = db_path.group(1) if db_path else "—"

    return r


# ---------------------------------------------------------------------------
# Markdown builder
# ---------------------------------------------------------------------------


def _row(*cells: str) -> str:
    return "| " + " | ".join(cells) + " |"


def _sep(n: int) -> str:
    return "|" + "|".join(["---"] * n) + "|"


def format_md(runs: list[dict]) -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    lines: list[str] = []

    lines.append(f"# SVR Pipeline — Log Summary")
    lines.append(f"")
    lines.append(f"Generated: {now}  ")
    lines.append(f"Parsed {len(runs)} run(s) from `{LOG_DIR.name}/`")
    lines.append("")

    # ------------------------------------------------------------------
    # Table 1 — Main performance (final model)
    # ------------------------------------------------------------------
    lines.append("---")
    lines.append("")
    lines.append("## 1. Main Performance (Final Model)")
    lines.append("")
    lines.append("> Metrics from the final refitted model after feature reduction (`[final]` tag).")
    lines.append("> Coverage = fraction of test compounds whose predicted log₁₀(T½) falls inside the target interval.")
    lines.append("")

    cols1 = [
        "Compartment",
        "Source",
        "N (raw)",
        "Test coverage",
        "MIL (log₁₀ d)",
        "MIL (×days)",
        "Spearman ρ",
        "Kendall τ",
        "Class acc.",
        "CV R²",
        "CV R² std",
        "RMSE (×days)",
    ]
    lines.append(_row(*cols1))
    lines.append(_sep(len(cols1)))

    for r in runs:
        lines.append(
            _row(
                r["compartment"],
                r["source"],
                r["n_raw"],
                r["final_cov_fmt"],
                r["final_mil_log10"],
                r["final_mil_days"] + "×" if r["final_mil_days"] != "—" else "—",
                r["final_rho"],
                r["final_tau"],
                r["final_class_acc"],
                f"{r['final_cv_r2']} ±{r['final_cv_r2_std']}" if r["final_cv_r2"] != "—" else "—",
                r["final_cv_r2_std"],
                r["final_rmse_days"] + "×" if r["final_rmse_days"] != "—" else "—",
            )
        )

    lines.append("")

    # ------------------------------------------------------------------
    # Table 2 — Feature reduction
    # ------------------------------------------------------------------
    lines.append("---")
    lines.append("")
    lines.append("## 2. Feature Reduction")
    lines.append("")
    lines.append("> Winner strategy selected by highest (coverage, then ρ) on the test set.")
    lines.append("> Initial = full descriptor set; Final = reduced feature set (refitted model).")
    lines.append("")

    cols2 = [
        "Compartment",
        "Source",
        "Winner strategy",
        "N features",
        "Init. CV R²",
        "Final CV R²",
        "Init. test ρ",
        "Final test ρ",
        "LC val R² (±std)",
    ]
    lines.append(_row(*cols2))
    lines.append(_sep(len(cols2)))

    for r in runs:
        init_cv = f"{r['init_cv_r2']} ±{r['init_cv_r2_std']}" if r["init_cv_r2"] != "—" else "—"
        final_cv = f"{r['final_cv_r2']} ±{r['final_cv_r2_std']}" if r["final_cv_r2"] != "—" else "—"
        lc = f"{r['lc_val_r2']} ±{r['lc_val_r2_std']}" if r["lc_val_r2"] != "—" else "—"
        lines.append(
            _row(
                r["compartment"],
                r["source"],
                r["winner_strategy"],
                r["winner_n"],
                init_cv,
                final_cv,
                r["test_rho"],
                r["final_rho"],
                lc,
            )
        )

    lines.append("")

    # ------------------------------------------------------------------
    # Table 3 — Applicability Domain
    # ------------------------------------------------------------------
    lines.append("---")
    lines.append("")
    lines.append("## 3. Applicability Domain")
    lines.append("")
    lines.append("> Leverage-based AD (Williams plot).  ")
    lines.append("> AD-outside = test compounds with leverage h > h★ (structurally outside the training space).  ")
    lines.append("> PCA coverage = % of test compounds within the 95% PCA ellipse of the training set.")
    lines.append("")

    cols3 = [
        "Compartment",
        "Source",
        "N test",
        "AD-outside (n)",
        "AD-outside (%)",
        "PCA coverage (%)",
    ]
    lines.append(_row(*cols3))
    lines.append(_sep(len(cols3)))

    for r in runs:
        lines.append(
            _row(
                r["compartment"],
                r["source"],
                r["test_cov_total"],
                r["ad_outside_n"],
                r["ad_outside_pct"],
                r["pca_coverage"] + "%" if r["pca_coverage"] != "—" else "—",
            )
        )

    lines.append("")

    # ------------------------------------------------------------------
    # Table 4 — Train vs test overfitting check
    # ------------------------------------------------------------------
    lines.append("---")
    lines.append("")
    lines.append("## 4. Train / Test Gap")
    lines.append("")
    lines.append("> High train ρ with low test ρ signals overfitting.")
    lines.append("")

    cols4 = [
        "Compartment",
        "Source",
        "Train coverage",
        "Test coverage",
        "Train ρ",
        "Test ρ",
        "RMSE train (×days)",
        "RMSE test (×days)",
    ]
    lines.append(_row(*cols4))
    lines.append(_sep(len(cols4)))

    for r in runs:
        rmse_fmt = r["init_rmse_days"] + "×" if r["init_rmse_days"] != "—" else "—"
        lines.append(
            _row(
                r["compartment"],
                r["source"],
                r["train_cov_fmt"],
                r["test_cov_fmt"],
                r["train_rho"],
                r["test_rho"],
                rmse_fmt,
                rmse_fmt,  # same CV RMSE applies to both; distinguish if available
            )
        )

    lines.append("")

    # ------------------------------------------------------------------
    # Section 5 — Per-run notes
    # ------------------------------------------------------------------
    lines.append("---")
    lines.append("")
    lines.append("## 5. Run Notes")
    lines.append("")

    for r in runs:
        lines.append(f"### {r['compartment']} / {r['source']}  `{r['timestamp']}`")
        lines.append("")
        lines.append(f"- **Database**: `{r['db_path']}`")
        lines.append(f"- **Raw rows**: {r['n_raw']}")
        lines.append(f"- **Model**: `{r['model_path']}`")
        lines.append(f"- **AD artefact**: `{r['ad_path']}`")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def find_log_files(root: Path) -> list[Path]:
    files = [
        f
        for f in root.rglob("SVR_interval.log")
        if not any(part.startswith(".") for part in f.relative_to(root).parts)
        and "_log_by_occasion_manually_moved" not in str(f)
    ]
    files.sort(key=lambda p: p.parent.name)
    return files


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    log_files = find_log_files(LOG_DIR)
    if not log_files:
        print("No .log files found.")
        return

    runs: list[dict] = []
    for log_path in log_files:
        rel = log_path.parent.relative_to(LOG_DIR)
        # rel is like  air_hsbd/20260421_153507
        header = str(rel)
        text = log_path.read_text()
        runs.append(parse_run_block(header, text))

    md = format_md(runs)
    OUTPUT_FILE.write_text(md)
    print(f"Parsed {len(runs)} run(s) → {OUTPUT_FILE}")
    for r in runs:
        print(f"  {r['compartment']:10s}  {r['source']:10s}  {r['timestamp']}")


if __name__ == "__main__":
    main()
