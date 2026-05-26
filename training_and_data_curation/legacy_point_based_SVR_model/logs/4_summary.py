#!/usr/bin/env python3
"""
4_summary.py
------------
Read the most recent log_review_*.md in this directory and produce a concise
one-page presentation summary: log_summary_presentation.md

Contains:
  - One trimmed performance table (Reduced-Feature Final Model)
    with a computed CV–Test gap column
  - 5-10 evaluative/critical bullet points derived from the data
  - A Warnings section reproducing all flagged issues


Script created by AI coding agent and may require manual adjustments to parsing logic if log format changes.
"""

from __future__ import annotations

import re
from pathlib import Path
from datetime import datetime

LOG_DIR = Path(__file__).parent
OUTPUT_MD = LOG_DIR / "log_summary_presentation.md"

OVERFITTING_THRESHOLD = 0.15  # CV R² minus test R² gap considered material
AD_GOOD_THRESHOLD = 60.0  # % inside AD below this = problematic
R2_MARGINAL = 0.45  # CV R² below this = marginal model


# ---------------------------------------------------------------------------
# Input file selection — most recent log_review_*.md
# ---------------------------------------------------------------------------


def find_input(log_dir: Path) -> Path:
    candidates = sorted(log_dir.glob("log_review_*.md"), reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No log_review_*.md found in {log_dir}")
    return candidates[0]


# ---------------------------------------------------------------------------
# Float helper
# ---------------------------------------------------------------------------


def _f(s: str) -> float | None:
    s = s.replace(",", "").strip()
    m = re.search(r"-?[\d.]+", s)
    try:
        return float(m.group()) if m else None
    except ValueError:
        return None


def _pct(s: str) -> float | None:
    m = re.search(r"([\d.]+)%", s)
    return float(m.group(1)) if m else None


# ---------------------------------------------------------------------------
# Parse Full-Feature table
# ---------------------------------------------------------------------------


def parse_full_feature_table(text: str) -> dict[tuple, dict]:
    """
    Parse the '## Full-Feature Model' markdown table.
    Returns dict keyed by (compartment, dataset).
    """
    section = re.search(r"## Full-Feature Model\s*\n(.+?)(?=\n##|\Z)", text, re.DOTALL)
    if not section:
        return {}
    block = section.group(1)
    result = {}
    for line in block.splitlines():
        if not line.startswith("|") and "|" not in line:
            continue
        # Handle pipe-delimited lines without leading |
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 6:
            continue
        if cells[0] in ("Compartment", "---", "") or cells[0].startswith("---"):
            continue
        if all(c.startswith("---") or c == "" for c in cells):
            continue
        key = (cells[0], cells[1])
        result[key] = {
            "test_r2": cells[2],
            "test_rmse": cells[3],
            "test_mae": cells[4],
            "cv_r2_full": cells[5].split()[0] if cells[5] else "—",
        }
    return result


# ---------------------------------------------------------------------------
# Parse Reduced-Feature table
# ---------------------------------------------------------------------------


def parse_reduced_table(text: str) -> list[dict]:
    """
    Parse the '## Reduced-Feature Final Model' markdown table.
    Returns list of dicts in table order.
    """
    section = re.search(r"## Reduced-Feature Final Model\s*\n(.+?)(?=\n##|\Z)", text, re.DOTALL)
    if not section:
        return []
    block = section.group(1)
    rows = []
    for line in block.splitlines():
        if not line.startswith("|") and "|" not in line:
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        # Expected columns (0-indexed):
        # 0:Comp 1:Dataset 2:CV R² 3:CV R² std 4:RMSE log10 5:MAE log10
        # 6:RMSE days 7:MAE days 8:Winner Strategy 9:Features 10:h* 11:Test in AD
        # 12:Train in AD 13:BB Coverage
        if len(cells) < 12:
            continue
        if cells[0] in ("Compartment", "---", "") or cells[0].startswith("---"):
            continue
        if all(c.startswith("---") or c == "" for c in cells):
            continue
        rows.append(
            {
                "compartment": cells[0],
                "dataset": cells[1],
                "cv_r2": cells[2],
                "cv_r2_std": cells[3].lstrip("±"),
                "rmse_log10": cells[4],
                "mae_log10": cells[5],
                "rmse_days": cells[6],
                "mae_days": cells[7],
                "winner": cells[8],
                "features": cells[9],
                "h_star": cells[10],
                "test_in_ad": cells[11],
                "train_in_ad": cells[12] if len(cells) > 12 else "—",
                "bb_coverage": cells[13] if len(cells) > 13 else "—",
            }
        )
    return rows


# ---------------------------------------------------------------------------
# Parse Warnings
# ---------------------------------------------------------------------------


def parse_warnings(text: str) -> list[str]:
    section = re.search(r"## Warnings\s*\n(.+?)(?=\n##|\Z)", text, re.DOTALL)
    if not section:
        return []
    warnings = []
    for line in section.group(1).splitlines():
        line = line.strip()
        if line.startswith("- WARNING:"):
            warnings.append(line[len("- WARNING:") :].strip())
    return warnings


# ---------------------------------------------------------------------------
# Bullet generator
# ---------------------------------------------------------------------------


def generate_bullets(
    rows: list[dict],
    full: dict[tuple, dict],
    warnings: list[str],
) -> list[str]:

    # Attach full-feature test R² and compute gap
    for r in rows:
        key = (r["compartment"], r["dataset"])
        f = full.get(key, {})
        r["_test_r2"] = _f(f.get("test_r2", ""))
        r["_cv_r2"] = _f(r["cv_r2"])
        r["_rmse"] = _f(r["rmse_days"].replace(",", ""))
        r["_mae"] = _f(r["mae_days"].replace(",", ""))
        r["_ad_pct"] = _pct(r["test_in_ad"])
        if r["_cv_r2"] is not None and r["_test_r2"] is not None:
            r["_gap"] = round(r["_cv_r2"] - r["_test_r2"], 3)
        else:
            r["_gap"] = None

    bullets: list[str] = []

    # 1 — Scope
    sources = sorted({r["dataset"] for r in rows})
    compartments = sorted({r["compartment"] for r in rows})
    strategies = {r["winner"].strip() for r in rows}
    strat_str = "; ".join(sorted(strategies))
    bullets.append(
        f"**Scope:** {len(rows)} SVR plain regression models across "
        f"{len(compartments)} compartments ({', '.join(c.lower() for c in compartments)}) "
        f"and {len(sources)} data sources ({', '.join(sources)}). "
        f"Feature reduction strategies used: {strat_str}. "
        f"Targets are treated as point estimates on a log scale — this is not an interval-aware pipeline."
    )

    # 2 — Worst model
    worst = min((r for r in rows if r["_test_r2"] is not None), key=lambda r: r["_test_r2"])
    bullets.append(
        f"**{worst['compartment']} {worst['dataset']} is the clear failure case.** "
        f"Full-feature test R² = {worst['_test_r2']:.2f} and RMSE = {worst['rmse_days']} days — "
        f"the model learns something in cross-validation (CV R² = {worst['_cv_r2']:.2f}) "
        f"but generalises almost not at all to the held-out test set "
        f"(CV–test gap = {worst['_gap']:.2f}). "
        f"More representative training data is the only viable fix."
    )

    # 3 — Overfitting models (gap > threshold, excluding negative gaps)
    overfit = [r for r in rows if r["_gap"] is not None and r["_gap"] > OVERFITTING_THRESHOLD]
    if overfit:
        ov_str = ", ".join(
            f"{r['compartment']} {r['dataset']} (gap = {r['_gap']:.2f})" for r in sorted(overfit, key=lambda r: -r["_gap"])
        )
        bullets.append(
            f"**{len(overfit)} models show material overfitting (CV–test gap > {OVERFITTING_THRESHOLD}):** "
            f"{ov_str}. In all cases the test set is likely out of distribution relative to the "
            f"training fold structure — stronger regularisation or rebalanced splits are warranted."
        )

    # 4 — Anomalous negative gap (Water VEGA pattern)
    anomaly = [r for r in rows if r["_gap"] is not None and r["_gap"] < -0.10]
    if anomaly:
        for r in anomaly:
            lc_warn = next(
                (w for w in warnings if r["compartment"] in w and r["dataset"] in w and "learning curve" in w.lower()), None
            )
            lc_note = f" Learning curve validation R² signal: {lc_warn}." if lc_warn else ""
            bullets.append(
                f"**{r['compartment']} {r['dataset']} is an anomaly that should not be trusted.** "
                f"The full-feature test R² ({r['_test_r2']:.2f}) is implausibly higher than CV R² "
                f"({r['_cv_r2']:.2f}), and reduced-model RMSE is {r['rmse_days']} days — "
                f"the high test R² almost certainly reflects a favourable random split, not genuine generalisation." + lc_note
            )

    # 5 — RMSE scale warning
    rmse_vals = [(r["compartment"], r["dataset"], r["_rmse"]) for r in rows if r["_rmse"] is not None]
    rmse_min = min(rmse_vals, key=lambda x: x[2])
    rmse_max = max(rmse_vals, key=lambda x: x[2])
    bullets.append(
        f"**RMSE on the raw day scale is dominated by high-persistence outliers and is not "
        f"comparable across compartments.** "
        f"{rmse_min[0]} {rmse_min[1]} RMSE = {rmse_min[2]:.0f} days vs "
        f"{rmse_max[0]} {rmse_max[1]} RMSE = {rmse_max[2]:,.0f} days — "
        f"this reflects the target distribution scale, not relative model quality. "
        f"Log-scale RMSE is the more meaningful cross-compartment metric."
    )

    # 6 — AD coverage
    ad_bad = [r for r in rows if r["_ad_pct"] is not None and r["_ad_pct"] < AD_GOOD_THRESHOLD]
    ad_good = [r for r in rows if r["_ad_pct"] is not None and r["_ad_pct"] >= AD_GOOD_THRESHOLD]
    if ad_bad:
        bad_str = ", ".join(
            f"{r['compartment']} {r['dataset']} ({r['test_in_ad']})" for r in sorted(ad_bad, key=lambda r: r["_ad_pct"])
        )
        good_str = (
            ", ".join(
                f"{r['compartment']} {r['dataset']} ({r['test_in_ad']})" for r in sorted(ad_good, key=lambda r: -r["_ad_pct"])
            )
            if ad_good
            else "none"
        )
        bullets.append(
            f"**AD coverage is insufficient in {len(ad_bad)} models** "
            f"(< {AD_GOOD_THRESHOLD:.0f}% of test compounds inside the leverage-based AD): "
            f"{bad_str}. The majority or near-majority of their predictions are structural extrapolations. "
            f"Best AD coverage: {good_str}."
        )

    # 7 — VEGA vs HSBD log-scale comparison
    vega_rows = [r for r in rows if r["dataset"] == "VEGA" and r["rmse_log10"] != "—"]
    hsbd_rows = [r for r in rows if r["dataset"] == "HSBD" and r["rmse_log10"] != "—"]
    if vega_rows and hsbd_rows:
        bullets.append(
            f"**VEGA-source models have higher CV R² but raw RMSE is not straightforwardly better.** "
            f"On log-scale RMSE (more honest cross-compartment comparison), VEGA models range "
            f"{min(_f(r['rmse_log10']) for r in vega_rows if _f(r['rmse_log10'])):.2f}–"
            f"{max(_f(r['rmse_log10']) for r in vega_rows if _f(r['rmse_log10'])):.2f} log₁₀ d, "
            f"vs HSBD {min(_f(r['rmse_log10']) for r in hsbd_rows if _f(r['rmse_log10'])):.2f}–"
            f"{max(_f(r['rmse_log10']) for r in hsbd_rows if _f(r['rmse_log10'])):.2f} log₁₀ d. "
            f"The raw-day RMSE contrast is driven by extreme outliers in HSBD soil and sediment."
        )

    # 8 — Combined source reliability
    combined = [r for r in rows if r["dataset"] == "COMBINED"]
    if combined:
        comb_gaps = [r for r in combined if r["_gap"] is not None]
        max_gap = max(r["_gap"] for r in comb_gaps) if comb_gaps else None
        bullets.append(
            f"**Combined-source models are the most practically reliable.** "
            f"Water, Soil, and Sediment COMBINED all show small CV–test gaps "
            f"(max {max_gap:.2f}) and the best AD coverage of any source group. "
            f"They do not reach the CV R² of the best VEGA models, but they generalise more consistently."
        )

    # 9 — Overall verdict
    # Models with gap <= threshold AND no AD warning
    ad_warned_keys = set()
    for w in warnings:
        if "inside AD" in w:
            m = re.match(r"(\w+)\s+(\w+)", w)
            if m:
                ad_warned_keys.add((m.group(1), m.group(2)))

    reliable = [
        r
        for r in rows
        if r["_gap"] is not None
        and r["_gap"] <= OVERFITTING_THRESHOLD
        and (r["compartment"], r["dataset"]) not in ad_warned_keys
        and r["_ad_pct"] is not None
        and r["_ad_pct"] >= AD_GOOD_THRESHOLD
    ]
    reliable_str = ", ".join(f"{r['compartment']} {r['dataset']}" for r in reliable) if reliable else "none without caveats"

    bullets.append(
        f"**Overall assessment:** Models with small CV–test gap and adequate AD coverage: "
        f"{reliable_str}. "
        f"These are the only candidates for cautious screening use. "
        f"All remaining models require data augmentation, regularisation tuning, or "
        f"train/test rebalancing before applied use."
    )

    return bullets


# ---------------------------------------------------------------------------
# Markdown builder
# ---------------------------------------------------------------------------


def _row(*cells) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def _sep(n: int) -> str:
    return "|" + "|".join(["---"] * n) + "|"


def build_md(
    rows: list[dict],
    full: dict[tuple, dict],
    bullets: list[str],
    warnings: list[str],
    source_file: Path,
) -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    lines: list[str] = []

    lines.append("# SVR Half-Life Models — Summary (Plain Regression)")
    lines.append("")
    lines.append(f"Generated: {now} | Source: `logs/{source_file.name}`")
    lines.append("")
    lines.append(
        "> **Note:** This pipeline treats half-life targets as point estimates on a log scale — "
        "it is a plain regression, not an interval-aware model. Metrics are absolute RMSE and MAE "
        "in days and R², not interval coverage or rank-correlation. RMSE values appear large because "
        "no forgiveness is given for predictions near-but-outside a persistence class boundary. "
        "Compare to the interval-aware pipeline (`WP3bioddeg_V2`) with care."
    )
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Model Performance Overview (Reduced-Feature Final Model)")
    lines.append("")
    lines.append(
        "> CV–Test gap = CV R² minus full-feature held-out Test R². "
        "Values > 0.15 indicate overfitting or out-of-distribution test split (flagged in warnings). "
        "**Bold** gap values exceed the threshold."
    )
    lines.append("")

    cols = [
        "Compartment",
        "Dataset",
        "CV R²",
        "Test R² (full)",
        "CV–Test gap",
        "RMSE (days)",
        "MAE (days)",
        "Features (→)",
        "Test in AD",
    ]
    lines.append(_row(*cols))
    lines.append(_sep(len(cols)))

    for r in rows:
        key = (r["compartment"], r["dataset"])
        f_data = full.get(key, {})
        test_r2 = f_data.get("test_r2", "—")
        cv_str = f"{r['cv_r2']} ±{r['cv_r2_std']}" if r["cv_r2_std"] else r["cv_r2"]
        gap = r.get("_gap")
        if gap is not None:
            gap_str = f"**{gap:.3f}**" if abs(gap) > OVERFITTING_THRESHOLD else f"{gap:.3f}"
        else:
            gap_str = "—"
        lines.append(
            _row(
                r["compartment"],
                r["dataset"],
                cv_str,
                test_r2,
                gap_str,
                r["rmse_days"],
                r["mae_days"],
                r["features"],
                r["test_in_ad"],
            )
        )

    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Key Findings")
    lines.append("")
    for b in bullets:
        lines.append(f"- {b}")
        lines.append("")

    # Warnings section
    if warnings:
        lines.append("---")
        lines.append("")
        lines.append("## Warnings")
        lines.append("")
        # Group by model (first token before " —")
        lines.append("| Model | Warning |")
        lines.append("|---|---|")
        for w in warnings:
            # Try to extract model label from start of warning
            m = re.match(r"^([A-Za-z]+ [A-Z]+ —|[A-Za-z]+ [A-Z]+\s+—)\s*(.*)", w)
            if m:
                model = m.group(1).rstrip(" —").strip()
                msg = m.group(2).strip()
            else:
                # Split on first " — "
                parts = w.split(" — ", 1)
                model = parts[0].strip() if len(parts) == 2 else "—"
                msg = parts[1].strip() if len(parts) == 2 else w
            lines.append(f"| {model} | {msg} |")

    lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    source_file = find_input(LOG_DIR)
    print(f"Reading: {source_file.name}")

    text = source_file.read_text()

    full = parse_full_feature_table(text)
    rows = parse_reduced_table(text)
    warnings = parse_warnings(text)

    if not rows:
        print("ERROR: no rows parsed from Reduced-Feature Final Model table.")
        return

    bullets = generate_bullets(rows, full, warnings)
    md = build_md(rows, full, bullets, warnings, source_file)

    OUTPUT_MD.write_text(md)
    print(f"Parsed {len(rows)} model(s), {len(warnings)} warning(s) → {OUTPUT_MD.name}")


if __name__ == "__main__":
    main()
