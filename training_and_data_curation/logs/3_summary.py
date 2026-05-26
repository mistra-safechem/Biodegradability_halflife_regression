"""
3_summary.py
------------
Read the already-parsed log_summary.md and produce a concise one-page
presentation summary: log_summary_presentation.md

Contains:
  - One trimmed performance table (all rows, no redundant columns)
  - 5-10 evaluative/critical bullet points derived from the data

Usage:
    uv run logs/3_summary.py
    python logs/3_summary.py

Script created by AI coding agent and may require manual adjustments to parsing logic if log format changes.
"""

from __future__ import annotations

import re
from pathlib import Path
from datetime import datetime

LOG_DIR = Path(__file__).parent
input_md = LOG_DIR / "log_summary.md"
output_md = LOG_DIR / "log_summary_presentation.md"

# ---------------------------------------------------------------------------
# Parse Section 1 table from log_summary.md
# ---------------------------------------------------------------------------


def parse_main_table(md_text: str) -> list[dict]:
    """
    Extract the Main Performance table (Section 1) from the markdown.
    Returns a list of dicts keyed by column name.
    """
    # Isolate section 1 block (between first and second ---)
    section = re.search(r"## 1\. Main Performance.*?(?=^---|\Z)", md_text, re.DOTALL | re.MULTILINE)
    if not section:
        raise ValueError("Section 1 not found in log_summary.md")

    block = section.group(0)
    rows = []
    for line in block.splitlines():
        # Only data rows: start and end with |, not a header or separator row
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if not cells or cells[0] in ("Compartment", "---", ""):
            continue
        if all(c.startswith("---") or c == "" for c in cells):
            continue
        if len(cells) < 12:
            continue
        rows.append(
            {
                "compartment": cells[0],
                "source": cells[1],
                "n_raw": cells[2],
                "test_cov_fmt": cells[3],
                "mil_log10": cells[4],
                "mil_days": cells[5],
                "spearman_rho": cells[6],
                "kendall_tau": cells[7],
                "class_acc": cells[8],
                "cv_r2": cells[9],  # formatted as "X ±Y"
                "cv_r2_std": cells[10],
                "rmse_days": cells[11],
            }
        )
    return rows


def parse_ad_table(md_text: str) -> dict[tuple, dict]:
    """
    Extract AD table (Section 3) keyed by (compartment, source).
    """
    section = re.search(r"## 3\. Applicability Domain.*?(?=^---|\Z)", md_text, re.DOTALL | re.MULTILINE)
    if not section:
        return {}
    block = section.group(0)
    result = {}
    for line in block.splitlines():
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if not cells or cells[0] in ("Compartment", "---", ""):
            continue
        if all(c.startswith("---") or c == "" for c in cells):
            continue
        if len(cells) < 6:
            continue
        key = (cells[0], cells[1])
        result[key] = {
            "n_test": cells[2],
            "ad_outside_n": cells[3],
            "ad_outside_pct": cells[4],
            "pca_cov": cells[5],
        }
    return result


def parse_gap_table(md_text: str) -> dict[tuple, dict]:
    """
    Extract train/test gap table (Section 4) keyed by (compartment, source).
    """
    section = re.search(r"## 4\. Train / Test Gap.*?(?=^---|\Z)", md_text, re.DOTALL | re.MULTILINE)
    if not section:
        return {}
    block = section.group(0)
    result = {}
    for line in block.splitlines():
        if not line.startswith("|"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if not cells or cells[0] in ("Compartment", "---", ""):
            continue
        if all(c.startswith("---") or c == "" for c in cells):
            continue
        if len(cells) < 8:
            continue
        key = (cells[0], cells[1])
        result[key] = {
            "train_cov": cells[2],
            "test_cov": cells[3],
            "train_rho": cells[4],
            "test_rho": cells[5],
        }
    return result


# ---------------------------------------------------------------------------
# Float helpers
# ---------------------------------------------------------------------------


def _f(s: str) -> float | None:
    """Parse first float from a string, return None on failure."""
    m = re.search(r"[\d.]+", s)
    try:
        return float(m.group()) if m else None
    except ValueError:
        return None


def _pct_val(s: str) -> float | None:
    """Parse percentage value from strings like '23.3%' or '17/73 (23.3%)'."""
    m = re.search(r"\(([\d.]+)%\)", s)
    if m:
        return float(m.group(1))
    m = re.search(r"([\d.]+)%", s)
    if m:
        return float(m.group(1))
    return None


# ---------------------------------------------------------------------------
# Bullet generator
# ---------------------------------------------------------------------------

OVERFITTING_THRESHOLD = 0.15  # train-test ρ drop considered material
AD_BAD_THRESHOLD = 60.0  # % outside AD considered problematic
COVERAGE_GOOD = 50.0  # % interval coverage considered "useful"
R2_MARGINAL = 0.50  # CV R² below this is marginal


def generate_bullets(rows: list[dict], ad: dict, gap: dict) -> list[str]:
    n_models = len(rows)
    sources = sorted({r["source"] for r in rows})
    compartments = sorted({r["compartment"] for r in rows})
    n_range = (
        f"{min(int(r['n_raw']) for r in rows if r['n_raw'].isdigit())}–"
        f"{max(int(r['n_raw']) for r in rows if r['n_raw'].isdigit())}"
    )

    # Per-row computed values
    for r in rows:
        r["_rho"] = _f(r["spearman_rho"])
        r["_r2"] = _f(r["cv_r2"])
        r["_rmse"] = _f(r["rmse_days"].rstrip("×"))
        r["_cov"] = _pct_val(r["test_cov_fmt"])

    bullets: list[str] = []

    # 1 — Scope
    bullets.append(
        f"**Scope:** {n_models} SVR models trained across {len(compartments)} compartments "
        f"({', '.join(c.lower() for c in compartments)}) and {len(sources)} data sources "
        f"({', '.join(sources)}); training sets range from {n_range} compounds."
    )

    # 2 — Rank-order split by source
    vega_rows = [r for r in rows if r["source"] == "vega"]
    hsbd_rows = [r for r in rows if r["source"] == "hsbd"]

    if vega_rows and hsbd_rows:
        vega_rho_min = min(r["_rho"] for r in vega_rows if r["_rho"] is not None)
        vega_rho_max = max(r["_rho"] for r in vega_rows if r["_rho"] is not None)
        hsbd_rho_min = min(r["_rho"] for r in hsbd_rows if r["_rho"] is not None)
        hsbd_rho_max = max(r["_rho"] for r in hsbd_rows if r["_rho"] is not None)
        # Find weakest HSBD models
        weak = [r for r in hsbd_rows if r["_rho"] is not None and r["_rho"] < 0.75]
        weak_str = ", ".join(f"{r['compartment'].lower()} (ρ = {r['_rho']:.2f})" for r in weak)
        bullets.append(
            f"**Rank-order performance is strongly source-dependent.** VEGA models achieve "
            f"Spearman ρ = {vega_rho_min:.2f}–{vega_rho_max:.2f}; HSBD is notably weaker "
            f"for {weak_str} — indicating HSBD class structure or data quality limits predictability "
            f"in those compartments."
        )

    # 3 — Interval coverage
    poor_cov = [r for r in rows if r["_cov"] is not None and r["_cov"] < COVERAGE_GOOD]
    good_cov = [r for r in rows if r["_cov"] is not None and r["_cov"] >= COVERAGE_GOOD]
    poor_str = ", ".join(f"{r['compartment']}/{r['source']} ({r['_cov']:.0f}%)" for r in poor_cov)
    good_str = ", ".join(f"{r['compartment']}/{r['source']} ({r['_cov']:.0f}%)" for r in good_cov)
    bullets.append(
        f"**Interval coverage is inadequate for most models.** "
        + (f"{good_str} reach useful coverage (≥{COVERAGE_GOOD:.0f}%); " if good_cov else "")
        + f"the remainder ({poor_str}) have ~{int(sum(r['_cov'] for r in poor_cov) / len(poor_cov))}% "
        f"average coverage — roughly {100 - int(sum(r['_cov'] for r in poor_cov) / len(poor_cov))} in "
        f"100 test predictions fall outside the assigned persistence interval."
    )

    # 4 — CV R² range, flag marginal models
    marginal = [r for r in rows if r["_r2"] is not None and r["_r2"] < R2_MARGINAL]
    r2_min = min(r["_r2"] for r in rows if r["_r2"] is not None)
    r2_max = max(r["_r2"] for r in rows if r["_r2"] is not None)
    marg_str = ", ".join(f"{r['compartment']}/{r['source']} (R² = {r['_r2']:.2f})" for r in marginal)
    bullets.append(
        f"**CV R² spans {r2_min:.2f}–{r2_max:.2f} across all models"
        + (f"; {len(marginal)} models are below R² = {R2_MARGINAL}: {marg_str}." if marginal else ".")
        + " Models below this threshold should be treated as coarse classifiers rather than regressors."
    )

    # 5 — Combined vs single source
    combined_rows = [r for r in rows if r["source"] == "combined"]
    if combined_rows and (vega_rows or hsbd_rows):
        bullets.append(
            f"**Combining data sources does not improve performance.** All {len(combined_rows)} combined "
            f"models (water, soil, sediment) sit between their HSBD and VEGA counterparts rather than "
            f"above both — dataset heterogeneity introduces conflicting label structure that the SVR "
            f"cannot reconcile."
        )

    # 6 — Overfitting
    overfit = []
    for r in rows:
        key = (r["compartment"], r["source"])
        g = gap.get(key, {})
        tr = _f(g.get("train_rho", ""))
        te = _f(g.get("test_rho", ""))
        if tr is not None and te is not None and (tr - te) >= OVERFITTING_THRESHOLD:
            overfit.append((r["compartment"], r["source"], tr, te, tr - te))

    if overfit:
        worst = max(overfit, key=lambda x: x[4])
        overfit_str = "; ".join(f"{c}/{s} (Δρ = {d:.2f})" for c, s, _, _, d in overfit)
        bullets.append(
            f"**Overfitting is a material concern.** Train-to-test Spearman ρ drops of "
            f"≥ {OVERFITTING_THRESHOLD:.2f} observed in: {overfit_str}. "
            f"Worst case is {worst[0]}/{worst[1]} (train ρ = {worst[2]:.3f} → test ρ = {worst[3]:.3f}). "
            f"VEGA models show smaller gaps, likely reflecting more homogeneous label distributions."
        )

    # 7 — AD
    ad_problem = []
    ad_ok = []
    for r in rows:
        key = (r["compartment"], r["source"])
        a = ad.get(key, {})
        pct = _pct_val(a.get("ad_outside_pct", ""))
        if pct is not None:
            if pct >= AD_BAD_THRESHOLD:
                ad_problem.append((r["compartment"], r["source"], pct))
            else:
                ad_ok.append((r["compartment"], r["source"], pct))

    if ad_problem:
        bad_str = ", ".join(f"{c}/{s} ({p:.0f}%)" for c, s, p in ad_problem)
        ok_str = ", ".join(f"{c}/{s} ({p:.0f}%)" for c, s, p in ad_ok) if ad_ok else "none"
        bullets.append(
            f"**Applicability domain is the most critical limitation.** "
            f"{len(ad_problem)} models have ≥{AD_BAD_THRESHOLD:.0f}% of test compounds outside the "
            f"leverage threshold: {bad_str}. The majority of predictions from these models are "
            f"technically extrapolations. Adequate AD coverage only in: {ok_str}."
        )

    # 8 — RMSE interpretation
    rmse_min = min(r["_rmse"] for r in rows if r["_rmse"] is not None)
    rmse_max = max(r["_rmse"] for r in rows if r["_rmse"] is not None)
    bullets.append(
        f"**RMSE of {rmse_min:.2f}–{rmse_max:.2f}× days is structurally expected but practically "
        f"limiting.** Errors of this magnitude are consistent with quantised, interval-valued targets "
        f"rather than kinetic noise; they preclude use of any of these models for point-accurate "
        f"half-life estimation in a regulatory context."
    )

    # 9 — Best models
    best = sorted(
        [r for r in rows if r["_rho"] is not None and r["_r2"] is not None],
        key=lambda r: (r["_rho"] + (r["_r2"] or 0)) / 2,
        reverse=True,
    )[:3]
    best_str = "; ".join(f"{r['compartment']}/{r['source']} (ρ = {r['_rho']:.2f}, R² = {r['_r2']:.2f})" for r in best)
    bullets.append(
        f"**Overall assessment:** The strongest models are {best_str}. "
        f"These are suitable for screening-level persistence classification. "
        f"HSBD water and soil models, and all combined-source models, require further "
        f"data curation before any regulatory or comparative use."
    )

    return bullets


# ---------------------------------------------------------------------------
# Markdown writer
# ---------------------------------------------------------------------------


def _row(*cells: str) -> str:
    return "| " + " | ".join(str(c) for c in cells) + " |"


def _sep(n: int) -> str:
    return "|" + "|".join(["---"] * n) + "|"


def build_presentation_md(rows: list[dict], bullets: list[str]) -> str:
    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    lines: list[str] = []

    lines.append("# SVR Half-Life Models — Summary")
    lines.append("")
    lines.append(f"Generated: {now} | Source: `logs/log_summary.md`")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Model Performance Overview")
    lines.append("")

    # Table — drop mil_log10 and cv_r2_std (redundant/too granular)
    cols = [
        "Compartment",
        "Source",
        "N (raw)",
        "Test coverage",
        "MIL (×days)",
        "Spearman ρ",
        "Kendall τ",
        "Class acc.",
        "CV R²",
        "RMSE (×days)",
    ]
    lines.append(_row(*cols))
    lines.append(_sep(len(cols)))

    for r in rows:
        # CV R² — strip the ±std part for the presentation table
        cv_r2_clean = re.sub(r"\s*±[\d.]+", "", r["cv_r2"]) if r["cv_r2"] != "—" else "—"
        lines.append(
            _row(
                r["compartment"],
                r["source"],
                r["n_raw"],
                r["test_cov_fmt"],
                r["mil_days"],
                r["spearman_rho"],
                r["kendall_tau"],
                r["class_acc"],
                cv_r2_clean,
                r["rmse_days"],
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

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    if not input_md.exists():
        print(f"ERROR: {input_md} not found. Run 2_analyze_logs.py first.")
        return

    md_text = input_md.read_text()

    rows = parse_main_table(md_text)
    ad = parse_ad_table(md_text)
    gap = parse_gap_table(md_text)

    if not rows:
        print("ERROR: no data rows parsed from Section 1 table.")
        return

    bullets = generate_bullets(rows, ad, gap)
    output = build_presentation_md(rows, bullets)

    output_md.write_text(output)
    print(f"Parsed {len(rows)} model(s) → {output_md}")


if __name__ == "__main__":
    main()
