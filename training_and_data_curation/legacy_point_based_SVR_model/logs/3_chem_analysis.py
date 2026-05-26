#!/usr/bin/env python3
"""Generate per-run chemistry_analysis_interpretation.md (in each sublog folder) files.

Recursively scans the logs directory for chemistry_analysis/ subdirectories
and writes an interpretation .md file into each one (skipping folders that
already have one).  Mirrors the structure of the sediment_vega reference file.

Script created by AI coding agent and may require manual adjustments to parsing logic if log format changes.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

LOG_DIR = Path(__file__).parent
OUTPUT_FILENAME = "chemistry_analysis_interpretation.md"
EXCLUDE_DIR = "_log_by_occasion_manually_moved"


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def find_chemistry_dirs(root: Path) -> list[Path]:
    """Return all chemistry_analysis/ subdirectories, excluding archived ones."""
    dirs = [p for p in root.rglob("chemistry_analysis") if p.is_dir()]
    dirs = [d for d in dirs if EXCLUDE_DIR not in str(d)]
    dirs.sort(key=lambda p: str(p))
    return dirs


# ---------------------------------------------------------------------------
# Safe readers
# ---------------------------------------------------------------------------


def _read_csv(path: Path) -> list[dict]:
    """Return list-of-dicts from CSV, or [] if file missing."""
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _read_json(path: Path) -> dict:
    """Return dict from JSON, or {} if file missing."""
    if not path.exists():
        return {}
    with path.open() as fh:
        return json.load(fh)


def _csv_to_fenced(rows: list[dict]) -> str:
    """Render list-of-dicts back to a fenced CSV code block."""
    if not rows:
        return "```\n(no data)\n```"
    headers = list(rows[0].keys())
    lines = [",".join(headers)]
    for row in rows:
        lines.append(",".join(row.get(h, "") for h in headers))
    return "```\n" + "\n".join(lines) + "\n```"


# ---------------------------------------------------------------------------
# Rating helpers
# ---------------------------------------------------------------------------


def _rate_r2(r2: float, n: int) -> tuple[str, str]:
    """Return (emoji, label) for R² value, considering sample size."""
    if n <= 2:
        return "⚠️", "UNRELIABLE (n≤2)"
    if r2 > 0.5:
        return "✅", "GOOD"
    if r2 >= 0.2:
        return "⚠️", "MEDIUM"
    return "❌", "BAD"


def _rate_rmse(rmse: float) -> tuple[str, str]:
    if rmse < 200:
        return "✅", "GOOD"
    if rmse <= 500:
        return "⚠️", "MEDIUM"
    return "❌", "BAD"


def _rate_ad(pct_outside: float) -> tuple[str, str]:
    if pct_outside < 20:
        return "✅", "GOOD"
    if pct_outside <= 50:
        return "⚠️", "CONCERNING"
    return "❌", "BAD"


def _rate_stability(stable: str, pred_std: float) -> tuple[str, str]:
    if stable.strip().lower() == "true" and pred_std < 200:
        return "✅", "GOOD"
    if stable.strip().lower() == "true":
        return "⚠️", "MEDIUM"
    if pred_std > 500:
        return "❌", "BAD"
    return "⚠️", "MEDIUM"


def _scaffold_label(smiles: str) -> str:
    """Return a human-readable label for common scaffolds.
    This is somewhat arbitrary but sufficient for this purpose.

    """
    _map = {
        "c1ccccc1": "benzene",
        "c1ccc(-c2ccccc2)cc1": "biphenyl",
        "c1ccc2c(c1)Oc1ccccc1O2": "diphenyl ether (fused)",
        "c1ccc(Oc2ccccc2)cc1": "diphenyl ether",
        "C1CCCCC1": "cyclohexane",
        "c1ccc2ccccc2c1": "naphthalene",
        "c1ccncc1": "pyridine",
        "c1ccc2sncc2c1": "benzothiazole",
        "c1ccsc1": "thiophene",
        "": "non-ring / acyclic",
    }
    return _map.get(smiles, smiles[:40] + ("..." if len(smiles) > 40 else ""))


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------


def _section_header(chem_dir: Path) -> str:
    """Build the file header with the run-specific path."""
    # Derive compartment/dataset/timestamp from path
    # Expected: .../logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/
    try:
        timestamp_dir = chem_dir.parent
        run_dir = timestamp_dir.parent
        rel = chem_dir.relative_to(LOG_DIR)
    except ValueError:
        rel = chem_dir
    _run_parts = rel.parts
    lines = [
        "# Chemistry Analysis Interpretation Guide",
        "## File Structure",
        "All chemistry analysis outputs are located in:",
        "```",
        "logs/{compartment}_{dataset}/{timestamp}/chemistry_analysis/",
        "```",
        f"For this run:",
        "```",
        str(f"logs/{_run_parts[0]}/{_run_parts[1]}/{_run_parts[2]}/"),
        "```",
        "---",
        "",
    ]
    return "\n".join(lines)


def _section1_scaffold(chem_dir: Path) -> str:
    dist_rows = _read_csv(chem_dir / "scaffold_distribution.csv")
    perf_rows = _read_csv(chem_dir / "scaffold_performance.csv")

    lines: list[str] = ["## Section 1: Scaffold Analysis", ""]

    # --- scaffold_distribution.csv ---
    lines += [
        "### `scaffold_distribution.csv`",
        "",
        "**Purpose:** Shows the distribution of Murcko scaffolds (core ring structures) across your dataset.",
        "",
        "**Columns:**",
        "| Column | Description | Type |",
        "|--------|-------------|------|",
        "| `scaffold_smiles` | SMILES representation of the scaffold | string |",
        "| `count` | Number of molecules containing this scaffold | int |",
        "| `coverage_pct` | Percentage of total dataset | float (0-100) |",
        "| `cumulative_pct` | Cumulative percentage (sorted by count) | float (0-100) |",
        "",
    ]

    if dist_rows:
        n_total = sum(int(r["count"]) for r in dist_rows)
        n_unique = len(dist_rows)
        top = dist_rows[0]
        top_label = _scaffold_label(top["scaffold_smiles"])
        top_pct = float(top["coverage_pct"])
        top3_cumulative = float(dist_rows[2]["cumulative_pct"]) if len(dist_rows) >= 3 else None
        has_empty = any(r["scaffold_smiles"].strip() == "" for r in dist_rows)

        lines += [
            "**Example data (top rows):**",
            _csv_to_fenced(dist_rows[:6]),
            "",
            "**Interpretation:**",
            "",
            "| Metric | Good | Bad | Action |",
            "|--------|------|-----|--------|",
            "| **Top scaffold coverage** | 20-40% | >60% or <15% | Broaden dataset if too narrow/diverse |",
            "| **Unique scaffolds** | 10-30 for 50-200 molecules | <5 or >50 | Consider stratification |",
            "| **Singleton scaffolds** | <30% appear once | >50% appear once | May need more data |",
            "",
            f"**Findings ({chem_dir.parent.parent.name}):**",
        ]

        # n_unique assessment
        if n_unique < 5:
            diversity_note = f"Very low scaffold diversity — only {n_unique} unique scaffolds."
        elif n_unique > 50:
            diversity_note = f"Very high scaffold diversity — {n_unique} unique scaffolds (many singletons expected)."
        else:
            diversity_note = f"{n_unique} unique scaffolds."

        lines += [
            f"- {n_total} total molecules, {diversity_note}",
            f"- Top scaffold: {top_label} (`{top['scaffold_smiles']}`) at {top_pct:.1f}%",
        ]
        if top3_cumulative is not None:
            lines.append(f"- Top 3 scaffolds cover {top3_cumulative:.1f}% of dataset")
        if has_empty:
            lines.append("- Empty scaffold present (non-ring / acyclic molecules in dataset)")

        # singleton count
        singletons = [r for r in dist_rows if int(r["count"]) == 1]
        singleton_pct = 100.0 * len(singletons) / n_unique if n_unique else 0
        if singleton_pct > 50:
            lines.append(
                f"- ⚠️ {len(singletons)}/{n_unique} scaffolds appear only once ({singleton_pct:.0f}%) "
                "→ long tail, limited generalisation per scaffold"
            )
        else:
            lines.append(f"- {len(singletons)}/{n_unique} scaffolds are singletons")
    else:
        lines.append("_`scaffold_distribution.csv` not found._")

    lines += ["", "---", ""]

    # --- scaffold_performance.csv ---
    lines += [
        "### `scaffold_performance.csv` ⭐ **KEY FILE**",
        "",
        "**Purpose:** Evaluates model prediction quality for each scaffold group.",
        "",
        "**Columns:**",
        "| Column | Description | Type | Unit |",
        "|--------|-------------|------|------|",
        "| `scaffold_smiles` | SMILES of the scaffold | string | - |",
        "| `n_samples` | Number of molecules with this scaffold | int | - |",
        "| `coverage_pct` | Percentage of dataset | float | % |",
        "| `rmse` | Root Mean Square Error | float | days |",
        "| `r2` | R² coefficient of determination | float | - |",
        "| `mean_error` | Average prediction bias | float | days |",
        "",
    ]

    if perf_rows:
        lines += [
            "**Data (top 10 by coverage):**",
            _csv_to_fenced(perf_rows),
            "",
            "**Interpretation:**",
            "",
            "| Metric | Good | Medium | Bad | Action |",
            "|--------|------|--------|-----|--------|",
            "| **R²** | >0.5 | 0.2-0.5 | <0 or negative | Don't trust predictions for bad scaffolds |",
            "| **RMSE** | <200 days | 200-500 days | >500 days | High uncertainty |",
            "| **Mean error** | ±50 days | ±50-200 days | >200 days | Systematic bias |",
            "| **n_samples** | ≥5 | 3-4 | 1-2 | Small groups unreliable |",
            "",
            f"**Findings ({chem_dir.parent.parent.name}):**",
            "",
            "```",
        ]

        for row in perf_rows:
            smiles = row["scaffold_smiles"]
            label = _scaffold_label(smiles)
            try:
                n = int(row["n_samples"])
                r2 = float(row["r2"])
                rmse = float(row["rmse"])
                mean_err = float(row["mean_error"])
                cov = float(row["coverage_pct"])
            except (ValueError, KeyError):
                continue

            r2_emoji, r2_label = _rate_r2(r2, n)
            rmse_emoji, rmse_label = _rate_rmse(rmse)

            # Overall rating: worst of R² and RMSE
            overall = (
                r2_emoji if r2_emoji == "❌" or rmse_emoji == "❌" else ("⚠️" if r2_emoji == "⚠️" or rmse_emoji == "⚠️" else "✅")
            )

            bias_note = ""
            if abs(mean_err) > 200:
                bias_dir = "over" if mean_err > 0 else "under"
                bias_note = f" | systematic {bias_dir}-prediction ({mean_err:+.0f} days)"

            lines.append(f"{overall} {r2_label} R²: scaffold={label} ({n} samples, {cov:.1f}%)")
            lines.append(f"   R²={r2:.3f}, RMSE={rmse:.0f} days, mean_error={mean_err:+.1f} days{bias_note}")
            lines.append("")

        lines += [
            "```",
            "",
            "**Actionable insights:**",
        ]

        good = [
            r
            for r in perf_rows
            if _rate_r2(float(r["r2"]), int(r["n_samples"]))[0] == "✅" and _rate_rmse(float(r["rmse"]))[0] == "✅"
        ]
        bad = [
            r
            for r in perf_rows
            if _rate_r2(float(r["r2"]), int(r["n_samples"]))[0] == "❌" or _rate_rmse(float(r["rmse"]))[0] == "❌"
        ]

        if good:
            labels = ", ".join(_scaffold_label(r["scaffold_smiles"]) for r in good)
            lines.append(f"- **Trust:** {labels} predictions")
        if bad:
            labels = ", ".join(_scaffold_label(r["scaffold_smiles"]) for r in bad)
            lines.append(f"- **Ignore:** {labels} predictions")
        med = [r for r in perf_rows if r not in good and r not in bad]
        if med:
            labels = ", ".join(_scaffold_label(r["scaffold_smiles"]) for r in med)
            lines.append(f"- **Question:** {labels} predictions")
    else:
        lines.append("_`scaffold_performance.csv` not found._")

    lines += ["", "---", ""]
    return "\n".join(lines)


def _section2_substituent(chem_dir: Path) -> str:
    rgroup_path = chem_dir / "rgroup_analysis.csv"
    lines: list[str] = [
        "## Section 2: Substituent Analysis",
        "",
        "### `rgroup_analysis.csv` (if generated)",
        "",
        "**Purpose:** Identifies R-groups (substituents) attached to dominant scaffolds and "
        "their effect on predicted activity.",
        "",
        "**Columns:**",
        "| Column | Description | Type |",
        "|--------|-------------|------|",
        "| `scaffold_smiles` | Core scaffold | string |",
        "| `rgroup_smiles` | R-group substituent | string |",
        "| `position` | Attachment position on scaffold | int |",
        "| `count` | Occurrences of this R-group | int |",
        "| `mean_activity` | Average half-life | float |",
        "| `activity_std` | Standard deviation | float |",
        "",
        "**Interpretation:**",
        "",
        "| Pattern | Chemical Meaning |",
        "|---------|------------------|",
        "| NO₂, CN → high activity | Electron-withdrawing groups increase persistence |",
        "| OH, NH₂ → low activity | Electron-donating groups decrease persistence |",
        "| Cl, Br → variable | Halogen effect depends on position/context |",
        "| High activity_std | R-group effect is context-dependent |",
        "",
    ]

    if rgroup_path.exists():
        rgroup_rows = _read_csv(rgroup_path)
        lines += [
            "**Data:**",
            _csv_to_fenced(rgroup_rows[:10]),
        ]
    else:
        lines.append(
            "**Note:** `rgroup_analysis.csv` was not generated for this run "
            "(requires multiple R-groups on the same scaffold with sufficient data)."
        )

    lines += ["", "---", ""]
    return "\n".join(lines)


def _section3_fragments(chem_dir: Path) -> str:
    freq_rows = _read_csv(chem_dir / "fragment_frequency.csv")
    outlier_rows = _read_csv(chem_dir / "outlier_fragments.csv")

    lines: list[str] = ["## Section 3: Fragment Analysis", ""]

    # --- fragment_frequency.csv ---
    lines += [
        "### `fragment_frequency.csv`",
        "",
        "**Purpose:** Lists BRICS fragments and their associated biodegradation activity.",
        "",
        "**Columns:**",
        "| Column | Description | Type | Unit |",
        "|--------|-------------|------|------|",
        "| `fragment_smiles` | SMILES of fragment (with `[*]` attachment points) | string | - |",
        "| `frequency` | Number of occurrences | int | - |",
        "| `mean_activity` | Mean half-life for molecules containing fragment | float | days |",
        "| `activity_std` | Standard deviation of activity | float | days |",
        "",
    ]

    if freq_rows:
        # Show top rows (freq >= 3 where possible, up to 10)
        usable = [r for r in freq_rows if int(r.get("frequency", 0)) >= 3]
        display_rows = usable[:10] if usable else freq_rows[:10]
        lines += [
            "**Example data (fragments with frequency ≥ 3, up to 10 shown):**",
            _csv_to_fenced(display_rows),
            "",
            "**Interpretation:**",
            "",
            "| Metric | Good | Bad | Action |",
            "|--------|------|-----|--------|",
            "| **Frequency** | ≥3 | 1-2 | Cannot generalise from singletons |",
            "| **activity_std** | <100 days | >500 days | Fragment too context-dependent |",
            "| **Trend clarity** | Chemically sensible | Contradicts literature | Check for data artifacts |",
            "",
            "**Chemical patterns to expect:**",
            "| Fragment type | Expected activity | Reason |",
            "|---------------|-------------------|--------|",
            "| Halogenated aromatics | High (persistent) | C-Cl bond stability, toxicity |",
            "| Nitro groups | High (persistent) | Electron-withdrawing, recalcitrant |",
            "| Esters, ethers | Low (labile) | Hydrolyzable bonds |",
            "| Carboxylic acids | Low-moderate | Readily metabolized |",
            "| Unsubstituted aromatics | Moderate | Ring cleavage required |",
            "",
            f"**Findings ({chem_dir.parent.parent.name}):**",
            "",
            "```",
        ]

        for row in usable[:15]:
            frag = row["fragment_smiles"]
            freq = int(row["frequency"])
            mean_act = float(row["mean_activity"])
            std_val = row.get("activity_std", "")
            try:
                std = float(std_val)
            except (ValueError, TypeError):
                std = None

            # Rating
            if std is not None and std > 500:
                emoji = "⚠️ WARNING"
            elif mean_act > 500:
                emoji = "⚠️ HIGH PERSISTENCE"
            elif mean_act < 100 and (std is None or std < 100):
                emoji = "✅ GOOD"
            else:
                emoji = "→"

            std_str = f", std={std:.1f}" if std is not None else ""
            lines.append(f"{emoji}: {frag}")
            lines.append(f"   freq={freq}, mean={mean_act:.1f} days{std_str}")
            lines.append("")

        lines.append("```")
    else:
        lines.append("_`fragment_frequency.csv` not found._")

    lines += ["", "---", ""]

    # --- outlier_fragments.csv ---
    lines += [
        "### `outlier_fragments.csv`",
        "",
        "**Purpose:** Identifies fragments associated with outlier predictions (large prediction errors).",
        "",
        "**Columns:**",
        "| Column | Description | Type |",
        "|--------|-------------|------|",
        "| `fragment_smiles` | SMILES of fragment | string |",
        "| `in_high_activity` | Present in high-activity molecules? | bool |",
        "| `in_outliers` | Present in outlier predictions? | bool |",
        "",
    ]

    if outlier_rows:
        lines += [
            "**Data:**",
            _csv_to_fenced(outlier_rows),
            "",
            "**Interpretation:**",
            "",
            "| Pattern | Meaning |",
            "|---------|---------|",
            "| `in_outliers=True, in_high_activity=False` | Fragment confuses model (not just extreme activity) |",
            "| `in_outliers=True, in_high_activity=True` | Fragment associated with extreme but predictable activity |",
            "| `in_outliers=False, in_high_activity=True` | Fragment well-handled by model |",
            "",
            f"**Findings ({chem_dir.parent.parent.name}):**",
            "",
            "```",
        ]

        for row in outlier_rows:
            frag = row["fragment_smiles"]
            in_high = row.get("in_high_activity", "").strip()
            in_out = row.get("in_outliers", "").strip()

            if in_out.lower() == "true" and in_high.lower() == "false":
                meaning = "Model doesn't know how to handle this chemistry → likely extrapolation"
            elif in_out.lower() == "true" and in_high.lower() == "true":
                meaning = "Extreme activity AND outlier → unusual but patterned chemistry"
            else:
                meaning = "Fragment well-handled by model"

            lines.append(f"Fragment: {frag}")
            lines.append(f"  in_outliers={in_out}, in_high_activity={in_high}")
            lines.append(f"  → {meaning}")
            lines.append("")

        lines += [
            "```",
            "",
            "**Actionable insights:**",
            "- Outlier fragments define the **boundary of model applicability**",
            "- Highly specific fragments (polychlorinated, fused cages) → likely extrapolation",
            "- Consider adding similar structures to training data if these are important",
        ]
    else:
        lines.append("_`outlier_fragments.csv` not found or empty._")

    lines += ["", "---", ""]
    return "\n".join(lines)


def _section4_consistency(chem_dir: Path) -> str:
    desc_map = _read_json(chem_dir / "descriptor_chemistry_map.json")
    motif_rows = _read_csv(chem_dir / "motif_stability.csv")
    ad_rows = _read_csv(chem_dir / "ad_scaffold_analysis.csv")

    lines: list[str] = ["## Section 4: Consistency Checks", ""]

    # --- descriptor_chemistry_map.json ---
    lines += [
        "### `descriptor_chemistry_map.json`",
        "",
        "**Purpose:** Maps molecular descriptors to their feature importance.",
        "",
        "**Structure:**",
        "```json",
        "{",
        '  "DESCRIPTOR_NAME": {',
        '    "importance": 0.00928,',
        '    "interpretation": "Structural fingerprint bit"',
        "  }",
        "}",
        "```",
        "",
        "**Descriptor types:**",
        "",
        "| Prefix | Type | Interpretability |",
        "|--------|------|------------------|",
        "| `fr_*` | Functional group count | High |",
        "| `Num*` | Count descriptor | High |",
        "| `PEOE_*` / `SlogP_*` | Partial charge / surface area | Medium |",
        "| `Min/MaxPartialCharge` | Partial charge extremes | Medium |",
        "| `NumH*` | H-bond counts | Medium |",
        "| `MACCS_*` | MACCS fingerprint bit | Low (structural pattern) |",
        "",
    ]

    if desc_map:
        # Sort by importance descending
        sorted_desc = sorted(desc_map.items(), key=lambda x: x[1].get("importance", 0), reverse=True)
        top5 = sorted_desc[:5]

        # Count feature types
        n_fr = sum(1 for k, _ in sorted_desc if k.startswith("fr_"))
        n_num = sum(1 for k, _ in sorted_desc if k.startswith("Num"))
        n_maccs = sum(1 for k, _ in sorted_desc if k.startswith("MACCS"))
        n_physchem = len(sorted_desc) - n_fr - n_num - n_maccs

        lines += [
            "**Example data (all features, sorted by importance):**",
            "```json",
            json.dumps(desc_map, indent=2),
            "```",
            "",
            f"**Findings ({chem_dir.parent.parent.name}):**",
            "",
            "```",
        ]

        # Top feature assessment
        top_name, top_info = top5[0]
        top_type = top_info.get("interpretation", "")
        if "fingerprint" in top_type.lower():
            lines.append(f"⚠️ WARNING: {top_name} is top feature")
            lines.append(f"   → Structural fingerprint bit — unclear chemical meaning")
        else:
            lines.append(f"✅ GOOD: {top_name} is top feature")
            lines.append(f"   → {top_type}")
        lines.append("")

        # Mix assessment
        if n_fr > 0 or n_num > 0:
            lines.append(f"✅ GOOD: Functional group / count descriptors present ({n_fr} fr_*, {n_num} Num*)")
        if n_maccs > 0 and n_fr == 0 and n_num == 0:
            lines.append(f"❌ BAD: Only MACCS fingerprint bits ({n_maccs}) — no interpretable functional group features")
        elif n_maccs > 0:
            lines.append(f"⚠️ Note: {n_maccs} MACCS fingerprint bit(s) in feature set")

        if n_physchem > 0:
            lines.append(f"✅ GOOD: Physicochemical/charge descriptors present ({n_physchem} features)")

        lines += ["```", ""]
    else:
        lines.append("_`descriptor_chemistry_map.json` not found._")

    lines += ["---", ""]

    # --- motif_stability.csv ---
    lines += [
        "### `motif_stability.csv`",
        "",
        "**Purpose:** Assesses prediction consistency within scaffold groups.",
        "",
        "**Columns:**",
        "| Column | Description | Type | Unit |",
        "|--------|-------------|------|------|",
        "| `scaffold_smiles` | Core scaffold | string | - |",
        "| `n_samples` | Number of samples | int | - |",
        "| `prediction_std` | Standard deviation of predictions | float | days |",
        "| `prediction_variance` | Variance of predictions | float | days² |",
        "| `mean_absolute_error` | MAE for this scaffold | float | days |",
        "| `stable` | Stability flag (std < threshold) | bool | - |",
        "",
    ]

    if motif_rows:
        lines += [
            "**Data:**",
            _csv_to_fenced(motif_rows),
            "",
            "**Interpretation:**",
            "",
            "| Metric | Good | Medium | Bad |",
            "|--------|------|--------|-----|",
            "| **stable** | True | - | False |",
            "| **prediction_std** | <200 days | 200-500 days | >500 days |",
            "| **mean_absolute_error** | <100 days | 100-300 days | >300 days |",
            "",
            f"**Findings ({chem_dir.parent.parent.name}):**",
            "",
            "```",
        ]

        for row in motif_rows:
            smiles = row["scaffold_smiles"]
            label = _scaffold_label(smiles)
            try:
                n = int(row["n_samples"])
                pred_std = float(row["prediction_std"])
                mae = float(row["mean_absolute_error"])
                stable = row["stable"]
            except (ValueError, KeyError):
                continue

            emoji, label_rating = _rate_stability(stable, pred_std)
            lines.append(f"{emoji} {label_rating}: {label} ({n} samples)")
            lines.append(f"   std={pred_std:.1f} days, MAE={mae:.1f} days, stable={stable}")
            lines.append("")

        lines.append("```")
        lines += [
            "",
            "**Key insight:** Compare with `scaffold_performance.csv`:",
            "- **Both good:** Reliable scaffold",
            "- **Good R², bad stability:** Model fits but predictions are sensitive",
            "- **Both bad:** Avoid this chemistry",
        ]
    else:
        lines.append("_`motif_stability.csv` not found._")

    lines += ["", "---", ""]

    # --- ad_scaffold_analysis.csv ---
    lines += [
        "### `ad_scaffold_analysis.csv` ⭐ **CRITICAL FILE**",
        "",
        "**Purpose:** Shows which scaffolds fall outside the model's Applicability Domain (AD).",
        "",
        "**Columns:**",
        "| Column | Description | Type | Unit |",
        "|--------|-------------|------|------|",
        "| `scaffold_smiles` | Core scaffold | string | - |",
        "| `n_samples` | Total samples | int | - |",
        "| `n_outside_ad` | Count outside AD | int | - |",
        "| `pct_outside_ad` | Percentage outside AD | float | % |",
        "",
    ]

    if ad_rows:
        # Sort by n_samples descending for display
        ad_sorted = sorted(ad_rows, key=lambda r: int(r.get("n_samples", 0)), reverse=True)
        lines += [
            "**Data (sorted by n_samples):**",
            _csv_to_fenced(ad_sorted[:15]),
            "",
            "**Interpretation:**",
            "",
            "| Metric | Good | Medium | Bad | Action |",
            "|--------|------|--------|-----|--------|",
            "| **pct_outside_ad** | <20% | 20-50% | >50% | Don't trust predictions if bad |",
            "| **Common scaffolds outside AD** | Rare | - | Major scaffold >50% | Retrain with more diverse data |",
            "",
            f"**Findings ({chem_dir.parent.parent.name}):**",
            "",
            "```",
        ]

        # Only report scaffolds with n_samples >= 2 to avoid clutter
        major = [r for r in ad_sorted if int(r.get("n_samples", 0)) >= 2]
        for row in major[:12]:
            smiles = row["scaffold_smiles"]
            label = _scaffold_label(smiles)
            try:
                n = int(row["n_samples"])
                n_out = int(row["n_outside_ad"])
                pct = float(row["pct_outside_ad"])
            except (ValueError, KeyError):
                continue

            emoji, ad_label = _rate_ad(pct)
            lines.append(f"{emoji} {ad_label}: {label} ({n} samples)")
            lines.append(f"   {n_out}/{n} outside AD ({pct:.1f}%)")
            lines.append("")

        lines += [
            "```",
            "",
            "**Critical insight:**",
            "```",
            "High R² + High % outside AD = Model is LUCKY, not RELIABLE",
            "",
            "Example pattern:",
            "  R²>0.8 but 100% outside AD",
            "  → Model may have memorized small sample, not learned generalizable patterns",
            "  → Do NOT trust predictions for new molecules with this scaffold",
            "```",
        ]
    else:
        lines.append("_`ad_scaffold_analysis.csv` not found._")

    lines += ["", "---", ""]
    return "\n".join(lines)


def _section_visualizations() -> str:
    lines = [
        "## Visualizations",
        "",
        "### `scaffold_bar_plot.png`",
        "",
        "**What it shows:** Bar chart of top 10 scaffolds by count",
        "",
        "**Use:** Quick visual of dataset composition",
        "",
        "### `fragment_activity_plot.png`",
        "",
        "**What it shows:** Scatter plot of fragment frequency vs. mean activity",
        "",
        "**Use:** Identify fragments that correlate with high/low persistence",
        "",
        "### `substituent_trends.png` (if generated)",
        "",
        "**What it shows:** R-group effect on activity for dominant scaffolds",
        "",
        "**Use:** Visualize substituent electronic/steric effects",
        "",
        "---",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------


def build_interpretation(chem_dir: Path) -> str:
    parts = [
        _section_header(chem_dir),
        _section1_scaffold(chem_dir),
        _section2_substituent(chem_dir),
        _section3_fragments(chem_dir),
        _section4_consistency(chem_dir),
        _section_visualizations(),
    ]
    return "\n".join(parts)


def main() -> None:
    chem_dirs = find_chemistry_dirs(LOG_DIR)
    if not chem_dirs:
        print("No chemistry_analysis/ directories found.")
        return

    written = 0
    skipped = 0
    for chem_dir in chem_dirs:
        output_path = chem_dir / OUTPUT_FILENAME
        if output_path.exists():
            print(f"  SKIP (already exists): {output_path.relative_to(LOG_DIR)}")
            skipped += 1
            continue

        md = build_interpretation(chem_dir)
        output_path.write_text(md, encoding="utf-8")
        print(f"  WROTE: {output_path.relative_to(LOG_DIR)}")
        written += 1

    print(f"\nDone. Written: {written}, Skipped: {skipped}")


if __name__ == "__main__":
    main()
