#!/usr/bin/env python3
"""
compare_abundance_files.py
==========================
Compare two transcript abundance TSV files from two experiment folders.

Usage (simple — just give experiment folder names):
    python compare_abundance_files.py exprmnt_2026_03_14__10_43_29 exprmnt_2026_03_15__10_42_32

    RESULTS_BASE and SAMPLE in the CONFIG section below tell the script where
    to find the abundance.tsv for each experiment.

Alternatively, override paths explicitly:
    python compare_abundance_files.py exprmnt_A exprmnt_B \
        --results-base /path/to/results \
        --sample ds_52_furtherDownsampled

What it does
------------
1. Loads abundance.tsv from each experiment folder.
2. Detects transcript-ID and TPM columns automatically.
3. Prints per-file stats: total transcripts, non-zero transcripts.
4. Reports overlap: how many transcripts are shared (non-zero in both).
5. Computes similarity metrics (Pearson, Spearman, MAE, RMSE, Bray-Curtis, JS).
6. Saves a timestamped comparison report into EACH experiment folder:
       {RESULTS_BASE}/{exp1}/comparison_vs_{exp2}_{timestamp}.txt
       {RESULTS_BASE}/{exp2}/comparison_vs_{exp1}_{timestamp}.txt

Inputs:
    exp1         : experiment folder name (e.g. exprmnt_2026_03_14__10_43_29)
    exp2         : experiment folder name (e.g. exprmnt_2026_03_15__10_42_32)
    RESULTS_BASE : base directory containing experiment folders (see CONFIG)
    SAMPLE       : subfolder within each experiment (see CONFIG)

Outputs:
    - Report printed to stdout (with timestamp header)
    - Report saved to each experiment folder as comparison_vs_{other}_{timestamp}.txt
"""

from __future__ import annotations

import argparse
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd


# ============================================================
# CONFIG — edit these before running; do not edit below
# ============================================================
RESULTS_BASE   = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"
SAMPLE         = "ds_52_furtherDownsampled"   # subfolder inside each experiment
ABUNDANCE_FILE = "abundance.tsv"              # filename within the sample folder
# Labels are auto-detected from running.log (see detect_pipeline_label()).
# Add entries here if you add new pipeline scripts.
SCRIPT_LABELS = {
    "run_joli_kallisto.sh": "JOLI-Kallisto",
    "run_lr_kallisto.sh":   "lr-kallisto",
}
# ============================================================


def detect_pipeline_label(exp_dir: Path) -> str:
    """
    Read the first line of running.log to identify which pipeline produced this
    experiment, then map the script filename to a human-readable label.

    The first line is written by the pipeline script as:
        Script: /full/path/to/run_joli_kallisto.sh

    Args:
        exp_dir : Path -- experiment root folder (contains running.log).

    Returns:
        str -- human-readable label (e.g. "JOLI-Kallisto" or "lr-kallisto"),
               or the raw script basename if not in SCRIPT_LABELS,
               or the experiment folder name if running.log is missing / has no
               Script: line (graceful fallback for old experiments).
    """
    log_path = exp_dir / "running.log"
    if not log_path.exists():
        print(f"  [label] No running.log in {exp_dir.name} — using folder name as label")
        return exp_dir.name

    with open(log_path) as fh:
        first_line = fh.readline().strip()   # "Script: /path/to/script.sh"

    if not first_line.startswith("Script:"):
        print(f"  [label] running.log in {exp_dir.name} has no 'Script:' first line "
              f"(old format) — using folder name as label")
        return exp_dir.name

    script_path = first_line[len("Script:"):].strip()
    script_name = Path(script_path).name   # e.g. "run_joli_kallisto.sh"

    label = SCRIPT_LABELS.get(script_name, script_name)
    print(f"  [label] {exp_dir.name} → '{label}' (from {script_name})")
    return label


DEFAULT_ID_CANDIDATES = [
    "transcript_id", "target_id", "Name", "name",
    "transcript", "target", "id", "ID",
]

DEFAULT_VALUE_CANDIDATES = [
    "tpm", "TPM", "est_counts", "NumReads",
    "count", "counts", "abundance", "expr", "expression",
]


# ---------- utility functions ----------

def detect_column(columns: Iterable[str], preferred_names: list, kind: str) -> str:
    cols   = list(columns)
    colset = set(cols)
    for name in preferred_names:
        if name in colset:
            return name
    lowered = {c.lower(): c for c in cols}
    for name in preferred_names:
        if name.lower() in lowered:
            return lowered[name.lower()]
    raise ValueError(
        f"Could not automatically detect {kind} column. Available columns: {cols}"
    )


def read_abundance_table(
    filepath: Path,
    id_col: Optional[str],
    value_col: Optional[str],
) -> Tuple[pd.DataFrame, str, str]:
    """
    Read an abundance TSV, returning a two-column DataFrame (transcript_id, value).

    Args:
        filepath  : Path  -- path to the abundance TSV file.
        id_col    : str   -- override for ID column name (None = auto-detect).
        value_col : str   -- override for value column name (None = auto-detect).

    Returns:
        Tuple of (DataFrame with columns [transcript_id, value],
                  detected ID column name,
                  detected value column name).
    """
    df = pd.read_csv(filepath, sep="\t", comment="#")
    if df.empty:
        raise ValueError(f"File is empty: {filepath}")

    chosen_id    = id_col    or detect_column(df.columns, DEFAULT_ID_CANDIDATES, "ID")
    chosen_value = value_col or detect_column(df.columns, DEFAULT_VALUE_CANDIDATES, "abundance")

    out = df[[chosen_id, chosen_value]].copy()
    out.columns = ["transcript_id", "value"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["value"]         = pd.to_numeric(out["value"], errors="coerce")
    out = out.dropna(subset=["transcript_id", "value"])
    out = out.groupby("transcript_id", as_index=False)["value"].sum()
    return out, chosen_id, chosen_value


# ---------- metric functions ----------

def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    xr = pd.Series(x).rank(method="average").to_numpy()
    yr = pd.Series(y).rank(method="average").to_numpy()
    return pearson_corr(xr, yr)


def rmse(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.sqrt(np.mean((x - y) ** 2)))


def mae(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean(np.abs(x - y)))


def bray_curtis_distance(x: np.ndarray, y: np.ndarray) -> float:
    denom = np.sum(np.abs(x + y))
    if denom == 0:
        return 0.0
    return float(np.sum(np.abs(x - y)) / denom)


def normalize_nonnegative(x: np.ndarray) -> np.ndarray:
    x = np.clip(np.asarray(x, dtype=float), a_min=0.0, a_max=None)
    s = x.sum()
    if s == 0:
        return np.full_like(x, 1.0 / len(x)) if len(x) > 0 else x
    return x / s


def kl_divergence(p: np.ndarray, q: np.ndarray, eps: float = 1e-12) -> float:
    p = np.clip(p, eps, None)
    q = np.clip(q, eps, None)
    return float(np.sum(p * np.log2(p / q)))


def jensen_shannon_distance(x: np.ndarray, y: np.ndarray) -> float:
    p = normalize_nonnegative(x)
    q = normalize_nonnegative(y)
    m = 0.5 * (p + q)
    js_div = 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)
    return float(np.sqrt(max(js_div, 0.0)))


def summarize_pair(x: np.ndarray, y: np.ndarray, prefix: str) -> dict:
    return {
        f"{prefix}_pearson":     pearson_corr(x, y),
        f"{prefix}_spearman":    spearman_corr(x, y),
        f"{prefix}_mae":         mae(x, y),
        f"{prefix}_rmse":        rmse(x, y),
        f"{prefix}_bray_curtis": bray_curtis_distance(x, y),
        f"{prefix}_js_distance": jensen_shannon_distance(x, y),
    }


def fmt(x: float) -> str:
    return "nan" if pd.isna(x) else f"{x:.6f}"


# ---------- main comparison ----------

def compare_tables(df1: pd.DataFrame, df2: pd.DataFrame) -> tuple:
    """
    Merge two abundance tables and compute all similarity metrics.

    Args:
        df1 : pd.DataFrame -- columns [transcript_id, value] for file 1.
        df2 : pd.DataFrame -- columns [transcript_id, value] for file 2.

    Returns:
        Tuple of (merged DataFrame, results dict).
    """
    merged = df1.merge(df2, on="transcript_id", how="outer", suffixes=("_1", "_2"))
    merged["present_in_file1"] = merged["value_1"].notna()
    merged["present_in_file2"] = merged["value_2"].notna()
    merged["shared"]           = merged["present_in_file1"] & merged["present_in_file2"]
    merged["value_1_filled"]   = merged["value_1"].fillna(0.0)
    merged["value_2_filled"]   = merged["value_2"].fillna(0.0)
    merged["abs_diff"]         = np.abs(merged["value_1_filled"] - merged["value_2_filled"])
    merged["log1p_value_1"]    = np.log1p(merged["value_1_filled"])
    merged["log1p_value_2"]    = np.log1p(merged["value_2_filled"])
    merged["log1p_abs_diff"]   = np.abs(merged["log1p_value_1"] - merged["log1p_value_2"])

    shared   = merged.loc[merged["shared"]].copy()
    x_shared = shared["value_1"].to_numpy(dtype=float)
    y_shared = shared["value_2"].to_numpy(dtype=float)
    x_all    = merged["value_1_filled"].to_numpy(dtype=float)
    y_all    = merged["value_2_filled"].to_numpy(dtype=float)

    # Non-zero counts per file and overlap
    nz1 = int((df1["value"] > 0).sum())
    nz2 = int((df2["value"] > 0).sum())
    # Non-zero in both: rows where both values are strictly > 0
    nz_both = int(((merged["value_1_filled"] > 0) & (merged["value_2_filled"] > 0)).sum())

    results = {
        # --- per-file totals ---
        "n_total_file1":          int(df1.shape[0]),
        "n_nonzero_file1":        nz1,
        "n_total_file2":          int(df2.shape[0]),
        "n_nonzero_file2":        nz2,
        # --- overlap ---
        "n_union_all":            int(merged.shape[0]),
        "n_shared_present":       int(shared.shape[0]),     # present in both files (may include zeros)
        "n_nonzero_both":         nz_both,                  # non-zero in both
        "n_only_file1":           int((merged["present_in_file1"] & ~merged["present_in_file2"]).sum()),
        "n_only_file2":           int((~merged["present_in_file1"] & merged["present_in_file2"]).sum()),
        "shared_fraction_file1":  float(shared.shape[0] / df1.shape[0]) if df1.shape[0] else float("nan"),
        "shared_fraction_file2":  float(shared.shape[0] / df2.shape[0]) if df2.shape[0] else float("nan"),
        "jaccard_ids":            float(shared.shape[0] / merged.shape[0]) if merged.shape[0] else float("nan"),
    }

    # Metrics on shared transcripts
    if len(x_shared) > 0:
        results.update(summarize_pair(x_shared,              y_shared,              prefix="shared_raw"))
        results.update(summarize_pair(np.log1p(x_shared),    np.log1p(y_shared),    prefix="shared_log1p"))
    else:
        for key in ["shared_raw_pearson", "shared_raw_spearman", "shared_raw_mae",
                    "shared_raw_rmse", "shared_raw_bray_curtis", "shared_raw_js_distance",
                    "shared_log1p_pearson", "shared_log1p_spearman", "shared_log1p_mae",
                    "shared_log1p_rmse", "shared_log1p_bray_curtis", "shared_log1p_js_distance"]:
            results[key] = float("nan")

    # Metrics on union (missing → 0)
    results.update(summarize_pair(x_all,              y_all,              prefix="union_zero_filled_raw"))
    results.update(summarize_pair(np.log1p(x_all),    np.log1p(y_all),    prefix="union_zero_filled_log1p"))

    return merged, results


# ---------- report generation ----------

def build_report(
    results:     dict,
    exp1:        str,
    exp2:        str,
    label1:      str,
    label2:      str,
    file1:       Path,
    file2:       Path,
    id_col1:     str,
    id_col2:     str,
    value_col1:  str,
    value_col2:  str,
    timestamp:   str,
) -> str:
    """
    Build the full text report string.

    Args:
        results      : dict   -- metric results from compare_tables().
        exp1/exp2    : str    -- experiment folder names.
        label1/label2: str    -- short human-readable labels (e.g. "LK", "JK").
        file1/file2  : Path   -- full paths to abundance files.
        id_col1/2    : str    -- detected ID column names.
        value_col1/2 : str    -- detected value column names.
        timestamp    : str    -- run timestamp string.

    Returns:
        str -- formatted report text.
    """
    # Display name = "LABEL (exprmnt_...)" when a label is provided, else just the folder name
    name1 = f"{label1} ({exp1})" if label1 != exp1 else exp1
    name2 = f"{label2} ({exp2})" if label2 != exp2 else exp2

    lines = []
    lines.append("=" * 80)
    lines.append(f"Abundance comparison report")
    lines.append(f"Timestamp : {timestamp}")
    lines.append(f"Comparing : {name1}")
    lines.append(f"      vs  : {name2}")
    lines.append("=" * 80)
    lines.append("")
    lines.append(f"File 1 — {name1}")
    lines.append(f"  Path        : {file1}")
    lines.append(f"  ID column   : {id_col1}")
    lines.append(f"  Value column: {value_col1}")
    lines.append("")
    lines.append(f"File 2 — {name2}")
    lines.append(f"  Path        : {file2}")
    lines.append(f"  ID column   : {id_col2}")
    lines.append(f"  Value column: {value_col2}")
    lines.append("")

    # --- Per-file counts (use labels as column headers, no truncation) ---
    col_w = max(14, len(label1), len(label2))   # column width fits the longest label
    lines.append("Transcript counts")
    lines.append("-" * 50)
    lines.append(f"  {'':30s}  {label1:>{col_w}s}  {label2:>{col_w}s}")
    lines.append(f"  {'Total transcripts':30s}  {results['n_total_file1']:>{col_w}d}  {results['n_total_file2']:>{col_w}d}")
    lines.append(f"  {'Non-zero transcripts':30s}  {results['n_nonzero_file1']:>{col_w}d}  {results['n_nonzero_file2']:>{col_w}d}")
    lines.append(f"  {'Zero (suppressed)':30s}  {results['n_total_file1']-results['n_nonzero_file1']:>{col_w}d}  {results['n_total_file2']-results['n_nonzero_file2']:>{col_w}d}")
    lines.append("")

    # --- Overlap ---
    lines.append("Overlap")
    lines.append("-" * 50)
    lines.append(f"  Union (all transcripts seen in either file) : {results['n_union_all']}")
    lines.append(f"  Shared (present in both files)              : {results['n_shared_present']}")
    lines.append(f"  Non-zero in both files                      : {results['n_nonzero_both']}")
    lines.append(f"  Only in {label1}                            : {results['n_only_file1']}")
    lines.append(f"  Only in {label2}                            : {results['n_only_file2']}")
    lines.append(f"  Shared fraction of {label1}                 : {fmt(results['shared_fraction_file1'])}")
    lines.append(f"  Shared fraction of {label2}                 : {fmt(results['shared_fraction_file2'])}")
    lines.append(f"  Jaccard (shared / union)                    : {fmt(results['jaccard_ids'])}")
    lines.append("")

    # Each section: (title, prefix_to_strip, [metric_keys])
    sections = [
        ("Shared transcripts — raw abundance",    "shared_raw_",
         ["shared_raw_pearson", "shared_raw_spearman", "shared_raw_mae",
          "shared_raw_rmse", "shared_raw_bray_curtis", "shared_raw_js_distance"]),
        ("Shared transcripts — log1p abundance",  "shared_log1p_",
         ["shared_log1p_pearson", "shared_log1p_spearman", "shared_log1p_mae",
          "shared_log1p_rmse", "shared_log1p_bray_curtis", "shared_log1p_js_distance"]),
        ("Union (missing→0) — raw abundance",     "union_zero_filled_raw_",
         ["union_zero_filled_raw_pearson", "union_zero_filled_raw_spearman",
          "union_zero_filled_raw_mae", "union_zero_filled_raw_rmse",
          "union_zero_filled_raw_bray_curtis", "union_zero_filled_raw_js_distance"]),
        ("Union (missing→0) — log1p abundance",   "union_zero_filled_log1p_",
         ["union_zero_filled_log1p_pearson", "union_zero_filled_log1p_spearman",
          "union_zero_filled_log1p_mae", "union_zero_filled_log1p_rmse",
          "union_zero_filled_log1p_bray_curtis", "union_zero_filled_log1p_js_distance"]),
    ]

    for title, prefix, keys in sections:
        lines.append(title)
        lines.append("-" * len(title))
        for key in keys:
            v     = results.get(key, float("nan"))
            label = key[len(prefix):]   # strip section prefix, e.g. "bray_curtis" not "curtis"
            lines.append(f"  {label:30s}: {fmt(v)}")
        lines.append("")

    lines.append("Interpretation guide")
    lines.append("-" * 20)
    lines.append("  Pearson/Spearman → 1 : more similar ranking / scale")
    lines.append("  MAE/RMSE         → 0 : smaller absolute differences")
    lines.append("  Bray-Curtis      → 0 : more similar composition")
    lines.append("  JS distance      → 0 : more similar normalized distributions")

    return "\n".join(lines)


# ---------- CLI ----------

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        argparse.Namespace -- parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Compare two abundance TSV files by experiment name.\n"
            "Example: python compare_abundance_files.py "
            "exprmnt_2026_03_14__10_43_29 exprmnt_2026_03_15__10_42_32"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "exp1",
        help="First experiment folder name (e.g. exprmnt_2026_03_14__10_43_29).",
    )
    parser.add_argument(
        "exp2",
        help="Second experiment folder name (e.g. exprmnt_2026_03_15__10_42_32).",
    )
    parser.add_argument(
        "--results-base", default=RESULTS_BASE,
        help=f"Base directory containing experiment folders. Default: {RESULTS_BASE}",
    )
    parser.add_argument(
        "--sample", default=SAMPLE,
        help=f"Sample subfolder within each experiment. Default: {SAMPLE}",
    )
    parser.add_argument(
        "--abundance-file", default=ABUNDANCE_FILE,
        help=f"Abundance filename within the sample folder. Default: {ABUNDANCE_FILE}",
    )
    parser.add_argument(
        "--label1", default=None,
        help="Human-readable name for exp1 (e.g. 'JOLI-Kallisto'). "
             "Auto-detected from running.log if omitted.",
    )
    parser.add_argument(
        "--label2", default=None,
        help="Human-readable name for exp2 (e.g. 'lr-kallisto'). "
             "Auto-detected from running.log if omitted.",
    )
    parser.add_argument("--id-col1",    default=None, help="Override ID column for file 1.")
    parser.add_argument("--id-col2",    default=None, help="Override ID column for file 2.")
    parser.add_argument("--value-col1", default=None, help="Override value column for file 1.")
    parser.add_argument("--value-col2", default=None, help="Override value column for file 2.")
    parser.add_argument(
        "--top-n", type=int, default=50,
        help="Number of most-different transcripts to write out. Default: 50.",
    )
    return parser.parse_args()


def main() -> None:
    """
    Main entry point: load files, compute metrics, print and save report.
    """
    args = parse_args()

    timestamp = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")

    # --- Resolve file paths ---
    base   = Path(args.results_base)
    file1  = base / args.exp1 / args.sample / args.abundance_file
    file2  = base / args.exp2 / args.sample / args.abundance_file

    for f, exp in [(file1, args.exp1), (file2, args.exp2)]:
        if not f.exists():
            raise FileNotFoundError(
                f"Abundance file not found for {exp}:\n  {f}\n"
                f"Check --results-base, --sample, --abundance-file."
            )

    # Auto-detect pipeline labels from running.log; --label1/--label2 override
    label1 = args.label1 if args.label1 else detect_pipeline_label(base / args.exp1)
    label2 = args.label2 if args.label2 else detect_pipeline_label(base / args.exp2)

    # --- Load ---
    df1, detected_id1, detected_val1 = read_abundance_table(file1, args.id_col1, args.value_col1)
    df2, detected_id2, detected_val2 = read_abundance_table(file2, args.id_col2, args.value_col2)

    # --- Compare ---
    merged, results = compare_tables(df1, df2)

    # --- Build report ---
    report = build_report(
        results=results,
        exp1=args.exp1,
        exp2=args.exp2,
        label1=label1,
        label2=label2,
        file1=file1,
        file2=file2,
        id_col1=detected_id1,
        id_col2=detected_id2,
        value_col1=detected_val1,
        value_col2=detected_val2,
        timestamp=timestamp,
    )

    # --- Print to stdout ---
    print(report)

    # --- Save report into each experiment folder ---
    for exp, other in [(args.exp1, args.exp2), (args.exp2, args.exp1)]:
        out_path = base / exp / f"comparison_vs_{other}__{timestamp}.txt"
        out_path.write_text(report)
        print(f"\nReport saved → {out_path}")

    # --- Save merged table and top-diff TSVs into exp1 folder (for detailed inspection) ---
    detail_dir = base / args.exp1 / f"comparison_detail_vs_{args.exp2}__{timestamp}"
    detail_dir.mkdir(parents=True, exist_ok=True)

    merged_sorted = merged.sort_values("abs_diff", ascending=False)
    merged.to_csv(detail_dir / "merged_abundance.tsv", sep="\t", index=False)
    merged_sorted.head(args.top_n).to_csv(detail_dir / "top_differences_raw.tsv", sep="\t", index=False)
    merged.sort_values("log1p_abs_diff", ascending=False).head(args.top_n).to_csv(
        detail_dir / "top_differences_log1p.tsv", sep="\t", index=False
    )
    pd.DataFrame([results]).to_csv(detail_dir / "summary_metrics.tsv", sep="\t", index=False)
    import json as _json
    (detail_dir / "summary_metrics.json").write_text(_json.dumps(results, indent=2))

    print(f"Detail tables → {detail_dir}/")


if __name__ == "__main__":
    main()
