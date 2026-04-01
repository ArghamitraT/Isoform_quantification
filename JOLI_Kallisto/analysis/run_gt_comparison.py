"""
run_gt_comparison.py
====================
Compare abundance.tsv outputs from one or more experiment folders against
ground-truth TPM files (simulation ground truth).

Auto-discovers samples inside the experiment folder. For each sample that has
a matching entry in SAMPLE_GT_MAP, loads both files and computes Spearman,
Pearson, MAE, RMSE, Bray-Curtis, and JS-distance metrics.

Three metric sets are reported per sample:
  1. All transcripts       — outer join of predicted and GT (includes TN pairs)
  2. Active universe       — GT transcripts ∪ non-zero predicted transcripts
                             (excludes true-negative (0,0) pairs; fair apples-to-
                             apples comparison between methods with different index
                             sizes, e.g. LK vs JK MS)
  3. GT non-zero only      — only transcripts present in ground truth

A contingency table (TP / FP / FN / TN counts) is printed before the metrics.

Usage:
    # Compare one experiment against ground truth:
    python analysis/run_gt_comparison.py exprmnt_2026_03_28__00_29_23

    # Compare multiple experiments (prints one results block per experiment):
    python analysis/run_gt_comparison.py exprmnt_A exprmnt_B exprmnt_C

    # Override results base:
    python analysis/run_gt_comparison.py exprmnt_A --results-base /path/to/results

Inputs:
    - One or more experiment folder names (CLI positional args).
    - SAMPLE_GT_MAP (CONFIG): maps sample subfolder name → ground truth TSV path.
    - Ground truth TSV format: CSV with columns (index, transcript_name, tpm).
    - Abundance TSV format: TSV with columns (transcript_id, tpm).

Outputs:
    - Per-sample metrics printed to stdout.
    - Summary saved to: {RESULTS_BASE}/{exp}/gt_comparison_{timestamp}.txt
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd


# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================

RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# Map sample subfolder name → absolute path to ground truth TSV.
# Ground truth TSV format: CSV with columns (unnamed_index, transcript_name, tpm)
SAMPLE_GT_MAP = {
    # Long-read (PacBio) simulated samples
    "sim1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv",
    "sim2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv",
    # Short-read (Illumina) simulated samples — same expression profile as sim1/sim2,
    # different read technology. Map to the same ground truth files.
    "sim_sr_s1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv",
    "sim_sr_s2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv",
}

ABUNDANCE_FILE = "abundance.tsv"

# ============================================================
# END CONFIG
# ============================================================


TIMESTAMP = time.strftime("%Y_%m_%d__%H_%M_%S")


# ============================================================
# I/O helpers
# ============================================================

def read_abundance(path: Path) -> pd.DataFrame:
    """
    Read an abundance TSV (tab-separated, columns: transcript_id, tpm).

    Args:
        path : Path -- path to abundance.tsv.

    Returns:
        pd.DataFrame with columns [transcript_id, tpm].
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    # Auto-detect ID column
    id_candidates = ["transcript_id", "target_id", "Name", "transcript"]
    id_col = next((c for c in id_candidates if c in df.columns), df.columns[0])
    # Auto-detect value column
    val_candidates = ["tpm", "TPM", "est_counts", "NumReads"]
    val_col = next((c for c in val_candidates if c in df.columns), df.columns[1])
    out = df[[id_col, val_col]].copy()
    out.columns = ["transcript_id", "tpm"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["tpm"] = pd.to_numeric(out["tpm"], errors="coerce").fillna(0.0)
    return out


def read_ground_truth(path: Path) -> pd.DataFrame:
    """
    Read a ground truth CSV (comma-separated, columns: unnamed_index, transcript_name, tpm).

    Args:
        path : Path -- path to ground truth TSV/CSV.

    Returns:
        pd.DataFrame with columns [transcript_id, tpm].
    """
    df = pd.read_csv(path, index_col=0)   # first column is unnamed integer index
    # Rename to standard names
    id_candidates  = ["transcript_name", "transcript_id", "target_id", "Name"]
    val_candidates = ["tpm", "TPM", "est_counts"]
    id_col  = next((c for c in id_candidates  if c in df.columns), df.columns[0])
    val_col = next((c for c in val_candidates if c in df.columns), df.columns[1])
    out = df[[id_col, val_col]].copy()
    out.columns = ["transcript_id", "tpm"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["tpm"] = pd.to_numeric(out["tpm"], errors="coerce").fillna(0.0)
    return out


# ============================================================
# Metric functions (same as compare_abundance_files.py)
# ============================================================

def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    xr = pd.Series(x).rank(method="average").to_numpy()
    yr = pd.Series(y).rank(method="average").to_numpy()
    return pearson_corr(xr, yr)


def bray_curtis(x: np.ndarray, y: np.ndarray) -> float:
    denom = np.sum(np.abs(x + y))
    return float(np.sum(np.abs(x - y)) / denom) if denom > 0 else 0.0


def mae(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean(np.abs(x - y)))


def rmse(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.sqrt(np.mean((x - y) ** 2)))


def fmt(x: float) -> str:
    return "nan" if pd.isna(x) else f"{x:.6f}"


def compute_metrics(x: np.ndarray, y: np.ndarray) -> dict:
    """
    Compute all scalar metrics between two aligned TPM arrays.

    Args:
        x : np.ndarray -- predicted TPM values.
        y : np.ndarray -- ground truth TPM values.

    Returns:
        dict with keys: spearman, pearson, mae, rmse, bray_curtis.
    """
    return {
        "spearman":    spearman_corr(x, y),
        "pearson":     pearson_corr(x, y),
        "mae":         mae(x, y),
        "rmse":        rmse(x, y),
        "bray_curtis": bray_curtis(x, y),
    }


# ============================================================
# Per-sample comparison
# ============================================================

def compare_sample(abund_path: Path, gt_path: Path, sample: str) -> dict:
    """
    Compare one sample's abundance.tsv against its ground truth.

    Computes three metric sets:
      1. all        — outer join of all predicted and GT transcripts (includes TN pairs)
      2. active     — GT transcripts ∪ non-zero predicted transcripts only
                      (excludes true-negative (0,0) pairs for fair cross-method comparison)
      3. gt_nonzero — only transcripts with non-zero ground truth TPM

    Also computes a contingency table:
      TP = pred > 0 AND gt > 0
      FP = pred > 0 AND gt = 0  (false positive — predicted non-zero, not in GT)
      FN = pred = 0 AND gt > 0  (false negative — in GT, predicted zero)
      TN = pred = 0 AND gt = 0  (true negative  — zero in both)

    Args:
        abund_path : Path -- path to abundance.tsv.
        gt_path    : Path -- path to ground truth CSV.
        sample     : str  -- sample name (for printing).

    Returns:
        dict with all computed metrics and contingency counts.
    """
    abund = read_abundance(abund_path)
    gt    = read_ground_truth(gt_path)

    # ----------------------------------------------------------
    # Universe 1: ALL transcripts (outer join, fill NA with 0)
    # Includes true-negative (0,0) pairs — inflated for methods
    # that output large reference files (e.g. JK MS with 218K).
    # ----------------------------------------------------------
    merged_all = abund.merge(gt, on="transcript_id", how="outer",
                             suffixes=("_pred", "_gt"))
    merged_all["tpm_pred"] = merged_all["tpm_pred"].fillna(0.0)
    merged_all["tpm_gt"]   = merged_all["tpm_gt"].fillna(0.0)

    pred_nz = merged_all["tpm_pred"] > 0
    gt_nz   = merged_all["tpm_gt"]   > 0

    # Contingency table counts
    n_tp = int(( pred_nz &  gt_nz).sum())   # true positive
    n_fp = int(( pred_nz & ~gt_nz).sum())   # false positive
    n_fn = int((~pred_nz &  gt_nz).sum())   # false negative
    n_tn = int((~pred_nz & ~gt_nz).sum())   # true negative
    n_total_all = len(merged_all)

    x_all = merged_all["tpm_pred"].to_numpy()
    y_all = merged_all["tpm_gt"].to_numpy()

    # ----------------------------------------------------------
    # Universe 2: ACTIVE — GT ∪ non-zero predicted
    # Excludes TN pairs. Fair comparison across methods regardless
    # of index size (LK already behaves this way since it only
    # writes non-zero transcripts to its output file).
    # ----------------------------------------------------------
    abund_nz = abund[abund["tpm"] > 0]
    merged_active = abund_nz.merge(gt, on="transcript_id", how="outer",
                                   suffixes=("_pred", "_gt"))
    merged_active["tpm_pred"] = merged_active["tpm_pred"].fillna(0.0)
    merged_active["tpm_gt"]   = merged_active["tpm_gt"].fillna(0.0)

    x_active    = merged_active["tpm_pred"].to_numpy()
    y_active    = merged_active["tpm_gt"].to_numpy()
    n_active    = len(merged_active)

    # ----------------------------------------------------------
    # Universe 3: GT non-zero only
    # Restricts to transcripts present in the simulation GT.
    # ----------------------------------------------------------
    gt_nonzero_mask = merged_all["tpm_gt"] > 0
    x_gt = merged_all.loc[gt_nonzero_mask, "tpm_pred"].to_numpy()
    y_gt = merged_all.loc[gt_nonzero_mask, "tpm_gt"].to_numpy()
    n_gt_nonzero = int(gt_nonzero_mask.sum())

    results = {
        "sample":       sample,
        # Contingency table
        "n_total_all":  n_total_all,
        "n_tp":         n_tp,
        "n_fp":         n_fp,
        "n_fn":         n_fn,
        "n_tn":         n_tn,
        # Universe sizes
        "n_active":     n_active,
        "n_gt_nonzero": n_gt_nonzero,
        # Metrics: all transcripts
        **{f"all_{k}": v for k, v in compute_metrics(x_all, y_all).items()},
        # Metrics: active universe (GT ∪ non-zero pred)
        **{f"active_{k}": v for k, v in compute_metrics(x_active, y_active).items()},
        # Metrics: GT non-zero only
        **{f"gt_{k}": v for k, v in compute_metrics(x_gt, y_gt).items()},
    }
    return results


# ============================================================
# Formatting
# ============================================================

def _metric_block(prefix: str, results: dict) -> list:
    """
    Build a list of formatted metric lines for a given prefix (all/active/gt).

    Args:
        prefix  : str  -- one of "all", "active", "gt".
        results : dict -- output of compare_sample().

    Returns:
        list[str] -- formatted lines (no trailing newline).
    """
    return [
        f"  Spearman    : {fmt(results[f'{prefix}_spearman'])}",
        f"  Pearson     : {fmt(results[f'{prefix}_pearson'])}",
        f"  MAE         : {fmt(results[f'{prefix}_mae'])}",
        f"  RMSE        : {fmt(results[f'{prefix}_rmse'])}",
        f"  Bray-Curtis : {fmt(results[f'{prefix}_bray_curtis'])}",
    ]


def format_results(results: dict) -> str:
    """
    Format a results dict into a human-readable string block.

    Includes contingency table and three metric sections, each labelled
    with the number of transcripts in that universe.

    Args:
        results : dict -- output of compare_sample().

    Returns:
        str -- formatted report block.
    """
    n_all    = results["n_total_all"]
    n_active = results["n_active"]
    n_gt     = results["n_gt_nonzero"]
    tp       = results["n_tp"]
    fp       = results["n_fp"]
    fn       = results["n_fn"]
    tn       = results["n_tn"]

    lines = [
        f"  Sample : {results['sample']}",
        "",
        # --- Contingency table ---
        "  Contingency table",
        "  " + "-" * 50,
        f"  True  positives  (pred>0, GT>0) : {tp:>10,}",
        f"  False positives  (pred>0, GT=0) : {fp:>10,}",
        f"  False negatives  (pred=0, GT>0) : {fn:>10,}",
        f"  True  negatives  (pred=0, GT=0) : {tn:>10,}",
        f"  Total (outer join universe)     : {n_all:>10,}",
        "",
        # --- Metrics: all ---
        f"  --- Metrics (all transcripts, N={n_all:,}) ---",
        *_metric_block("all", results),
        "",
        # --- Metrics: active universe ---
        f"  --- Metrics (active universe: GT ∪ non-zero pred, N={n_active:,}) ---",
        "  (excludes true-negative (0,0) pairs — fair cross-method comparison)",
        *_metric_block("active", results),
        "",
        # --- Metrics: GT non-zero only ---
        f"  --- Metrics (GT non-zero transcripts only, N={n_gt:,}) ---",
        "  (only transcripts present in the simulation ground truth)",
        *_metric_block("gt", results),
    ]
    return "\n".join(lines)


# ============================================================
# Main
# ============================================================

def run_for_experiment(exp: str, results_base: Path, sample_gt_map: dict) -> None:
    """
    Run ground-truth comparison for all matching samples in one experiment folder.

    Args:
        exp            : str         -- experiment folder name.
        results_base   : Path        -- base directory containing experiment folders.
        sample_gt_map  : dict        -- {sample_name: gt_file_path}.
    """
    exp_dir = results_base / exp
    if not exp_dir.is_dir():
        print(f"ERROR: experiment folder not found: {exp_dir}")
        return

    # Discover samples in this experiment that have both abundance.tsv and a GT mapping
    samples = sorted([
        d.name for d in exp_dir.iterdir()
        if d.is_dir() and (d / ABUNDANCE_FILE).exists() and d.name in sample_gt_map
    ])

    skipped = sorted([
        d.name for d in exp_dir.iterdir()
        if d.is_dir() and (d / ABUNDANCE_FILE).exists() and d.name not in sample_gt_map
    ])

    print("=" * 70)
    print(f"GT comparison — {exp}")
    print(f"  Samples with GT : {samples}")
    if skipped:
        print(f"  Skipped (no GT) : {skipped}")
    print("=" * 70)

    all_results = []
    for sample in samples:
        abund_path = exp_dir / sample / ABUNDANCE_FILE
        gt_path    = Path(sample_gt_map[sample])

        if not gt_path.exists():
            print(f"\n[SKIP] {sample} — GT file not found: {gt_path}")
            continue

        print(f"\n{'=' * 60}")
        print(f"{sample}")
        print("=" * 60)
        results = compare_sample(abund_path, gt_path, sample)
        block   = format_results(results)
        print(block)
        all_results.append((sample, block))

    if not all_results:
        print("No samples compared.")
        return

    # Save report
    report_path = exp_dir / f"gt_comparison_{TIMESTAMP}.txt"
    with open(report_path, "w") as fh:
        fh.write(f"GT comparison — {exp}\n")
        fh.write(f"Timestamp: {TIMESTAMP}\n")
        fh.write(f"GT map:\n")
        for s, g in sample_gt_map.items():
            fh.write(f"  {s} -> {g}\n")
        fh.write("\n")
        for sample, block in all_results:
            fh.write(f"\n{'=' * 60}\n{sample}\n{'=' * 60}\n")
            fh.write(block + "\n")
    print(f"\nReport saved: {report_path}")


def main() -> None:
    """
    Parse CLI args and run GT comparison for each provided experiment folder.
    """
    parser = argparse.ArgumentParser(
        description="Compare experiment abundance.tsv files against ground truth TPMs."
    )
    parser.add_argument("experiments", nargs="+",
                        help="One or more experiment folder names.")
    parser.add_argument("--results-base", default=RESULTS_BASE,
                        help=f"Base results directory. Default: {RESULTS_BASE}")
    args = parser.parse_args()

    results_base = Path(args.results_base)

    for exp in args.experiments:
        run_for_experiment(exp, results_base, SAMPLE_GT_MAP)


if __name__ == "__main__":
    main()
