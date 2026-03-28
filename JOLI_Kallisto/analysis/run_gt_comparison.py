"""
run_gt_comparison.py
====================
Compare abundance.tsv outputs from one or more experiment folders against
ground-truth TPM files (simulation ground truth).

Auto-discovers samples inside the experiment folder. For each sample that has
a matching entry in SAMPLE_GT_MAP, loads both files and computes Spearman,
Pearson, MAE, RMSE, Bray-Curtis, and JS-distance metrics.

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
    "sim1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv",
    "sim2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv",
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


# ============================================================
# Per-sample comparison
# ============================================================

def compare_sample(abund_path: Path, gt_path: Path, sample: str) -> dict:
    """
    Compare one sample's abundance.tsv against its ground truth.

    Merges on transcript_id (outer join, zeros for missing), then computes
    metrics on (a) all transcripts and (b) transcripts non-zero in ground truth.

    Args:
        abund_path : Path -- path to abundance.tsv.
        gt_path    : Path -- path to ground truth CSV.
        sample     : str  -- sample name (for printing).

    Returns:
        dict with all computed metrics.
    """
    abund = read_abundance(abund_path)
    gt    = read_ground_truth(gt_path)

    merged = abund.merge(gt, on="transcript_id", how="outer",
                         suffixes=("_pred", "_gt"))
    merged["tpm_pred"] = merged["tpm_pred"].fillna(0.0)
    merged["tpm_gt"]   = merged["tpm_gt"].fillna(0.0)

    x_all = merged["tpm_pred"].to_numpy()
    y_all = merged["tpm_gt"].to_numpy()

    # Restrict to transcripts with non-zero ground truth
    gt_nonzero = merged["tpm_gt"] > 0
    x_gt = merged.loc[gt_nonzero, "tpm_pred"].to_numpy()
    y_gt = merged.loc[gt_nonzero, "tpm_gt"].to_numpy()

    # Transcript count stats
    n_pred_nonzero = int((merged["tpm_pred"] > 0).sum())
    n_gt_nonzero   = int((merged["tpm_gt"]   > 0).sum())
    n_both_nonzero = int(((merged["tpm_pred"] > 0) & (merged["tpm_gt"] > 0)).sum())
    n_total        = len(merged)

    results = {
        "sample":          sample,
        "n_total":         n_total,
        "n_pred_nonzero":  n_pred_nonzero,
        "n_gt_nonzero":    n_gt_nonzero,
        "n_both_nonzero":  n_both_nonzero,
        # All-transcript metrics
        "all_spearman":    spearman_corr(x_all, y_all),
        "all_pearson":     pearson_corr(x_all, y_all),
        "all_mae":         mae(x_all, y_all),
        "all_rmse":        rmse(x_all, y_all),
        "all_bray_curtis": bray_curtis(x_all, y_all),
        # GT-nonzero metrics (more meaningful — only transcripts present in simulation)
        "gt_spearman":     spearman_corr(x_gt, y_gt),
        "gt_pearson":      pearson_corr(x_gt, y_gt),
        "gt_mae":          mae(x_gt, y_gt),
        "gt_rmse":         rmse(x_gt, y_gt),
        "gt_bray_curtis":  bray_curtis(x_gt, y_gt),
    }
    return results


def format_results(results: dict) -> str:
    """
    Format a results dict into a human-readable string block.

    Args:
        results : dict -- output of compare_sample().

    Returns:
        str -- formatted report block.
    """
    lines = [
        f"  Sample            : {results['sample']}",
        f"  Total transcripts : {results['n_total']}",
        f"  Non-zero pred     : {results['n_pred_nonzero']}",
        f"  Non-zero GT       : {results['n_gt_nonzero']}",
        f"  Non-zero in both  : {results['n_both_nonzero']}",
        "",
        "  --- Metrics (all transcripts) ---",
        f"  Spearman    : {fmt(results['all_spearman'])}",
        f"  Pearson     : {fmt(results['all_pearson'])}",
        f"  MAE         : {fmt(results['all_mae'])}",
        f"  RMSE        : {fmt(results['all_rmse'])}",
        f"  Bray-Curtis : {fmt(results['all_bray_curtis'])}",
        "",
        "  --- Metrics (GT non-zero transcripts only) ---",
        f"  Spearman    : {fmt(results['gt_spearman'])}",
        f"  Pearson     : {fmt(results['gt_pearson'])}",
        f"  MAE         : {fmt(results['gt_mae'])}",
        f"  RMSE        : {fmt(results['gt_rmse'])}",
        f"  Bray-Curtis : {fmt(results['gt_bray_curtis'])}",
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

        print(f"\n--- {sample} ---")
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
