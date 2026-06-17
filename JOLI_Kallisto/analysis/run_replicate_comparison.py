"""
run_replicate_comparison.py
===========================
Compare abundance.tsv outputs between two biological replicates.

Since neither replicate is a ground truth, TP/FP/FN/TN are replaced by
concordance categories based on which replicate detects each transcript:

    Concordant positive  (both > 0)        — both replicates detect the transcript
    Rep1-only            (rep1>0, rep2=0)  — FP(rep1) / FN(rep2)
    Rep2-only            (rep1=0, rep2>0)  — FP(rep2) / FN(rep1)
    Concordant negative  (both = 0)        — neither detects the transcript

Precision/recall/F1 are reported from both perspectives.

Three metric universes (Spearman, Pearson, MAE, RMSE, Bray-Curtis):
  1. All transcripts     — outer join (includes concordant negatives)
  2. Active universe     — rep1-nonzero ∪ rep2-nonzero (excludes 0,0 pairs)
  3. Concordant positive — only transcripts detected in both replicates

Inputs:
    Two abundance.tsv paths (passed as CLI positional args).
    Optional --names to label rep1 and rep2.

Outputs:
    - Report printed to stdout.
    - Report saved inside the experiment folder(s) of each abundance.tsv:
        <exp_dir>/replicate_comparison_{timestamp}.txt  (one file per exp dir)

Usage:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5

    # Same experiment folder:
    python analysis/run_replicate_comparison.py \\
        /path/to/exprmnt_ABC/flnc_31/abundance.tsv \\
        /path/to/exprmnt_ABC/flnc_32/abundance.tsv \\
        --names flnc_31 flnc_32

    # Different experiment folders:
    python analysis/run_replicate_comparison.py \\
        /path/to/exprmnt_ABC/flnc_31/abundance.tsv \\
        /path/to/exprmnt_XYZ/flnc_31/abundance.tsv \\
        --names exprmnt_ABC/flnc_31 exprmnt_XYZ/flnc_31
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# CONFIG — edit before running; do not edit below
# ============================================================
RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"
# ============================================================
# END CONFIG
# ============================================================

TIMESTAMP = time.strftime("%Y_%m_%d__%H_%M_%S")


# ============================================================
# I/O
# ============================================================

def read_abundance(path: Path) -> pd.DataFrame:
    """
    Read an abundance TSV (tab-separated).

    Args:
        path : Path -- path to abundance.tsv.

    Returns:
        pd.DataFrame with columns [transcript_id, tpm].
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    id_candidates  = ["transcript_id", "target_id", "Name", "transcript"]
    val_candidates = ["tpm", "TPM", "est_counts", "NumReads"]
    id_col  = next((c for c in id_candidates  if c in df.columns), df.columns[0])
    val_col = next((c for c in val_candidates if c in df.columns), df.columns[1])
    out = df[[id_col, val_col]].copy()
    out.columns = ["transcript_id", "tpm"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["tpm"] = pd.to_numeric(out["tpm"], errors="coerce").fillna(0.0)
    return out


# ============================================================
# Metric functions
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
        x : np.ndarray -- TPM values for rep1.
        y : np.ndarray -- TPM values for rep2.

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


def precision_recall_f1(tp: int, fp: int, fn: int) -> tuple:
    """
    Compute precision, recall, and F1 from TP/FP/FN counts.

    Args:
        tp : int -- true positives (concordant positives).
        fp : int -- false positives (rep-only detections from this rep's perspective).
        fn : int -- false negatives (rep-only detections from the other rep's perspective).

    Returns:
        (precision, recall, f1) : tuple of floats.
    """
    precision = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
    recall    = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    if pd.isna(precision) or pd.isna(recall) or (precision + recall) == 0:
        f1 = float("nan")
    else:
        f1 = 2 * precision * recall / (precision + recall)
    return precision, recall, f1


# ============================================================
# Core comparison
# ============================================================

def compare_replicates(path1: Path, path2: Path,
                       name1: str, name2: str) -> dict:
    """
    Compare two abundance.tsv files (biological replicates).

    Computes:
      - Concordance table: concordant positive, rep1-only, rep2-only, concordant negative
      - Precision/recall/F1 from both perspectives
      - Three metric sets: all / active / concordant-positive

    Args:
        path1 : Path -- abundance.tsv for replicate 1.
        path2 : Path -- abundance.tsv for replicate 2.
        name1 : str  -- label for replicate 1.
        name2 : str  -- label for replicate 2.

    Returns:
        dict with all concordance counts and metrics.
    """
    rep1 = read_abundance(path1)
    rep2 = read_abundance(path2)

    # Outer join — fill missing with 0
    merged = rep1.merge(rep2, on="transcript_id", how="outer",
                        suffixes=("_r1", "_r2"))
    merged["tpm_r1"] = merged["tpm_r1"].fillna(0.0)
    merged["tpm_r2"] = merged["tpm_r2"].fillna(0.0)

    r1_nz = merged["tpm_r1"] > 0
    r2_nz = merged["tpm_r2"] > 0

    # Concordance categories
    n_conc_pos  = int(( r1_nz &  r2_nz).sum())   # both detected
    n_rep1_only = int(( r1_nz & ~r2_nz).sum())   # rep1>0, rep2=0  → FP(rep1) / FN(rep2)
    n_rep2_only = int((~r1_nz &  r2_nz).sum())   # rep2>0, rep1=0  → FP(rep2) / FN(rep1)
    n_conc_neg  = int((~r1_nz & ~r2_nz).sum())   # neither detected
    n_total     = len(merged)

    # Jaccard: concordant positive / union of all detected
    n_union = n_conc_pos + n_rep1_only + n_rep2_only
    jaccard = n_conc_pos / n_union if n_union > 0 else float("nan")

    # Precision / recall / F1 from rep1's perspective
    #   rep1 as "predictor", rep2 as "reference":
    #     FP(rep1) = rep1_only, FN(rep1) = rep2_only
    prec1, rec1, f1_1 = precision_recall_f1(n_conc_pos, n_rep1_only, n_rep2_only)
    # From rep2's perspective:
    #     FP(rep2) = rep2_only, FN(rep2) = rep1_only
    prec2, rec2, f1_2 = precision_recall_f1(n_conc_pos, n_rep2_only, n_rep1_only)

    x_all = merged["tpm_r1"].to_numpy()
    y_all = merged["tpm_r2"].to_numpy()

    # Active universe: rep1-nonzero ∪ rep2-nonzero (excludes 0,0 pairs)
    active_mask = r1_nz | r2_nz
    x_active = merged.loc[active_mask, "tpm_r1"].to_numpy()
    y_active = merged.loc[active_mask, "tpm_r2"].to_numpy()
    n_active = int(active_mask.sum())

    # Concordant positive only
    conc_mask = r1_nz & r2_nz
    x_conc = merged.loc[conc_mask, "tpm_r1"].to_numpy()
    y_conc = merged.loc[conc_mask, "tpm_r2"].to_numpy()

    return {
        "name1": name1,
        "name2": name2,
        # Concordance table
        "n_total":     n_total,
        "n_conc_pos":  n_conc_pos,
        "n_rep1_only": n_rep1_only,
        "n_rep2_only": n_rep2_only,
        "n_conc_neg":  n_conc_neg,
        "n_active":    n_active,
        "jaccard":     jaccard,
        # Precision / recall / F1
        "prec1": prec1, "rec1": rec1, "f1_1": f1_1,
        "prec2": prec2, "rec2": rec2, "f1_2": f1_2,
        # Metrics: all transcripts
        **{f"all_{k}": v for k, v in compute_metrics(x_all, y_all).items()},
        # Metrics: active universe
        **{f"active_{k}": v for k, v in compute_metrics(x_active, y_active).items()},
        # Metrics: concordant positive only
        **{f"conc_{k}": v for k, v in compute_metrics(x_conc, y_conc).items()},
    }


# ============================================================
# Formatting
# ============================================================

def _metric_block(prefix: str, results: dict) -> list:
    """
    Build formatted metric lines for a given universe prefix.

    Args:
        prefix  : str  -- one of "all", "active", "conc".
        results : dict -- output of compare_replicates().

    Returns:
        list[str] -- formatted lines.
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
    Format a results dict into a human-readable report block.

    Args:
        results : dict -- output of compare_replicates().

    Returns:
        str -- formatted report.
    """
    n1 = results["name1"]
    n2 = results["name2"]
    cp = results["n_conc_pos"]
    r1 = results["n_rep1_only"]
    r2 = results["n_rep2_only"]
    cn = results["n_conc_neg"]
    nt = results["n_total"]
    na = results["n_active"]
    nc = cp

    lines = [
        f"  Rep1 : {n1}",
        f"  Rep2 : {n2}",
        "",
        # --- Concordance table ---
        "  Concordance table",
        "  " + "-" * 55,
        f"  Concordant positive  (both > 0)       : {cp:>10,}",
        f"  Rep1-only  (rep1>0, rep2=0) FP(rep1)/FN(rep2) : {r1:>10,}",
        f"  Rep2-only  (rep1=0, rep2>0) FP(rep2)/FN(rep1) : {r2:>10,}",
        f"  Concordant negative  (both = 0)       : {cn:>10,}",
        f"  Total (outer join universe)            : {nt:>10,}",
        f"  Jaccard (conc_pos / union detected)    : {fmt(results['jaccard'])}",
        "",
        # --- Precision / recall / F1 ---
        "  Precision / Recall / F1",
        "  " + "-" * 55,
        f"  From {n1}'s perspective (rep1=pred, rep2=ref):",
        f"    Precision : {fmt(results['prec1'])}   (conc_pos / (conc_pos + rep1_only))",
        f"    Recall    : {fmt(results['rec1'])}   (conc_pos / (conc_pos + rep2_only))",
        f"    F1        : {fmt(results['f1_1'])}",
        f"  From {n2}'s perspective (rep2=pred, rep1=ref):",
        f"    Precision : {fmt(results['prec2'])}   (conc_pos / (conc_pos + rep2_only))",
        f"    Recall    : {fmt(results['rec2'])}   (conc_pos / (conc_pos + rep1_only))",
        f"    F1        : {fmt(results['f1_2'])}",
        "",
        # --- Metrics: all ---
        f"  --- Metrics (all transcripts, N={nt:,}) ---",
        "  (includes concordant negatives — inflated by index size)",
        *_metric_block("all", results),
        "",
        # --- Metrics: active ---
        f"  --- Metrics (active universe: rep1 ∪ rep2 non-zero, N={na:,}) ---",
        "  (excludes concordant negatives — fair cross-method comparison)",
        *_metric_block("active", results),
        "",
        # --- Metrics: concordant positive ---
        f"  --- Metrics (concordant positive only, N={nc:,}) ---",
        "  (only transcripts detected in both replicates — strictest reproducibility view)",
        *_metric_block("conc", results),
    ]
    return "\n".join(lines)


# ============================================================
# Output path helpers
# ============================================================

def get_exp_dir(abundance_path: Path) -> Path:
    """
    Infer the experiment folder from an abundance.tsv path.
    Walks up the path tree to find the first directory named exprmnt_*.

    Args:
        abundance_path : Path -- path to abundance.tsv.

    Returns:
        Path -- experiment directory, or parent of the abundance file if not found.
    """
    for parent in abundance_path.parents:
        if parent.name.startswith("exprmnt_"):
            return parent
    return abundance_path.parent


# ============================================================
# Main
# ============================================================

def main() -> None:
    """
    Parse CLI args, run replicate comparison, print and save report.

    Run:
        cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
        conda activate NanoCount_5

        python analysis/run_replicate_comparison.py \\
            /path/to/exprmnt_ABC/flnc_31/abundance.tsv \\
            /path/to/exprmnt_ABC/flnc_32/abundance.tsv \\
            --names flnc_31 flnc_32
    """
    parser = argparse.ArgumentParser(
        description="Compare abundance.tsv between two biological replicates."
    )
    parser.add_argument("rep1", type=Path, help="Path to abundance.tsv for replicate 1.")
    parser.add_argument("rep2", type=Path, help="Path to abundance.tsv for replicate 2.")
    parser.add_argument("--names", nargs=2, default=None,
                        metavar=("NAME1", "NAME2"),
                        help="Labels for rep1 and rep2 (default: derived from paths).")
    args = parser.parse_args()

    path1 = args.rep1.resolve()
    path2 = args.rep2.resolve()

    for p in [path1, path2]:
        if not p.exists():
            print(f"ERROR: file not found: {p}")
            raise SystemExit(1)

    name1 = args.names[0] if args.names else path1.parent.name
    name2 = args.names[1] if args.names else path2.parent.name

    print("=" * 65)
    print("Replicate comparison")
    print(f"  Rep1 : {name1}  ({path1})")
    print(f"  Rep2 : {name2}  ({path2})")
    print("=" * 65)

    results = compare_replicates(path1, path2, name1, name2)
    block   = format_results(results)
    print(block)

    # Save report inside each unique experiment folder
    exp_dirs = list({get_exp_dir(path1), get_exp_dir(path2)})
    for exp_dir in exp_dirs:
        exp_dir.mkdir(parents=True, exist_ok=True)
        report_path = exp_dir / f"replicate_comparison_{TIMESTAMP}.txt"
        with open(report_path, "w") as fh:
            fh.write(f"Replicate comparison\n")
            fh.write(f"Timestamp : {TIMESTAMP}\n")
            fh.write(f"Rep1      : {name1}  ({path1})\n")
            fh.write(f"Rep2      : {name2}  ({path2})\n\n")
            fh.write(block + "\n")
        print(f"\nReport saved: {report_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
