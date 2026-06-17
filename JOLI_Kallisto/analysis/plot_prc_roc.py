"""
plot_prc_roc.py
===============
Plot Precision-Recall (PRC) and ROC curves for multiple experiments,
comparing transcript-level detection performance across methods.

For each transcript:
  - Label  : 1 if GT TPM > 0 (truly expressed), 0 otherwise
  - Score  : predicted TPM (higher = more confident it is expressed)

Curves are computed by sweeping the TPM threshold from high to low.
One figure per sample, each with two subplots: PRC (left) and ROC (right).

Inputs:
    - abundance.tsv from each experiment/{sample}/
    - Ground truth TSV from GT_BASE

Outputs:
    {OUTPUT_EXP}/figures/prc_roc_{sample}_{timestamp}.png

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate Joli_kallisto
    python analysis/plot_prc_roc.py
"""

from __future__ import annotations

import sys
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve, auc

# ============================================================
# CONFIG — edit here; do not edit below
# ============================================================

RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"
GT_BASE      = "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths"

# Experiments to compare — (label, experiment_folder)
EXPERIMENTS = [
    ("LK",    "exprmnt_2026_05_18__12_42_22"),
    ("JK SS", "exprmnt_2026_05_18__12_46_25"),
    ("JK MS", "exprmnt_2026_05_20__13_04_38"),
]

# Figures saved into this experiment's figures/ folder
OUTPUT_EXP = "exprmnt_2026_05_20__13_04_38"

SAMPLE_GT_MAP = {
    "sim1": "PB_sample1_gt.tsv",
    "sim2": "PB_sample2_gt.tsv",
}

# ============================================================
# END CONFIG
# ============================================================

TIMESTAMP = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
COLORS    = ["#1565C0", "#C62828", "#2E7D32"]   # blue, red, green


def read_ground_truth(path: Path) -> pd.DataFrame:
    """
    Read ground truth TSV/CSV, return DataFrame with columns [tid, tpm_gt].

    Args:
        path : Path -- ground truth file path.

    Returns:
        pd.DataFrame -- columns [tid, tpm_gt].
    """
    with open(path) as f:
        first = f.readline()
    sep = "," if first.count(",") > first.count("\t") else "\t"
    df  = pd.read_csv(path, sep=sep, index_col=0)
    id_col  = next(c for c in ["transcript_name","transcript_id","target_id","Name"]  if c in df.columns)
    val_col = next(c for c in ["tpm","TPM","est_counts"] if c in df.columns)
    out = df[[id_col, val_col]].rename(columns={id_col:"tid", val_col:"tpm_gt"})
    out["tpm_gt"] = pd.to_numeric(out["tpm_gt"], errors="coerce").fillna(0.0)
    return out


def read_abundance(path: Path) -> pd.DataFrame:
    """
    Read abundance.tsv, return DataFrame with columns [tid, tpm_pred].

    Args:
        path : Path -- abundance.tsv file path.

    Returns:
        pd.DataFrame -- columns [tid, tpm_pred].
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    id_col  = "target_id" if "target_id" in df.columns else df.columns[0]
    tpm_col = "tpm"       if "tpm"       in df.columns else df.columns[-1]
    out = df[[id_col, tpm_col]].rename(columns={id_col:"tid", tpm_col:"tpm_pred"})
    out["tpm_pred"] = pd.to_numeric(out["tpm_pred"], errors="coerce").fillna(0.0)
    return out


def compute_curves(
    tpm_pred: np.ndarray,
    labels:   np.ndarray,
) -> dict:
    """
    Compute PR and ROC curves plus summary scalars.

    Args:
        tpm_pred : np.ndarray -- predicted TPM scores.
        labels   : np.ndarray -- binary labels (1 = GT > 0, 0 = GT = 0).

    Returns:
        dict with keys: precision, recall, pr_auc, fpr, tpr, roc_auc,
                        operating_precision, operating_recall, operating_fpr, operating_tpr.
    """
    # PR curve
    precision, recall, pr_thresh = precision_recall_curve(labels, tpm_pred)
    pr_auc_val = auc(recall, precision)

    # ROC curve
    fpr, tpr, roc_thresh = roc_curve(labels, tpm_pred)
    roc_auc_val = auc(fpr, tpr)

    # Operating point: threshold = 0 (any pred > 0 is called expressed)
    op_mask = tpm_pred > 0
    tp = int(( op_mask & (labels == 1)).sum())
    fp = int(( op_mask & (labels == 0)).sum())
    fn = int((~op_mask & (labels == 1)).sum())
    tn = int((~op_mask & (labels == 0)).sum())
    op_prec   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    op_recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    op_fpr    = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    op_tpr    = op_recall

    return {
        "precision":          precision,
        "recall":             recall,
        "pr_auc":             pr_auc_val,
        "fpr":                fpr,
        "tpr":                tpr,
        "roc_auc":            roc_auc_val,
        "operating_precision": op_prec,
        "operating_recall":    op_recall,
        "operating_fpr":       op_fpr,
        "operating_tpr":       op_tpr,
    }


def plot_sample(
    sample:   str,
    results:  list[tuple[str, dict]],
    out_path: Path,
) -> None:
    """
    Plot PRC and ROC side by side for one sample, all methods overlaid.

    Args:
        sample   : str                       -- sample name (for title).
        results  : list[(label, curves_dict)] -- per-method curve data.
        out_path : Path                       -- output figure path.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(f"PRC & ROC — {sample}", fontsize=13)

    ax_pr, ax_roc = axes

    for (label, curves), color in zip(results, COLORS):
        # PRC
        ax_pr.plot(
            curves["recall"], curves["precision"],
            color=color, linewidth=1.8, label=f"{label} (AUC={curves['pr_auc']:.3f})"
        )
        # Operating point on PRC
        ax_pr.scatter(
            curves["operating_recall"], curves["operating_precision"],
            color=color, s=60, zorder=5, marker="o"
        )

        # ROC
        ax_roc.plot(
            curves["fpr"], curves["tpr"],
            color=color, linewidth=1.8, label=f"{label} (AUC={curves['roc_auc']:.3f})"
        )
        # Operating point on ROC
        ax_roc.scatter(
            curves["operating_fpr"], curves["operating_tpr"],
            color=color, s=60, zorder=5, marker="o"
        )

    # PRC formatting
    ax_pr.set_xlabel("Recall", fontsize=10)
    ax_pr.set_ylabel("Precision", fontsize=10)
    ax_pr.set_title("Precision-Recall Curve", fontsize=11)
    ax_pr.legend(fontsize=9, loc="lower left")
    ax_pr.set_xlim([0, 1]); ax_pr.set_ylim([0, 1.02])
    ax_pr.grid(True, alpha=0.3)
    ax_pr.annotate("● = operating point (TPM > 0)", xy=(0.02, 0.04),
                   xycoords="axes fraction", fontsize=7, color="gray")

    # ROC formatting
    ax_roc.plot([0, 1], [0, 1], "k--", linewidth=0.8, alpha=0.4)
    ax_roc.set_xlabel("False Positive Rate", fontsize=10)
    ax_roc.set_ylabel("True Positive Rate", fontsize=10)
    ax_roc.set_title("ROC Curve", fontsize=11)
    ax_roc.legend(fontsize=9, loc="lower right")
    ax_roc.set_xlim([0, 1]); ax_roc.set_ylim([0, 1.02])
    ax_roc.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def run() -> None:
    """
    Load abundance and GT for each experiment and sample, compute curves, save figures.
    """
    out_dir = Path(RESULTS_BASE) / OUTPUT_EXP / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    for sample, gt_file in SAMPLE_GT_MAP.items():
        print(f"\n{'='*60}\nSample: {sample}\n{'='*60}")

        gt_df = read_ground_truth(Path(GT_BASE) / gt_file)
        print(f"  GT non-zero: {(gt_df.tpm_gt > 0).sum():,}")

        results = []
        for label, exp in EXPERIMENTS:
            ab_path = Path(RESULTS_BASE) / exp / sample / "abundance.tsv"
            ab_df   = read_abundance(ab_path)

            merged  = ab_df.merge(gt_df, on="tid", how="outer").fillna(0.0)
            scores  = merged["tpm_pred"].to_numpy(dtype=np.float64)
            labels  = (merged["tpm_gt"] > 0).astype(int).to_numpy()

            print(f"\n  [{label}]  N={len(merged):,}  positives={labels.sum():,}  negatives={(1-labels).sum():,}")
            curves = compute_curves(scores, labels)
            print(f"    PR-AUC={curves['pr_auc']:.4f}  ROC-AUC={curves['roc_auc']:.4f}")
            print(f"    Operating point: precision={curves['operating_precision']:.4f}  "
                  f"recall={curves['operating_recall']:.4f}  FPR={curves['operating_fpr']:.4f}")
            results.append((label, curves))

        out_path = out_dir / f"prc_roc_{sample}_{TIMESTAMP}.png"
        plot_sample(sample, results, out_path)

    print(f"\nDone. Figures saved to: {out_dir}")


if __name__ == "__main__":
    run()
