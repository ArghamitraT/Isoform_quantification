"""
plot_fp_tpm_zoom.py
===================
Focused figure: TPM distribution of JK MS false-positive (leaked) transcripts
for sim1, zoomed to the range [1e-3, 10] with a fine linear x-axis.

Reads data from the most-recently created EDA output folder inside JK_EXP.
Saves the figure directly into that EDA folder's figures/ directory.

Inputs:
    - JK MS abundance.tsv  : JK_EXP/sim1/abundance.tsv
    - Ground truth TSV      : SAMPLE_GT_MAP["sim1"]

Outputs:
    - {JK_EXP}/eda_{latest}/figures/fig_fp_tpm_zoom_sim1.png

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/plot_fp_tpm_zoom.py
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================

RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"
JK_EXP       = "exprmnt_2026_03_28__00_51_02"   # JOLI multi-sample experiment
SAMPLE        = "sim1"

GT_PATH = (
    "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/"
    "sim_real_data/ground_truths/PB_sample1_gt.tsv"
)

ABUNDANCE_FILE = "abundance.tsv"

# Zoom window (linear x-axis)
TPM_MIN = 1e-3
TPM_MAX = 1.0

# Fine tick spacing on the linear x-axis
MAJOR_TICK_STEP = 0.1    # major gridlines/ticks every 1 TPM
MINOR_TICK_STEP = 0.01    # minor ticks every 0.1 TPM

# ============================================================
# END CONFIG
# ============================================================

C_FP = "#d01c8b"   # magenta — false positive


# ──────────────────────────────────────────────────────────────
# I/O helpers
# ──────────────────────────────────────────────────────────────

def read_abundance(path: Path) -> pd.DataFrame:
    """Read abundance.tsv and return DataFrame with [transcript_id, tpm]."""
    df = pd.read_csv(path, sep="\t", comment="#")
    id_col  = next((c for c in ["transcript_id","target_id","Name","transcript"]
                    if c in df.columns), df.columns[0])
    val_col = next((c for c in ["tpm","TPM","est_counts","NumReads"]
                    if c in df.columns), df.columns[1])
    out = df[[id_col, val_col]].copy()
    out.columns = ["transcript_id", "tpm"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["tpm"] = pd.to_numeric(out["tpm"], errors="coerce").fillna(0.0)
    return out


def read_gt(path: Path) -> pd.DataFrame:
    """Read ground truth CSV and return DataFrame with [transcript_id, tpm]."""
    df = pd.read_csv(path, index_col=0)
    id_col  = next((c for c in ["transcript_name","transcript_id","target_id","Name"]
                    if c in df.columns), df.columns[0])
    val_col = next((c for c in ["tpm","TPM","est_counts"]
                    if c in df.columns), df.columns[1])
    out = df[[id_col, val_col]].copy()
    out.columns = ["transcript_id", "tpm"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["tpm"] = pd.to_numeric(out["tpm"], errors="coerce").fillna(0.0)
    return out


def classify(abund: pd.DataFrame, gt: pd.DataFrame) -> pd.DataFrame:
    """
    Classify each predicted transcript as TP / FP / FN.

    Args:
        abund : DataFrame [transcript_id, tpm] — predicted abundances (non-zero only used for TP/FP)
        gt    : DataFrame [transcript_id, tpm] — ground truth

    Returns:
        DataFrame with columns [transcript_id, tpm_pred, tpm_gt, label]
        where label is one of "TP", "FP", "FN".
    """
    abund_nz = abund[abund["tpm"] > 0].copy()
    merged = abund_nz.merge(
        gt[["transcript_id","tpm"]].rename(columns={"tpm":"tpm_gt"}),
        on="transcript_id", how="outer"
    )
    merged["tpm"] = merged["tpm"].fillna(0.0)
    merged["tpm_gt"] = merged["tpm_gt"].fillna(0.0)
    merged = merged.rename(columns={"tpm": "tpm_pred"})

    pred_nz = merged["tpm_pred"] > 0
    gt_nz   = merged["tpm_gt"]   > 0
    merged["label"] = np.where(
        pred_nz & gt_nz,  "TP",
        np.where(pred_nz & ~gt_nz, "FP", "FN")
    )
    return merged


# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────

def main():
    jk_dir  = Path(RESULTS_BASE) / JK_EXP
    abund_jk = read_abundance(jk_dir / SAMPLE / ABUNDANCE_FILE)
    gt       = read_gt(Path(GT_PATH))

    cls_jk = classify(abund_jk, gt)
    fp = cls_jk[cls_jk["label"] == "FP"]
    fp_tpm = fp["tpm_pred"].values

    # Restrict to zoom window
    mask   = (fp_tpm >= TPM_MIN) & (fp_tpm <= TPM_MAX)
    fp_zoom = fp_tpm[mask]

    print(f"FP transcripts total          : {len(fp_tpm):,}")
    print(f"FP in [{TPM_MIN}, {TPM_MAX}]  : {len(fp_zoom):,}  "
          f"({100*len(fp_zoom)/len(fp_tpm):.1f}%)")
    print(f"TPM range in window           : [{fp_zoom.min():.4f}, {fp_zoom.max():.4f}]")

    # Find most recent EDA folder to save output
    eda_dirs = sorted(jk_dir.glob("eda_*"))
    if not eda_dirs:
        raise RuntimeError(f"No eda_* folder found in {jk_dir}")
    eda_dir = eda_dirs[-1]
    figures_dir = eda_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output folder: {eda_dir}")

    # ── Figure: histogram + CDF side by side ─────────────────
    # Fine linear bins over [TPM_MIN, TPM_MAX]
    bins = np.arange(TPM_MIN, TPM_MAX + MINOR_TICK_STEP, MINOR_TICK_STEP)

    fig, (ax_hist, ax_cdf) = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        f"FP (leaked) transcript TPM distribution — sim1\n"
        f"Zoom: [{TPM_MIN}, {TPM_MAX}]  |  "
        f"N (in window) = {len(fp_zoom):,} / {len(fp_tpm):,} total FP",
        fontsize=11, fontweight="bold"
    )

    # ── Left: histogram ───────────────────────────────────────
    ax_hist.hist(fp_zoom, bins=bins, color=C_FP, edgecolor="none", density=True)
    ax_hist.set_xlabel("Predicted TPM (JK MS)")
    ax_hist.set_ylabel("Density")
    ax_hist.set_title("Histogram (density)")
    ax_hist.set_xlim(TPM_MIN, TPM_MAX)

    # Major ticks every MAJOR_TICK_STEP, minor ticks every MINOR_TICK_STEP
    ax_hist.xaxis.set_major_locator(ticker.MultipleLocator(MAJOR_TICK_STEP))
    ax_hist.xaxis.set_minor_locator(ticker.MultipleLocator(MINOR_TICK_STEP))
    ax_hist.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax_hist.grid(which="major", alpha=0.4)
    ax_hist.grid(which="minor", alpha=0.15, linestyle=":")
    ax_hist.tick_params(axis="x", rotation=45)

    # ── Right: CDF of the zoomed window ──────────────────────
    fp_zoom_sorted = np.sort(fp_zoom)
    cdf_y = np.arange(1, len(fp_zoom_sorted) + 1) / len(fp_tpm)  # fraction of ALL FP

    ax_cdf.plot(fp_zoom_sorted, cdf_y, color=C_FP, linewidth=1.5)
    ax_cdf.set_xlabel("Predicted TPM (JK MS)")
    ax_cdf.set_ylabel(f"Cumulative fraction of all {len(fp_tpm):,} FP")
    ax_cdf.set_title("CDF (relative to total FP)")
    ax_cdf.set_xlim(TPM_MIN, TPM_MAX)

    ax_cdf.xaxis.set_major_locator(ticker.MultipleLocator(MAJOR_TICK_STEP))
    ax_cdf.xaxis.set_minor_locator(ticker.MultipleLocator(MINOR_TICK_STEP))
    ax_cdf.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax_cdf.grid(which="major", alpha=0.4)
    ax_cdf.grid(which="minor", alpha=0.15, linestyle=":")
    ax_cdf.tick_params(axis="x", rotation=45)

    plt.tight_layout()
    out_path = figures_dir / "fig_fp_tpm_zoom_sim1.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
