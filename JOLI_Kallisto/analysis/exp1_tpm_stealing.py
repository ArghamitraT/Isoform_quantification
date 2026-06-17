"""
exp1_tpm_stealing.py
====================
Experiment 1: Does the Dirichlet prior cause FP transcripts to steal TPM
from true positive (TP) transcripts?

Hypothesis: In JK MS (MAP EM with Dirichlet prior), FP transcripts are kept
alive by the prior. Since total TPM = 1M, any TPM assigned to FPs must come
from somewhere — the claim is it comes from TPs (shared EC competition).

Test: Compare sum(TPM for TPs) in JK MS vs JK single (plain EM, no prior).
If the hypothesis holds:
    sum(TPM for TPs) in JK MS  <  sum(TPM for TPs) in JK single
    sum(TPM for FPs) in JK MS  >  sum(TPM for FPs) in JK single (≈ 0)

Inputs:
    - JK MS abundance.tsv   : exprmnt_2026_03_30__22_37_55 / sim1 & sim2
    - JK single abundance.tsv: exprmnt_2026_03_30__22_41_19 / sim1 & sim2
    - Ground truth TSVs      : sample1_gt.tsv, sample2_gt.tsv

Outputs (saved inside exprmnt_2026_03_30__22_37_55/figures/):
    - exp1_tpm_stealing_{timestamp}.png  : bar chart + scatter
    - exp1_tpm_stealing_{timestamp}.txt  : numerical summary

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate Joli_kallisto
    python analysis/exp1_tpm_stealing.py
"""

import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ============================================================
# CONFIG — edit here; do not edit below
# ============================================================

RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

JKMS_EXPRMNT   = "exprmnt_2026_03_30__22_37_55"   # JK multi-sample MAP EM
JKSINGLE_EXPRMNT = "exprmnt_2026_03_30__22_41_19"  # JK single-sample plain EM

GT_BASE = "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths"

SAMPLES = {
    "sim1": "PB_sample1_gt.tsv",
    "sim2": "PB_sample2_gt.tsv",
}

# TPM threshold to call a transcript "expressed" (nonzero)
TPM_THRESHOLD = 0.0

# ============================================================
# END CONFIG
# ============================================================


def load_abundance(path: str) -> pd.Series:
    """
    Load abundance.tsv and return a Series: transcript_id -> tpm.

    Args:
        path : str -- Path to abundance.tsv.

    Returns:
        pd.Series -- Index = transcript_id, values = tpm (float).
    """
    df = pd.read_csv(path, sep="\t")
    return df.set_index("target_id")["tpm"]


def load_ground_truth(path: str) -> pd.Series:
    """
    Load ground truth TSV and return a Series: transcript_name -> tpm.

    Args:
        path : str -- Path to ground truth TSV (columns: transcript_name, tpm).

    Returns:
        pd.Series -- Index = transcript_name, values = tpm (float).
    """
    # PB_sample*_gt.tsv files are comma-separated despite the .tsv extension;
    # sample*_gt.tsv files are tab-separated. Sniff the separator automatically.
    with open(path) as fh:
        first_line = fh.readline()
    sep = "," if first_line.count(",") > first_line.count("\t") else "\t"
    df = pd.read_csv(path, sep=sep, index_col=0)
    # Keep only expressed transcripts (tpm > 0)
    df = df[df["tpm"] > 0]
    return df.set_index("transcript_name")["tpm"]


def classify_transcripts(pred: pd.Series, gt: pd.Series, threshold: float = 0.0):
    """
    Classify transcripts into TP, FP, FN, TN.

    Args:
        pred      : pd.Series -- Predicted TPM (transcript_id -> tpm).
        gt        : pd.Series -- Ground truth TPM (transcript_id -> tpm, GT nonzero only).
        threshold : float     -- TPM threshold to call a transcript expressed.

    Returns:
        dict with keys: tp_ids, fp_ids, fn_ids
                        tp_tpm (predicted TPM for TPs),
                        fp_tpm (predicted TPM for FPs),
                        fn_tpm (GT TPM for FNs)
    """
    pred_nonzero = set(pred[pred > threshold].index)
    gt_nonzero   = set(gt.index)

    tp_ids = pred_nonzero & gt_nonzero
    fp_ids = pred_nonzero - gt_nonzero
    fn_ids = gt_nonzero - pred_nonzero

    return {
        "tp_ids":  tp_ids,
        "fp_ids":  fp_ids,
        "fn_ids":  fn_ids,
        "tp_tpm":  pred[list(tp_ids)].sum(),
        "fp_tpm":  pred[list(fp_ids)].sum(),
        "fn_tpm":  gt[list(fn_ids)].sum(),
        "tp_count": len(tp_ids),
        "fp_count": len(fp_ids),
        "fn_count": len(fn_ids),
    }


def run_experiment():
    """
    Main experiment: compare TP/FP TPM sums between JK MS and JK single.

    Run:
        cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
        python analysis/exp1_tpm_stealing.py
    """
    timestamp = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")

    # Output goes into JK MS experiment figures/ folder (Rule 3: single figure → figures/)
    out_dir = Path(RESULTS_BASE) / JKMS_EXPRMNT / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    lines = []   # collect text output
    lines.append(f"Experiment 1: TPM Stealing by FP Transcripts")
    lines.append(f"Timestamp: {timestamp}")
    lines.append(f"JK MS    : {JKMS_EXPRMNT}")
    lines.append(f"JK Single: {JKSINGLE_EXPRMNT}")
    lines.append("")

    results = {}   # sample -> dict of metrics

    for sample, gt_file in SAMPLES.items():
        print(f"\n{'='*60}")
        print(f"Sample: {sample}")
        print(f"{'='*60}")

        # Load data
        jkms_path   = Path(RESULTS_BASE) / JKMS_EXPRMNT   / sample / "abundance.tsv"
        jksingle_path = Path(RESULTS_BASE) / JKSINGLE_EXPRMNT / sample / "abundance.tsv"
        gt_path     = Path(GT_BASE) / gt_file

        jkms   = load_abundance(str(jkms_path))
        jksingle = load_abundance(str(jksingle_path))
        gt     = load_ground_truth(str(gt_path))

        print(f"  JK MS   nonzero: {(jkms > TPM_THRESHOLD).sum()}")
        print(f"  JK single nonzero: {(jksingle > TPM_THRESHOLD).sum()}")
        print(f"  GT nonzero:      {len(gt)}")

        # Classify for each method
        ms  = classify_transcripts(jkms, gt, TPM_THRESHOLD)
        sng = classify_transcripts(jksingle, gt, TPM_THRESHOLD)

        # Total TPM check
        total_ms  = jkms.sum()
        total_sng = jksingle.sum()

        lines.append(f"{'='*60}")
        lines.append(f"Sample: {sample}")
        lines.append(f"{'='*60}")
        lines.append(f"{'Metric':<35} {'JK MS':>12} {'JK Single':>12} {'Diff (MS-Single)':>18}")
        lines.append(f"{'-'*80}")

        rows = [
            ("Total TPM (sanity check)",   total_ms,    total_sng),
            ("TP count",                   ms["tp_count"],  sng["tp_count"]),
            ("FP count",                   ms["fp_count"],  sng["fp_count"]),
            ("FN count",                   ms["fn_count"],  sng["fn_count"]),
            ("sum(TPM) for TPs",           ms["tp_tpm"],    sng["tp_tpm"]),
            ("sum(TPM) for FPs",           ms["fp_tpm"],    sng["fp_tpm"]),
            ("sum(TPM) for FNs (GT TPM)",  ms["fn_tpm"],    sng["fn_tpm"]),
            ("TP TPM as % of total",       100*ms["tp_tpm"]/total_ms,  100*sng["tp_tpm"]/total_sng),
            ("FP TPM as % of total",       100*ms["fp_tpm"]/total_ms,  100*sng["fp_tpm"]/total_sng),
        ]

        for label, val_ms, val_sng in rows:
            diff = val_ms - val_sng
            lines.append(f"{label:<35} {val_ms:>12.2f} {val_sng:>12.2f} {diff:>+18.2f}")
            print(f"  {label:<35} MS={val_ms:.2f}  Single={val_sng:.2f}  diff={diff:+.2f}")

        lines.append("")

        results[sample] = {
            "ms": ms, "sng": sng,
            "total_ms": total_ms, "total_sng": total_sng,
            "jkms": jkms, "jksingle": jksingle, "gt": gt,
        }

    # ---- Figure: 2 panels per sample (bar chart of TPM sums) ----
    # Layout: 2 rows (sim1, sim2) × 2 cols (TPM breakdown bar + TP scatter)
    fig = plt.figure(figsize=(14, 8))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.35)
    fig.suptitle("Experiment 1: Does the Dirichlet Prior Steal TPM from TPs?", fontsize=13)

    for row_idx, (sample, gt_file) in enumerate(SAMPLES.items()):
        r      = results[sample]
        ms     = r["ms"]
        sng    = r["sng"]
        total  = r["total_ms"]   # ~1M for both

        # --- Panel left: stacked bar of TP / FP / other TPM ---
        ax_bar = fig.add_subplot(gs[row_idx, 0])

        categories  = ["JK MS", "JK Single"]
        tp_vals  = [ms["tp_tpm"],    sng["tp_tpm"]]
        fp_vals  = [ms["fp_tpm"],    sng["fp_tpm"]]
        other_vals = [total - ms["tp_tpm"] - ms["fp_tpm"],
                      r["total_sng"] - sng["tp_tpm"] - sng["fp_tpm"]]

        x = np.arange(len(categories))
        w = 0.5
        b1 = ax_bar.bar(x, tp_vals,    w, label="TP",    color="#2196F3")
        b2 = ax_bar.bar(x, fp_vals,    w, bottom=tp_vals, label="FP", color="#F44336")
        b3 = ax_bar.bar(x, other_vals, w,
                        bottom=[tp_vals[i]+fp_vals[i] for i in range(2)],
                        label="FN/other", color="#9E9E9E")

        ax_bar.set_title(f"{sample} — TPM breakdown", fontsize=10)
        ax_bar.set_ylabel("TPM")
        ax_bar.set_xticks(x)
        ax_bar.set_xticklabels(categories)
        ax_bar.legend(fontsize=8)

        # Annotate TP bar with exact value and % diff
        pct_diff = 100 * (ms["tp_tpm"] - sng["tp_tpm"]) / sng["tp_tpm"]
        ax_bar.text(0, ms["tp_tpm"] / 2,
                    f"{ms['tp_tpm']:.0f}", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        ax_bar.text(1, sng["tp_tpm"] / 2,
                    f"{sng['tp_tpm']:.0f}", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")
        ax_bar.set_title(
            f"{sample} — TPM breakdown\nTP TPM diff: {pct_diff:+.2f}%", fontsize=9)

        # --- Panel right: scatter of TP TPM (JK MS vs JK Single) ---
        ax_sc = fig.add_subplot(gs[row_idx, 1])

        tp_ids   = list(ms["tp_ids"] & sng["tp_ids"])   # TPs detected by BOTH
        ms_tpm   = r["jkms"][tp_ids].values
        sng_tpm  = r["jksingle"][tp_ids].values

        # Ratio: JK MS / JK Single — values < 1 mean MS is lower (stolen)
        ratio = ms_tpm / np.maximum(sng_tpm, 1e-10)
        pct_below_1 = 100 * (ratio < 1).mean()

        ax_sc.scatter(sng_tpm, ms_tpm, s=1, alpha=0.3, color="#555")
        lim = max(ms_tpm.max(), sng_tpm.max()) * 1.05
        ax_sc.plot([0, lim], [0, lim], "r--", lw=1, label="y=x (equal)")
        ax_sc.set_xlabel("JK Single TPM", fontsize=9)
        ax_sc.set_ylabel("JK MS TPM", fontsize=9)
        ax_sc.set_title(
            f"{sample} — TP scatter (TPs in both)\n"
            f"{pct_below_1:.1f}% of TPs: JK MS < JK Single",
            fontsize=9)
        ax_sc.legend(fontsize=8)
        ax_sc.set_xscale("symlog", linthresh=1)
        ax_sc.set_yscale("symlog", linthresh=1)

    fig_path = out_dir / f"exp1_tpm_stealing_{timestamp}.png"
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nFigure saved: {fig_path}")

    # Save text summary
    txt_path = out_dir / f"exp1_tpm_stealing_{timestamp}.txt"
    with open(txt_path, "w") as fh:
        fh.write("\n".join(lines))
    print(f"Summary saved: {txt_path}")


if __name__ == "__main__":
    run_experiment()
