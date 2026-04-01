"""
run_eda.py
==========
EDA for simulated long-read RNA-seq data: Priority 1 and Priority 6.

Priority 1 — Index size discrepancy:
  Compares the transcript universes of LK (lr-kallisto) and JK MS (JOLI multi-sample)
  against ground truth. Shows how the metric (Spearman) changes with universe definition:
  all transcripts vs active universe (GT ∪ non-zero pred) vs GT-nonzero only.

Priority 6 — Leaked transcript (false positive) analysis:
  Classifies JK MS predicted transcripts as TP / FP / FN / TN relative to GT.
  Investigates predicted TPM distribution, rank position, EC membership, and
  comparison with LK for the 16K false positive transcripts.
  Also produces a continuous-color theta scatter plot (GT theta vs predicted theta,
  colored by GT theta magnitude: white = near zero, red = high).

Inputs:
    - LK experiment folder (abundance.tsv per sample subfolder)
    - JK MS experiment folder (abundance.tsv per sample subfolder)
    - Ground truth TSV files (CSV format: index, transcript_name, tpm)
    - bustools matrix.ec files per sample (for EC size figure)

Outputs (saved inside each experiment folder's eda_{timestamp}/):
    figures/                       — main GT vs predicted TPM scatter plots
      fig_gt_vs_pred_{sample}.png          — GT TPM vs predicted TPM (all GT, log-log)
      fig_gt_vs_pred_classified_{sample}.png — same, colored by TP (green) / FN (orange)
    extra_fig/                     — all supporting / diagnostic figures
      fig_p1_{sample}.png
      fig_p6a_tpm_{sample}.png
      fig_p6b_rank_lk_{sample}.png
      fig_p6c_ec_sizes_{sample}.png
      fig_p6d_theta_scatter_{sample}.png
    priority1_universe_table.txt
    priority6_leaked_stats.txt
    running.log
    runtime.txt
"""

import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from scipy import stats


# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================

RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# Experiment folder names (inside RESULTS_BASE)
LK_EXP = "exprmnt_2026_03_30__22_41_19"    # lr-kallisto baseline
JK_EXP =  "exprmnt_2026_03_30__22_37_55"  # "exprmnt_2026_03_29__18_14_15"    # JOLI multi-sample

# Ground truth files
SAMPLE_GT_MAP = {
    "sim1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv",
    "sim2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv",
}

# bustools output dirs (contain matrix.ec) — for EC size figure (Priority 6)
SAMPLE_BUSTOOLS_DIRS = {
    "sim1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/kallisto_output/ds_100_num1_aln_01_long",
    "sim2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/kallisto_output/ds_100_num1_aln_21_long",
}

# Which priorities to run
RUN_PRIORITY_1 = True
RUN_PRIORITY_6 = True

ABUNDANCE_FILE = "abundance.tsv"

# ============================================================
# END CONFIG
# ============================================================


# ── color palette ──────────────────────────────────────────────
WHITE_TO_RED = LinearSegmentedColormap.from_list("white_red", ["#ffffff", "#cc0000"])
C_LK   = "#2166ac"   # blue for LK
C_JK   = "#d6604d"   # red-orange for JK MS
C_TP   = "#4dac26"   # green  — true positive
C_FP   = "#d01c8b"   # magenta — false positive (leaked)
C_FN   = "#f1a340"   # amber  — false negative


# ============================================================
# I/O helpers
# ============================================================

def read_abundance(path: Path) -> pd.DataFrame:
    """
    Read abundance.tsv (tab-separated). Auto-detects transcript_id and tpm columns.

    Args:
        path : Path -- path to abundance.tsv file.

    Returns:
        pd.DataFrame with columns [transcript_id, tpm].
    """
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
    """
    Read ground truth CSV (comma-separated: unnamed_index, transcript_name, tpm).

    Args:
        path : Path -- path to ground truth TSV/CSV.

    Returns:
        pd.DataFrame with columns [transcript_id, tpm].
    """
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


def read_matrix_ec(bustools_dir: str) -> dict:
    """
    Parse bustools matrix.ec file into a dict mapping EC index → list of transcript indices.

    Args:
        bustools_dir : str -- path to bustools output directory containing matrix.ec.

    Returns:
        dict {ec_id (int): [tx_idx (int), ...]}
    """
    ec_path = Path(bustools_dir) / "matrix.ec"
    ec_map = {}
    with open(ec_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            _, tx_str = line.split("\t")
            ec_idx = len(ec_map)
            ec_map[ec_idx] = [int(x) for x in tx_str.split(",")]
    print(f"  [read_matrix_ec] Loaded {len(ec_map)} ECs from {ec_path}")
    return ec_map


def read_transcripts_txt(bustools_dir: str) -> list:
    """
    Read transcripts.txt (one transcript name per line).

    Args:
        bustools_dir : str -- path to bustools output directory.

    Returns:
        list[str] of transcript names, indexed by transcript position.
    """
    tx_path = Path(bustools_dir) / "transcripts.txt"
    with open(tx_path) as fh:
        names = [line.strip() for line in fh if line.strip()]
    print(f"  [read_transcripts_txt] {len(names)} transcripts from {tx_path}")
    return names


# ============================================================
# Metric helpers
# ============================================================

def spearman(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman correlation between x and y. Returns nan if degenerate."""
    if len(x) < 2:
        return float("nan")
    r, _ = stats.spearmanr(x, y)
    return float(r)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation. Returns nan if degenerate."""
    if len(x) < 2 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return float("nan")
    r, _ = stats.pearsonr(x, y)
    return float(r)


def build_gt_all(abund: pd.DataFrame, gt: pd.DataFrame) -> pd.DataFrame:
    """
    Left-join from GT to predictions: all GT transcripts, fill pred=0 for FN.

    "GT all" universe — every transcript in the GT file is included regardless
    of whether the method predicted it (FN get tpm_pred=0).  FP transcripts
    (pred>0, GT=0) are excluded because they are not present in the GT file.

    Args:
        abund : pd.DataFrame -- predicted abundance (transcript_id, tpm).
        gt    : pd.DataFrame -- ground truth (transcript_id, tpm).

    Returns:
        pd.DataFrame with columns [transcript_id, tpm_gt, tpm_pred, label]
        where label is "TP" (pred>0) or "FN" (pred=0).
    """
    merged = gt.rename(columns={"tpm": "tpm_gt"}).merge(
        abund[["transcript_id","tpm"]].rename(columns={"tpm":"tpm_pred"}),
        on="transcript_id", how="left"
    )
    merged["tpm_pred"] = merged["tpm_pred"].fillna(0.0)
    merged["label"] = np.where(merged["tpm_pred"] > 0, "TP", "FN")
    return merged


def build_merged(abund: pd.DataFrame, gt: pd.DataFrame,
                 active_only: bool = False) -> pd.DataFrame:
    """
    Outer-join predicted and GT on transcript_id, fill NA with 0.

    Args:
        abund       : pd.DataFrame -- predicted abundance (transcript_id, tpm).
        gt          : pd.DataFrame -- ground truth (transcript_id, tpm).
        active_only : bool -- if True, filter abund to non-zero rows first
                             (removes true-negative (0,0) pairs).

    Returns:
        pd.DataFrame with columns [transcript_id, tpm_pred, tpm_gt].
    """
    if active_only:
        abund = abund[abund["tpm"] > 0]
    merged = abund.merge(gt, on="transcript_id", how="outer", suffixes=("_pred","_gt"))
    merged["tpm_pred"] = merged["tpm_pred"].fillna(0.0)
    merged["tpm_gt"]   = merged["tpm_gt"].fillna(0.0)
    return merged


def compute_all_metrics(abund: pd.DataFrame, gt: pd.DataFrame) -> dict:
    """
    Compute Spearman for three universes: all, active, GT-nonzero.

    Args:
        abund : pd.DataFrame -- predicted (transcript_id, tpm).
        gt    : pd.DataFrame -- ground truth (transcript_id, tpm).

    Returns:
        dict with keys all_spearman, active_spearman, gt_spearman,
        plus universe sizes n_all, n_active, n_gt.
    """
    m_all    = build_merged(abund, gt, active_only=False)
    m_active = build_merged(abund, gt, active_only=True)

    gt_mask = m_all["tpm_gt"] > 0

    return {
        "all_spearman":    spearman(m_all["tpm_pred"].values, m_all["tpm_gt"].values),
        "active_spearman": spearman(m_active["tpm_pred"].values, m_active["tpm_gt"].values),
        "gt_spearman":     spearman(m_all.loc[gt_mask,"tpm_pred"].values,
                                    m_all.loc[gt_mask,"tpm_gt"].values),
        "n_all":    len(m_all),
        "n_active": len(m_active),
        "n_gt":     int(gt_mask.sum()),
    }


def classify(abund: pd.DataFrame, gt: pd.DataFrame) -> pd.DataFrame:
    """
    Classify each transcript in the outer join as TP / FP / FN / TN.

    Args:
        abund : pd.DataFrame -- predicted abundance (transcript_id, tpm).
        gt    : pd.DataFrame -- ground truth (transcript_id, tpm).

    Returns:
        pd.DataFrame with columns [transcript_id, tpm_pred, tpm_gt, label]
        where label ∈ {"TP","FP","FN","TN"}.
    """
    merged = build_merged(abund, gt, active_only=False)
    pred_nz = merged["tpm_pred"] > 0
    gt_nz   = merged["tpm_gt"]   > 0
    conditions = [
        pred_nz &  gt_nz,    # TP
         pred_nz & ~gt_nz,   # FP
        ~pred_nz &  gt_nz,   # FN
        ~pred_nz & ~gt_nz,   # TN
    ]
    merged["label"] = np.select(conditions, ["TP","FP","FN","TN"], default="TN")
    return merged


# ============================================================
# Priority 1 — Index size discrepancy
# ============================================================

def run_priority_1(samples, lk_dir, jk_dir, gt_map, figures_dir, out_dir, log):
    """
    Priority 1: compare transcript universes and metric sensitivity across methods.

    Figures produced:
      fig_metric_sensitivity_{sample}.png  — Spearman by universe for LK vs JK MS
      fig_universe_overlap_{sample}.png    — bar chart of TP/FP/FN/TN counts

    Args:
        samples     : list[str] -- sample names.
        lk_dir      : Path      -- LK experiment directory.
        jk_dir      : Path      -- JK MS experiment directory.
        gt_map      : dict      -- {sample: gt_path}.
        figures_dir : Path      -- output figures directory.
        out_dir     : Path      -- output text directory.
        log         : file      -- open log file handle.
    """
    def _w(msg):
        print(msg)
        log.write(msg + "\n")

    _w("\n" + "="*70)
    _w("PRIORITY 1 — Index size discrepancy")
    _w("="*70)

    table_rows = []

    for sample in samples:
        _w(f"\n--- {sample} ---")

        abund_lk = read_abundance(lk_dir / sample / ABUNDANCE_FILE)
        abund_jk = read_abundance(jk_dir / sample / ABUNDANCE_FILE)
        gt       = read_gt(Path(gt_map[sample]))

        # Contingency tables
        cls_lk = classify(abund_lk, gt)
        cls_jk = classify(abund_jk, gt)

        def _counts(cls):
            return {l: int((cls["label"] == l).sum()) for l in ["TP","FP","FN","TN"]}

        lk_c = _counts(cls_lk)
        jk_c = _counts(cls_jk)

        _w(f"\n  {'Pair type':<40} {'LK':>10} {'JK MS':>10}")
        _w(f"  {'-'*60}")
        for label, desc in [("TP","True  positives (pred>0, GT>0)"),
                             ("FP","False positives (pred>0, GT=0)"),
                             ("FN","False negatives (pred=0, GT>0)"),
                             ("TN","True  negatives (pred=0, GT=0)")]:
            _w(f"  {desc:<40} {lk_c[label]:>10,} {jk_c[label]:>10,}")
        _w(f"  {'Total (outer join)':<40} {sum(lk_c.values()):>10,} {sum(jk_c.values()):>10,}")

        # Metric sensitivity
        m_lk = compute_all_metrics(abund_lk, gt)
        m_jk = compute_all_metrics(abund_jk, gt)

        _w(f"\n  {'Universe':<45} {'LK Spearman':>12} {'JK MS Spearman':>15}")
        _w(f"  {'-'*72}")
        for key, label in [("all_spearman",    f"All transcripts (N_lk={m_lk['n_all']:,} / N_jk={m_jk['n_all']:,})"),
                            ("active_spearman", f"Active: GT∪nonzero pred (N_lk={m_lk['n_active']:,} / N_jk={m_jk['n_active']:,})"),
                            ("gt_spearman",     f"GT all only (N={m_lk['n_gt']:,})")]:
            _w(f"  {label:<45} {m_lk[key]:>12.4f} {m_jk[key]:>15.4f}")

        table_rows.append((sample, lk_c, jk_c, m_lk, m_jk))

        # ── Fig P1: Metric sensitivity + Universe overlap (1×2 subplots) ──
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
        fig.suptitle(f"Priority 1 — Index size discrepancy | {sample}", fontsize=12, fontweight="bold")

        # Left: Spearman by universe definition
        universe_labels = ["All\ntranscripts", "Active universe\n(GT∪nonzero pred)", "GT all\nonly"]
        lk_vals = [m_lk["all_spearman"], m_lk["active_spearman"], m_lk["gt_spearman"]]
        jk_vals = [m_jk["all_spearman"], m_jk["active_spearman"], m_jk["gt_spearman"]]
        x = np.arange(len(universe_labels))
        width = 0.35
        bars_lk = ax1.bar(x - width/2, lk_vals, width, label="LK",    color=C_LK, alpha=0.85)
        bars_jk = ax1.bar(x + width/2, jk_vals, width, label="JK MS", color=C_JK, alpha=0.85)
        for bar in list(bars_lk) + list(bars_jk):
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.003,
                     f"{bar.get_height():.3f}", ha="center", va="bottom", fontsize=8)
        ax1.set_xticks(x)
        ax1.set_xticklabels(universe_labels, fontsize=9)
        ax1.set_ylabel("Spearman correlation")
        ax1.set_ylim(0, 1.08)
        ax1.set_title("Spearman by universe definition")
        ax1.legend(fontsize=9)
        ax1.grid(axis="y", alpha=0.3)

        # Right: TP/FP/FN/TN counts for LK vs JK MS
        categories = ["TP\n(pred>0,GT>0)", "FP\n(pred>0,GT=0)", "FN\n(pred=0,GT>0)", "TN\n(pred=0,GT=0)"]
        lk_vals2 = [lk_c["TP"], lk_c["FP"], lk_c["FN"], lk_c["TN"]]
        jk_vals2 = [jk_c["TP"], jk_c["FP"], jk_c["FN"], jk_c["TN"]]
        x2 = np.arange(len(categories))
        bars_lk2 = ax2.bar(x2 - width/2, lk_vals2, width, label="LK",    color=C_LK, alpha=0.85)
        bars_jk2 = ax2.bar(x2 + width/2, jk_vals2, width, label="JK MS", color=C_JK, alpha=0.85)
        for bar in list(bars_lk2) + list(bars_jk2):
            h = bar.get_height()
            if h > 0:
                ax2.text(bar.get_x() + bar.get_width()/2, h + 200,
                         f"{h:,}", ha="center", va="bottom", fontsize=7, rotation=45)
        ax2.set_xticks(x2)
        ax2.set_xticklabels(categories, fontsize=9)
        ax2.set_ylabel("Number of transcripts")
        ax2.set_title("Transcript classification (LK vs JK MS)")
        ax2.legend(fontsize=9)
        ax2.grid(axis="y", alpha=0.3)

        fig.tight_layout()
        out_path = figures_dir / f"fig_p1_{sample}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        _w(f"\n  [P1] Saved: {out_path.name}")

    # Save table to file
    table_path = out_dir / "priority1_universe_table.txt"
    with open(table_path, "w") as fh:
        fh.write("Priority 1 — Universe comparison table\n")
        fh.write(f"LK:    {LK_EXP}\n")
        fh.write(f"JK MS: {JK_EXP}\n\n")
        for sample, lk_c, jk_c, m_lk, m_jk in table_rows:
            fh.write(f"\n{'='*70}\n{sample}\n{'='*70}\n")
            fh.write(f"  {'Pair type':<40} {'LK':>10} {'JK MS':>10}\n")
            for label, desc in [("TP","True  positives (pred>0, GT>0)"),
                                 ("FP","False positives (pred>0, GT=0)"),
                                 ("FN","False negatives (pred=0, GT>0)"),
                                 ("TN","True  negatives (pred=0, GT=0)")]:
                fh.write(f"  {desc:<40} {lk_c[label]:>10,} {jk_c[label]:>10,}\n")
            fh.write(f"\n  Spearman by universe:\n")
            fh.write(f"  {'Universe':<45} {'LK':>10} {'JK MS':>10}\n")
            for key, label in [("all_spearman","All"), ("active_spearman","Active"), ("gt_spearman","GT-nonzero")]:
                fh.write(f"  {label:<45} {m_lk[key]:>10.4f} {m_jk[key]:>10.4f}\n")
    _w(f"\n[P1] Saved table: {table_path}")


# ============================================================
# Priority 6 — Leaked transcript analysis
# ============================================================

def run_priority_6(samples, lk_dir, jk_dir, gt_map, bustools_dirs,
                   figures_dir, out_dir, log):
    """
    Priority 6: characterize false positive (leaked) transcripts in BOTH LK and JK MS.

    Figures produced per sample (all show LK and JK MS side-by-side):
      fig_p6a_tpm_{sample}.png       — 2×3: [LK row | JK MS row] × [TP hist | FP hist | FP CDF]
      fig_p6b_rank_lk_{sample}.png   — 3×3 gridspec: LK rank (col 0) | JK rank (col 1) |
                                         LK vs JK FP scatter (col 2, full height)
      fig_p6c_ec_sizes_{sample}.png  — 2×2: LK TP EC | LK FP EC (row 0)
                                             JK TP EC | JK FP EC (row 1), shared x
                                       (only if bustools dir is set)
      fig_p6d_theta_scatter_{sample}.png — 1×2: LK theta scatter | JK MS theta scatter

    Args:
        samples      : list[str]      -- sample names.
        lk_dir       : Path           -- LK experiment directory.
        jk_dir       : Path           -- JK MS experiment directory.
        gt_map       : dict           -- {sample: gt_path}.
        bustools_dirs: dict           -- {sample: bustools_output_dir or None}.
        figures_dir  : Path           -- output figures directory.
        out_dir      : Path           -- output text directory.
        log          : file           -- open log file handle.
    """
    def _w(msg):
        print(msg)
        log.write(msg + "\n")

    _w("\n" + "="*70)
    _w("PRIORITY 6 — Leaked transcript analysis (LK and JK MS)")
    _w("="*70)

    stats_rows = []

    for sample in samples:
        _w(f"\n--- {sample} ---")

        abund_lk = read_abundance(lk_dir / sample / ABUNDANCE_FILE)
        abund_jk = read_abundance(jk_dir / sample / ABUNDANCE_FILE)
        gt       = read_gt(Path(gt_map[sample]))

        # Classify BOTH LK and JK MS transcripts relative to GT
        cls_lk = classify(abund_lk, gt)
        cls_jk = classify(abund_jk, gt)

        fp_lk = cls_lk[cls_lk["label"] == "FP"]
        tp_lk = cls_lk[cls_lk["label"] == "TP"]
        fn_lk = cls_lk[cls_lk["label"] == "FN"]

        fp_jk = cls_jk[cls_jk["label"] == "FP"]
        tp_jk = cls_jk[cls_jk["label"] == "TP"]
        fn_jk = cls_jk[cls_jk["label"] == "FN"]

        # TPM budget consumed by FP transcripts — both methods
        total_lk      = cls_lk["tpm_pred"].sum()
        total_jk      = cls_jk["tpm_pred"].sum()
        leaked_lk     = fp_lk["tpm_pred"].sum()
        leaked_jk     = fp_jk["tpm_pred"].sum()
        leaked_lk_pct = 100.0 * leaked_lk / total_lk if total_lk > 0 else 0.0
        leaked_jk_pct = 100.0 * leaked_jk / total_jk if total_jk > 0 else 0.0

        _w(f"\n  Transcript counts ({sample}):")
        _w(f"  {'':35} {'LK':>10} {'JK MS':>10}")
        _w(f"  {'TP (pred>0, GT>0)':<35} {len(tp_lk):>10,} {len(tp_jk):>10,}")
        _w(f"  {'FP (pred>0, GT=0)  ← leaked':<35} {len(fp_lk):>10,} {len(fp_jk):>10,}")
        _w(f"  {'FN (pred=0, GT>0)':<35} {len(fn_lk):>10,} {len(fn_jk):>10,}")
        _w(f"  {'Leaked TPM total':<35} {leaked_lk:>10.2f} {leaked_jk:>10.2f}")
        _w(f"  {'Leaked TPM %':<35} {leaked_lk_pct:>9.2f}% {leaked_jk_pct:>9.2f}%")

        stats_rows.append({
            "sample": sample,
            "lk_tp": len(tp_lk), "lk_fp": len(fp_lk), "lk_fn": len(fn_lk),
            "jk_tp": len(tp_jk), "jk_fp": len(fp_jk), "jk_fn": len(fn_jk),
            "lk_leaked_tpm": leaked_lk, "lk_total_tpm": total_lk, "lk_leaked_pct": leaked_lk_pct,
            "jk_leaked_tpm": leaked_jk, "jk_total_tpm": total_jk, "jk_leaked_pct": leaked_jk_pct,
        })

        # JK MS FP transcripts cross-referenced with LK predictions (for scatter)
        fp_with_lk = fp_jk[["transcript_id","tpm_pred"]].rename(columns={"tpm_pred":"jk_tpm"})
        fp_with_lk = fp_with_lk.merge(
            abund_lk[["transcript_id","tpm"]].rename(columns={"tpm":"lk_tpm"}),
            on="transcript_id", how="left"
        )
        fp_with_lk["lk_tpm"] = fp_with_lk["lk_tpm"].fillna(0.0)

        lk_also_nz = int((fp_with_lk["lk_tpm"] > 0).sum())
        lk_is_zero = int((fp_with_lk["lk_tpm"] == 0).sum())
        _w(f"\n  Among {len(fp_jk):,} JK MS FP transcripts:")
        _w(f"    LK also predicts non-zero : {lk_also_nz:>8,}  (leakage in both)")
        _w(f"    LK predicts zero          : {lk_is_zero:>8,}  (JK-only, stronger prior effect)")

        eps = 1e-7
        from matplotlib.patches import Patch

        tp_lk_tpm = tp_lk["tpm_pred"].values
        fp_lk_tpm = fp_lk["tpm_pred"].values
        tp_jk_tpm = tp_jk["tpm_pred"].values
        fp_jk_tpm = fp_jk["tpm_pred"].values

        # ── Fig P6a: TPM distribution — 2×3 subplots ─────────────────
        # Row 0: LK  — [TP hist | FP hist | FP CDF]
        # Row 1: JK MS — [TP hist | FP hist | FP CDF]
        shared_max = max(tp_lk_tpm.max(), fp_lk_tpm.max(),
                         tp_jk_tpm.max(), fp_jk_tpm.max()) * 2
        bins = np.logspace(np.log10(eps), np.log10(shared_max), 60)

        fp_lk_sorted = np.sort(fp_lk_tpm)
        cdf_lk = np.arange(1, len(fp_lk_sorted) + 1) / len(fp_lk_sorted)
        fp_jk_sorted = np.sort(fp_jk_tpm)
        cdf_jk = np.arange(1, len(fp_jk_sorted) + 1) / len(fp_jk_sorted)

        # ── Fig P6a (density): 2×2 — TP density | FP density per row ─
        fig, axes = plt.subplots(2, 2, figsize=(14, 9))
        fig.suptitle(f"P6 — TPM distribution (density): LK vs JK MS | {sample}",
                     fontsize=12, fontweight="bold")

        for row, method, c_method, tp_tpm, fp_tpm in [
            (0, "LK",    C_LK, tp_lk_tpm, fp_lk_tpm),
            (1, "JK MS", C_JK, tp_jk_tpm, fp_jk_tpm),
        ]:
            axes[row, 0].hist(tp_tpm, bins=bins, color=C_TP, density=True)
            axes[row, 0].set_xscale("log")
            axes[row, 0].set_ylabel(f"{method}\nDensity")
            axes[row, 0].set_title(f"True positives — {method}  (N={len(tp_tpm):,})")
            axes[row, 0].grid(alpha=0.3)
            axes[row, 0].set_xlim(eps, shared_max)

            axes[row, 1].hist(fp_tpm, bins=bins, color=c_method, density=True)
            axes[row, 1].set_xscale("log")
            axes[row, 1].set_title(f"False positives — {method}  (N={len(fp_tpm):,})")
            axes[row, 1].grid(alpha=0.3)
            axes[row, 1].set_xlim(eps, shared_max)

        for col in range(2):
            axes[1, col].set_xlabel("Predicted TPM")

        fig.tight_layout()
        out_path = figures_dir / f"fig_p6a_tpm_{sample}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        _w(f"\n  [P6] Saved: {out_path.name}")

        # ── Fig P6a (counts + CDF): 2×3 — TP count | FP count | FP CDF ──
        # Row 0: LK  — [TP count | FP count | FP CDF]
        # Row 1: JK MS — [TP count | FP count | FP CDF]
        fig, axes_cnt = plt.subplots(2, 3, figsize=(20, 9))
        fig.suptitle(f"P6 — TPM distribution (counts + CDF): LK vs JK MS | {sample}",
                     fontsize=12, fontweight="bold")

        for row, method, c_method, tp_tpm, fp_tpm, fp_sorted, cdf in [
            (0, "LK",    C_LK, tp_lk_tpm, fp_lk_tpm, fp_lk_sorted, cdf_lk),
            (1, "JK MS", C_JK, tp_jk_tpm, fp_jk_tpm, fp_jk_sorted, cdf_jk),
        ]:
            axes_cnt[row, 0].hist(tp_tpm, bins=bins, color=C_TP)
            axes_cnt[row, 0].set_xscale("log")
            axes_cnt[row, 0].set_ylabel(f"{method}\nCount")
            axes_cnt[row, 0].set_title(f"True positives — {method}  (N={len(tp_tpm):,})")
            axes_cnt[row, 0].grid(alpha=0.3)
            axes_cnt[row, 0].set_xlim(eps, shared_max)

            axes_cnt[row, 1].hist(fp_tpm, bins=bins, color=c_method)
            axes_cnt[row, 1].set_xscale("log")
            axes_cnt[row, 1].set_title(f"False positives — {method}  (N={len(fp_tpm):,})")
            axes_cnt[row, 1].grid(alpha=0.3)
            axes_cnt[row, 1].set_xlim(eps, shared_max)

            axes_cnt[row, 2].plot(fp_sorted, cdf, color=c_method, linewidth=2)
            for thr, ls in [(1.0,"--"), (5.0,":"), (10.0,"-.")]:
                frac = float(np.mean(fp_sorted <= thr))
                axes_cnt[row, 2].axvline(thr, color="gray", linestyle=ls, alpha=0.7,
                                         label=f"TPM≤{thr}: {frac*100:.1f}%")
            axes_cnt[row, 2].set_xscale("log")
            axes_cnt[row, 2].set_title(f"CDF — FP transcripts — {method}")
            axes_cnt[row, 2].legend(fontsize=7)
            axes_cnt[row, 2].grid(alpha=0.3)

        for col, lbl in enumerate(["Predicted TPM", "Predicted TPM",
                                    "Predicted TPM threshold"]):
            axes_cnt[1, col].set_xlabel(lbl)

        fig.tight_layout()
        out_path = figures_dir / f"fig_p6a_tpm_counts_{sample}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        _w(f"  [P6] Saved: {out_path.name}")

        # ── Fig P6b: Rank plots — 3×3 gridspec ───────────────────────
        # Col 0: LK  rank panels (TP / FP / FN stacked, shared axes)
        # Col 1: JK MS rank panels (TP / FP / FN stacked, shared axes)
        # Col 2: LK vs JK MS scatter for JK MS FP transcripts (full height)
        nz_lk = cls_lk[cls_lk["tpm_pred"] > 0].copy()
        nz_lk = nz_lk.sort_values("tpm_pred", ascending=False).reset_index(drop=True)
        nz_lk["rank"] = np.arange(1, len(nz_lk) + 1)

        nz_jk = cls_jk[cls_jk["tpm_pred"] > 0].copy()
        nz_jk = nz_jk.sort_values("tpm_pred", ascending=False).reset_index(drop=True)
        nz_jk["rank"] = np.arange(1, len(nz_jk) + 1)

        jk_fp_arr = fp_with_lk["jk_tpm"].values
        lk_fp_arr = fp_with_lk["lk_tpm"].values
        jk_plot   = np.where(jk_fp_arr > 0, jk_fp_arr, eps)
        lk_plot   = np.where(lk_fp_arr > 0, lk_fp_arr, eps)
        colors_pt = np.where(lk_fp_arr > 0, "#e66101", "#5e3c99")

        fig = plt.figure(figsize=(22, 12))
        fig.suptitle(f"P6 — Rank plots: LK vs JK MS | {sample}",
                     fontsize=12, fontweight="bold")
        gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)

        # LK rank panels (col 0)
        ax_lk_tp = fig.add_subplot(gs[0, 0])
        ax_lk_fp = fig.add_subplot(gs[1, 0], sharex=ax_lk_tp, sharey=ax_lk_tp)
        ax_lk_fn = fig.add_subplot(gs[2, 0], sharex=ax_lk_tp, sharey=ax_lk_tp)
        lk_rank_max = nz_lk["rank"].max()
        for ax, lbl, col, desc in [
            (ax_lk_tp, "TP", C_TP, f"LK TP (N={len(tp_lk):,})"),
            (ax_lk_fp, "FP", C_LK, f"LK FP (N={len(fp_lk):,})"),
            (ax_lk_fn, "FN", C_FN, f"LK FN in non-zero (N={int((nz_lk['label']=='FN').sum())})"),
        ]:
            sub = nz_lk[nz_lk["label"] == lbl]
            ax.scatter(sub["rank"], sub["tpm_pred"], c=col, s=2, alpha=0.5, rasterized=True)
            ax.set_yscale("log")
            ax.set_xlim(0, lk_rank_max)
            ax.set_title(desc, fontsize=9)
            ax.set_ylabel("Predicted TPM")
            ax.grid(alpha=0.3)
        ax_lk_fn.set_xlabel("Rank (LK)")
        plt.setp(ax_lk_tp.get_xticklabels(), visible=False)
        plt.setp(ax_lk_fp.get_xticklabels(), visible=False)

        # JK MS rank panels (col 1)
        ax_jk_tp = fig.add_subplot(gs[0, 1])
        ax_jk_fp = fig.add_subplot(gs[1, 1], sharex=ax_jk_tp, sharey=ax_jk_tp)
        ax_jk_fn = fig.add_subplot(gs[2, 1], sharex=ax_jk_tp, sharey=ax_jk_tp)
        jk_rank_max = nz_jk["rank"].max()
        for ax, lbl, col, desc in [
            (ax_jk_tp, "TP", C_TP, f"JK MS TP (N={len(tp_jk):,})"),
            (ax_jk_fp, "FP", C_FP, f"JK MS FP (N={len(fp_jk):,})"),
            (ax_jk_fn, "FN", C_FN, f"JK MS FN in non-zero (N={int((nz_jk['label']=='FN').sum())})"),
        ]:
            sub = nz_jk[nz_jk["label"] == lbl]
            ax.scatter(sub["rank"], sub["tpm_pred"], c=col, s=2, alpha=0.5, rasterized=True)
            ax.set_yscale("log")
            ax.set_xlim(0, jk_rank_max)
            ax.set_title(desc, fontsize=9)
            ax.set_ylabel("Predicted TPM")
            ax.grid(alpha=0.3)
        ax_jk_fn.set_xlabel("Rank (JK MS)")
        plt.setp(ax_jk_tp.get_xticklabels(), visible=False)
        plt.setp(ax_jk_fp.get_xticklabels(), visible=False)

        # Scatter: JK MS FP — LK vs JK MS predicted TPM (col 2, full height)
        ax_scat = fig.add_subplot(gs[:, 2])
        ax_scat.scatter(lk_plot, jk_plot, c=colors_pt, s=3, alpha=0.4, rasterized=True)
        ax_scat.axvline(eps, color="gray", linestyle=":", alpha=0.5, linewidth=0.8)
        ax_scat.set_xscale("log")
        ax_scat.set_yscale("log")
        ax_scat.set_xlabel("LK predicted TPM (0 shown at 1e-7)")
        ax_scat.set_ylabel("JK MS predicted TPM")
        ax_scat.set_title(f"JK MS FP transcripts\nLK vs JK MS predicted TPM")
        ax_scat.legend(handles=[
            Patch(color="#e66101", label=f"LK also non-zero ({lk_also_nz:,})"),
            Patch(color="#5e3c99", label=f"LK = 0 ({lk_is_zero:,})"),
        ], fontsize=9)
        ax_scat.grid(alpha=0.3)

        out_path = figures_dir / f"fig_p6b_rank_lk_{sample}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        _w(f"  [P6] Saved: {out_path.name}")

        # ── Fig P6c: EC size distribution — 2×2 subplots ─────────────
        # Row 0: LK  — [TP EC sizes | FP EC sizes]
        # Row 1: JK MS — [TP EC sizes | FP EC sizes]   shared x-axis
        bdir = bustools_dirs.get(sample)
        if bdir is not None:
            _w(f"\n  [P6] Loading EC data for EC size figure from: {bdir}")
            ec_map    = read_matrix_ec(bdir)
            tx_names  = read_transcripts_txt(bdir)
            tx_to_idx = {name: idx for idx, name in enumerate(tx_names)}

            tx_ec_sizes: dict = {i: [] for i in range(len(tx_names))}
            for _, tx_list in ec_map.items():
                ec_size = len(tx_list)
                for tx_idx in tx_list:
                    tx_ec_sizes[tx_idx].append(ec_size)

            def get_min_ec_size(transcript_id: str) -> int:
                """Return the minimum EC size this transcript appears in (0 if not found)."""
                idx = tx_to_idx.get(transcript_id, None)
                if idx is None or not tx_ec_sizes[idx]:
                    return 0
                return min(tx_ec_sizes[idx])

            tp_lk_ec = np.array([get_min_ec_size(t) for t in tp_lk["transcript_id"]])
            fp_lk_ec = np.array([get_min_ec_size(t) for t in fp_lk["transcript_id"]])
            tp_jk_ec = np.array([get_min_ec_size(t) for t in tp_jk["transcript_id"]])
            fp_jk_ec = np.array([get_min_ec_size(t) for t in fp_jk["transcript_id"]])

            tp_lk_ec_nz = tp_lk_ec[tp_lk_ec > 0]
            fp_lk_ec_nz = fp_lk_ec[fp_lk_ec > 0]
            tp_jk_ec_nz = tp_jk_ec[tp_jk_ec > 0]
            fp_jk_ec_nz = fp_jk_ec[fp_jk_ec > 0]

            all_nz  = np.concatenate([tp_lk_ec_nz, fp_lk_ec_nz, tp_jk_ec_nz, fp_jk_ec_nz])
            bins_ec = np.arange(0.5, min(all_nz.max(), 50) + 1.5, 1)

            # 2×4: cols 0-1 = density (TP | FP), cols 2-3 = count (TP | FP), rows = LK | JK MS
            fig, axes_ec = plt.subplots(2, 4, figsize=(22, 9), sharex=True)
            fig.suptitle(f"P6 — EC size distribution: LK vs JK MS | {sample}",
                         fontsize=12, fontweight="bold")

            ec_panels = [
                # (row, col_density, col_count, data, color, method_label)
                (0, tp_lk_ec_nz, C_TP, fp_lk_ec_nz, C_LK, "LK"),
                (1, tp_jk_ec_nz, C_TP, fp_jk_ec_nz, C_FP, "JK MS"),
            ]
            for row, tp_data, tp_col, fp_data, fp_col, method in ec_panels:
                # Density columns (0, 1)
                axes_ec[row, 0].hist(tp_data, bins=bins_ec, color=tp_col, density=True)
                axes_ec[row, 0].set_title(f"{method} — TP density  (N={len(tp_data):,})", fontsize=9)
                axes_ec[row, 0].set_ylabel(f"{method}\nDensity")
                axes_ec[row, 0].grid(alpha=0.3)

                axes_ec[row, 1].hist(fp_data, bins=bins_ec, color=fp_col, density=True)
                axes_ec[row, 1].set_title(f"{method} — FP density  (N={len(fp_data):,})", fontsize=9)
                axes_ec[row, 1].grid(alpha=0.3)

                # Count columns (2, 3)
                axes_ec[row, 2].hist(tp_data, bins=bins_ec, color=tp_col)
                axes_ec[row, 2].set_title(f"{method} — TP count  (N={len(tp_data):,})", fontsize=9)
                axes_ec[row, 2].set_ylabel(f"{method}\nCount")
                axes_ec[row, 2].grid(alpha=0.3)

                axes_ec[row, 3].hist(fp_data, bins=bins_ec, color=fp_col)
                axes_ec[row, 3].set_title(f"{method} — FP count  (N={len(fp_data):,})", fontsize=9)
                axes_ec[row, 3].grid(alpha=0.3)

            x_label = "Min EC size transcript appears in\n(1 = single-tx EC; larger = more multi-mapping)"
            for col in range(4):
                axes_ec[1, col].set_xlabel(x_label)

            fig.tight_layout()
            out_path = figures_dir / f"fig_p6c_ec_sizes_{sample}.png"
            fig.savefig(out_path, dpi=150)
            plt.close(fig)
            _w(f"  [P6] Saved: {out_path.name}")

            _w(f"\n  EC size breakdown for FP transcripts:")
            _w(f"  {'Category':<25} {'LK FP':>12} {'JK FP':>12}")
            for size_label, lk_mask, jk_mask in [
                ("Single-tx EC (size=1)",
                 fp_lk_ec_nz == 1,   fp_jk_ec_nz == 1),
                ("Size 2",
                 fp_lk_ec_nz == 2,   fp_jk_ec_nz == 2),
                ("Size 3-5",
                 (fp_lk_ec_nz >= 3) & (fp_lk_ec_nz <= 5),
                 (fp_jk_ec_nz >= 3) & (fp_jk_ec_nz <= 5)),
                ("Size 6-10",
                 (fp_lk_ec_nz >= 6) & (fp_lk_ec_nz <= 10),
                 (fp_jk_ec_nz >= 6) & (fp_jk_ec_nz <= 10)),
                ("Size >10",
                 fp_lk_ec_nz > 10,   fp_jk_ec_nz > 10),
            ]:
                _w(f"  {size_label:<25} {lk_mask.sum():>6,} ({100*lk_mask.mean():.1f}%)"
                   f"  {jk_mask.sum():>6,} ({100*jk_mask.mean():.1f}%)")
        else:
            _w(f"  [P6] Skipping EC size figure (bustools dir not set for {sample})")

        # ── Fig P6d: Theta scatter — 1×2: LK (left) | JK MS (right) ─
        def _theta_panel(ax, abund, method_label, c_fp_dot):
            """Draw GT theta vs predicted theta scatter on a single axis."""
            active = build_merged(abund, gt, active_only=True)
            total_pred_sum = active["tpm_pred"].sum()
            total_gt_sum   = active["tpm_gt"].sum()
            active = active.copy()
            active["theta_pred"] = (active["tpm_pred"] / total_pred_sum
                                    if total_pred_sum > 0 else 0.0)
            active["theta_gt"]   = (active["tpm_gt"]   / total_gt_sum
                                    if total_gt_sum > 0 else 0.0)

            color_val = active["theta_gt"].values
            x_p = np.where(active["theta_gt"]   > 0, active["theta_gt"],   eps)
            y_p = np.where(active["theta_pred"] > 0, active["theta_pred"], eps)

            is_fp_t    = active["theta_gt"] == 0
            gt_nz_min  = active.loc[~is_fp_t, "theta_gt"].min()
            gt_max     = active["theta_gt"].max()
            lnorm      = LogNorm(vmin=gt_nz_min, vmax=gt_max)

            ax.scatter(x_p[is_fp_t],  y_p[is_fp_t],
                       c=c_fp_dot, s=1, alpha=0.3, rasterized=True,
                       label=f"FP (GT=0, N={is_fp_t.sum():,})")
            sc = ax.scatter(x_p[~is_fp_t], y_p[~is_fp_t],
                            c=color_val[~is_fp_t], cmap=WHITE_TO_RED, norm=lnorm,
                            s=2, alpha=0.6, rasterized=True,
                            label=f"GT>0 (N={(~is_fp_t).sum():,})")
            lims = [min(x_p.min(), y_p.min()), max(x_p.max(), y_p.max())]
            ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5, label="y=x")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("GT theta (log scale)")
            ax.set_ylabel(f"{method_label} predicted theta (log scale)")
            sp_val = spearman(active.loc[~is_fp_t, "theta_gt"].values,
                              active.loc[~is_fp_t, "theta_pred"].values)
            ax.set_title(f"{method_label} | {sample}\n"
                         f"Color = GT theta  (red=high, white=near zero)\n"
                         f"Spearman (GT>0) = {sp_val:.4f}")
            ax.legend(markerscale=5, fontsize=8)
            ax.grid(alpha=0.3)
            return sc

        fig, (ax_lk_sc, ax_jk_sc) = plt.subplots(1, 2, figsize=(16, 7))
        sc_lk = _theta_panel(ax_lk_sc, abund_lk, "LK",    c_fp_dot="#aaaaaa")
        sc_jk = _theta_panel(ax_jk_sc, abund_jk, "JK MS", c_fp_dot="#dddddd")
        fig.colorbar(sc_lk, ax=ax_lk_sc, label="GT theta (log scale)")
        fig.colorbar(sc_jk, ax=ax_jk_sc, label="GT theta (log scale)")
        fig.suptitle(f"P6 — GT theta vs predicted theta | {sample}",
                     fontsize=12, fontweight="bold")
        fig.tight_layout()
        out_path = figures_dir / f"fig_p6d_theta_scatter_{sample}.png"
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        _w(f"  [P6] Saved: {out_path.name}")

    # Save stats table (both methods)
    stats_path = out_dir / "priority6_leaked_stats.txt"
    with open(stats_path, "w") as fh:
        fh.write("Priority 6 — Leaked transcript statistics (LK and JK MS)\n")
        fh.write(f"LK:    {LK_EXP}\n")
        fh.write(f"JK MS: {JK_EXP}\n\n")
        for row in stats_rows:
            fh.write(f"\n{'='*60}\n{row['sample']}\n{'='*60}\n")
            fh.write(f"  {'':35} {'LK':>10} {'JK MS':>10}\n")
            fh.write(f"  {'TP':<35} {row['lk_tp']:>10,} {row['jk_tp']:>10,}\n")
            fh.write(f"  {'FP (leaked)':<35} {row['lk_fp']:>10,} {row['jk_fp']:>10,}\n")
            fh.write(f"  {'FN':<35} {row['lk_fn']:>10,} {row['jk_fn']:>10,}\n")
            fh.write(f"  {'Leaked TPM total':<35} {row['lk_leaked_tpm']:>10.2f}"
                     f" {row['jk_leaked_tpm']:>10.2f}\n")
            fh.write(f"  {'Leaked TPM %':<35} {row['lk_leaked_pct']:>9.2f}%"
                     f" {row['jk_leaked_pct']:>9.2f}%\n")
    _w(f"\n[P6] Saved stats: {stats_path}")


# ============================================================
# Main figures — GT vs predicted TPM
# ============================================================

def run_main_figures(samples, lk_dir, jk_dir, gt_map, figures_dir, log):
    """
    Produce the two primary GT vs predicted TPM scatter figures saved to figures/.

    Both figures are 1×2 panels (LK left | JK MS right) so both methods are
    always visible.  Only GT transcripts are shown (FP transcripts, which have
    no GT entry, are excluded).

    Fig 1 — plain scatter (fig_gt_vs_pred_{sample}.png):
      x = GT TPM, y = predicted TPM (log-log).  FN transcripts shown at eps
      on the y-axis so they are visible on the log scale.

    Fig 2 — classified scatter (fig_gt_vs_pred_classified_{sample}.png):
      Same axes, but points colored TP (green, pred>0) or FN (orange, pred=0).

    Args:
        samples     : list[str] -- sample names.
        lk_dir      : Path      -- LK experiment directory.
        jk_dir      : Path      -- JK MS experiment directory.
        gt_map      : dict      -- {sample: gt_path}.
        figures_dir : Path      -- output figures directory (figures/, NOT extra_fig/).
        log         : file      -- open log file handle.
    """
    def _w(msg):
        print(msg)
        log.write(msg + "\n")

    _w("\n" + "="*70)
    _w("MAIN FIGURES — GT all vs predicted TPM")
    _w("="*70)

    eps = 1e-7   # floor for log-scale (FN transcripts with pred=0)

    for sample in samples:
        _w(f"\n--- {sample} ---")

        abund_lk = read_abundance(lk_dir / sample / ABUNDANCE_FILE)
        abund_jk = read_abundance(jk_dir / sample / ABUNDANCE_FILE)
        gt       = read_gt(Path(gt_map[sample]))

        # GT all: all GT transcripts, left-joined to predictions
        gt_lk = build_gt_all(abund_lk, gt)
        gt_jk = build_gt_all(abund_jk, gt)

        _w(f"  GT all size           : {len(gt):,}")
        _w(f"  LK  TP={int((gt_lk['label']=='TP').sum()):,}  "
           f"FN={int((gt_lk['label']=='FN').sum()):,}")
        _w(f"  JK  TP={int((gt_jk['label']=='TP').sum()):,}  "
           f"FN={int((gt_jk['label']=='FN').sum()):,}")

        def _prep(df):
            """Return (x_gt, y_pred) arrays ready for log-scale plotting."""
            x = df["tpm_gt"].values.astype(float)
            y = df["tpm_pred"].values.astype(float)
            # FN transcripts have pred=0 → place at eps so they show on log scale
            y = np.where(y > 0, y, eps)
            # GT TPM should already be > 0 (GT file only has expressed transcripts)
            x = np.where(x > 0, x, eps)
            return x, y

        x_lk, y_lk = _prep(gt_lk)
        x_jk, y_jk = _prep(gt_jk)

        sp_lk = spearman(gt_lk["tpm_gt"].values, gt_lk["tpm_pred"].values)
        sp_jk = spearman(gt_jk["tpm_gt"].values, gt_jk["tpm_pred"].values)

        # ── Fig 1: plain scatter ──────────────────────────────────────
        fig, (ax_lk, ax_jk) = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(f"GT all vs predicted TPM | {sample}  "
                     f"(FN shown at y=1e-7)",
                     fontsize=11, fontweight="bold")

        for ax, x, y, c, method, sp in [
            (ax_lk, x_lk, y_lk, C_LK, "LK",    sp_lk),
            (ax_jk, x_jk, y_jk, C_JK, "JK MS", sp_jk),
        ]:
            ax.scatter(x, y, c=c, s=2, alpha=0.4, rasterized=True)
            # y=x reference line
            lims = [min(x.min(), y.min()), max(x.max(), y.max())]
            ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5, label="y=x")
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("GT TPM (log scale)")
            ax.set_ylabel("Predicted TPM (log scale)")
            ax.set_title(f"{method}  |  N={len(x):,}  |  Spearman={sp:.4f}")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)

        fig.tight_layout()
        out_path = figures_dir / f"fig_gt_vs_pred_{sample}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        _w(f"  Saved: {out_path.name}")

        # ── Fig 2: classified scatter (TP=green / FN=orange / FP=magenta) ──
        # FP transcripts have GT=0 → shown at x=eps (left edge of log scale).
        # FN transcripts have pred=0 → shown at y=eps (bottom edge).
        cls_lk_full = classify(abund_lk, gt)
        cls_jk_full = classify(abund_jk, gt)

        fig, (ax_lk, ax_jk) = plt.subplots(1, 2, figsize=(14, 6))
        fig.suptitle(f"GT all vs predicted TPM — classified | {sample}\n"
                     f"FN shown at y=1e-7  |  FP shown at x=1e-7",
                     fontsize=11, fontweight="bold")

        from matplotlib.patches import Patch

        for ax, df_gt, cls, method, sp in [
            (ax_lk, gt_lk, cls_lk_full, "LK",    sp_lk),
            (ax_jk, gt_jk, cls_jk_full, "JK MS", sp_jk),
        ]:
            # TP and FN come from GT-centric df (x = GT TPM, y = pred TPM / eps)
            is_tp = (df_gt["label"] == "TP").values
            is_fn = ~is_tp
            x_gt  = np.where(df_gt["tpm_gt"].values   > 0, df_gt["tpm_gt"].values,   eps)
            y_pred = np.where(df_gt["tpm_pred"].values > 0, df_gt["tpm_pred"].values, eps)

            # FP: pred>0, GT=0 → x=eps, y=pred TPM
            fp_rows = cls[cls["label"] == "FP"]
            x_fp = np.full(len(fp_rows), eps)
            y_fp = fp_rows["tpm_pred"].values

            # Draw order: FP bottom, FN middle, TP on top
            ax.scatter(x_fp,        y_fp,        c=C_FP, s=2, alpha=0.3,
                       rasterized=True, zorder=1)
            ax.scatter(x_gt[is_fn], y_pred[is_fn], c=C_FN, s=3, alpha=0.5,
                       rasterized=True, zorder=2)
            ax.scatter(x_gt[is_tp], y_pred[is_tp], c=C_TP, s=2, alpha=0.4,
                       rasterized=True, zorder=3)

            all_x = np.concatenate([x_gt, x_fp])
            all_y = np.concatenate([y_pred, y_fp])
            lims = [min(all_x.min(), all_y.min()), max(all_x.max(), all_y.max())]
            ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("GT TPM (log scale)")
            ax.set_ylabel("Predicted TPM (log scale)")
            n_tp = int(is_tp.sum())
            n_fn = int(is_fn.sum())
            n_fp = len(fp_rows)
            ax.set_title(f"{method}  |  Spearman={sp:.4f}")
            ax.legend(handles=[
                Patch(color=C_TP, label=f"TP ({n_tp:,})"),
                Patch(color=C_FN, label=f"FN ({n_fn:,}) — y=1e-7"),
                Patch(color=C_FP, label=f"FP ({n_fp:,}) — x=1e-7"),
            ], fontsize=8)
            ax.grid(alpha=0.3)

        fig.tight_layout()
        out_path = figures_dir / f"fig_gt_vs_pred_classified_{sample}.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        _w(f"  Saved: {out_path.name}")


# ============================================================
# Main
# ============================================================

def _setup_out_dir(exp_dir: Path, timestamp: str) -> tuple:
    """
    Create eda_{timestamp}/, figures/, and extra_fig/ inside an experiment folder.

    figures/   — main GT vs predicted TPM scatter plots.
    extra_fig/ — all diagnostic / supporting figures (P1, P6).

    Args:
        exp_dir   : Path -- experiment folder (e.g. exprmnt_2026_03_28__00_51_02/).
        timestamp : str  -- run timestamp string.

    Returns:
        (out_dir, fig_dir, extra_fig_dir) as Path objects.
    """
    out_dir       = exp_dir / f"eda_{timestamp}"
    fig_dir       = out_dir / "figures"
    extra_fig_dir = out_dir / "extra_fig"
    os.makedirs(fig_dir,       exist_ok=True)
    os.makedirs(extra_fig_dir, exist_ok=True)
    return out_dir, fig_dir, extra_fig_dir


def main():
    """
    Run EDA priorities for both experiment folders independently.

    Each experiment folder gets its own eda_{timestamp}/ subfolder with the
    full set of figures and tables. This mirrors how run_gt_comparison.py and
    run_batch_comparison.py save outputs — each folder being analysed gets its
    own independent copy of the analysis results.

    Output structure (same in both experiment folders):
        exprmnt_.../
        └── eda_{timestamp}/
            ├── running.log
            ├── runtime.txt
            ├── priority1_universe_table.txt
            ├── priority6_leaked_stats.txt
            └── figures/
                └── fig_*.png

    Run:
        cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
        conda activate NanoCount_5
        python analysis/run_eda.py
    """
    run_start = time.time()
    timestamp = time.strftime("%Y_%m_%d__%H_%M_%S")

    lk_dir = Path(RESULTS_BASE) / LK_EXP
    jk_dir = Path(RESULTS_BASE) / JK_EXP
    samples = sorted(SAMPLE_GT_MAP.keys())

    # Validate experiment dirs exist
    for name, d in [("LK", lk_dir), ("JK MS", jk_dir)]:
        if not d.is_dir():
            print(f"ERROR: {name} experiment dir not found: {d}")
            sys.exit(1)

    # Create independent eda_{timestamp}/ inside EACH experiment folder
    lk_out, lk_fig, lk_extra = _setup_out_dir(lk_dir, timestamp)
    jk_out, jk_fig, jk_extra = _setup_out_dir(jk_dir, timestamp)
    print(f"LK output  : {lk_out}")
    print(f"JK output  : {jk_out}")
    print(f"Samples    : {samples}")

    # Open one log per experiment folder
    def _open_log(out_dir: Path) -> object:
        log = open(out_dir / "running.log", "w")
        log.write(f"Script: {__file__}\n")
        log.write(f"Timestamp: {timestamp}\n")
        log.write(f"LK exp : {LK_EXP}\n")
        log.write(f"JK exp : {JK_EXP}\n\n")
        return log

    lk_log = _open_log(lk_out)
    jk_log = _open_log(jk_out)

    # Combined log writer — writes to both experiment logs simultaneously
    class DualLog:
        def write(self, msg):
            lk_log.write(msg)
            jk_log.write(msg)

    dual_log = DualLog()

    import shutil

    # ── Main figures (figures/) — GT vs predicted TPM, both methods ───
    # Run once, save to lk_fig, copy to jk_fig (same data, independent copies)
    run_main_figures(samples, lk_dir, jk_dir, SAMPLE_GT_MAP, lk_fig, dual_log)
    for f in lk_fig.glob("fig_*.png"):
        shutil.copy2(f, jk_fig / f.name)

    if RUN_PRIORITY_1:
        # P1 figures → extra_fig/ in both folders
        run_priority_1(samples, lk_dir, jk_dir, SAMPLE_GT_MAP, lk_extra, lk_out, dual_log)
        for f in lk_extra.glob("fig_*.png"):
            shutil.copy2(f, jk_extra / f.name)
        shutil.copy2(lk_out / "priority1_universe_table.txt",
                     jk_out / "priority1_universe_table.txt")

    if RUN_PRIORITY_6:
        # P6 figures → extra_fig/ in both folders
        run_priority_6(samples, lk_dir, jk_dir, SAMPLE_GT_MAP, SAMPLE_BUSTOOLS_DIRS,
                       lk_extra, lk_out, dual_log)
        for f in lk_extra.glob("fig_p6*.png"):
            shutil.copy2(f, jk_extra / f.name)
        shutil.copy2(lk_out / "priority6_leaked_stats.txt",
                     jk_out / "priority6_leaked_stats.txt")

    # Save runtime to both folders
    elapsed = time.time() - run_start
    for out_dir in [lk_out, jk_out]:
        with open(out_dir / "runtime.txt", "w") as fh:
            fh.write(f"{elapsed:.2f} seconds\n")

    print(f"\nDone. Runtime: {elapsed:.1f}s")
    print(f"LK results : {lk_out}")
    print(f"JK results : {jk_out}")

    for log in [lk_log, jk_log]:
        log.write(f"\nRuntime: {elapsed:.2f}s\n")
        log.close()


if __name__ == "__main__":
    main()
