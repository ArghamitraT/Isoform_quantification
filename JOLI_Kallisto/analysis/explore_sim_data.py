"""
explore_sim_data.py — EDA on two simulation samples from lr-kallisto output.

Answers four questions:
  1. How many reads and transcripts?
  2. How many ECs? How many single-tx vs multi-tx?
  3. Spearman/Pearson correlation between the two ground-truth TPM vectors.
  4. Theta (TPM) distribution for each sample.

Inputs (read from DATA_BASE):
  ds_100_num1_aln_21_long/  → Sample 2  (gt: PB_sample2_gt.tsv)
  ds_100_num1_aln_01_long/  → Sample 1  (gt: PB_sample1_gt.tsv)

Outputs:
  JOLI_Kallisto/analysis/figures/sim_eda_{timestamp}/
    reads_transcripts.txt       — Q1 summary table
    ec_breakdown.txt            — Q2 summary table
    ec_size_distribution.png    — Q2 EC-size histograms
    gt_correlation.txt          — Q3 Spearman/Pearson numbers
    gt_correlation_scatter.png  — Q3 scatter plot
    theta_distribution.png      — Q4 TPM histograms

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/explore_sim_data.py
"""

import json
import os
import sys
from collections import Counter
from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================
DATA_BASE   = "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data"
GT_DIR      = os.path.join(DATA_BASE, "ground_truths")
KAL_BASE    = os.path.join(DATA_BASE, "kallisto_output")

SAMPLE_21   = "ds_100_num1_aln_21_long"   # Sample 2
SAMPLE_01   = "ds_100_num1_aln_01_long"   # Sample 1
GT_21       = "PB_sample2_gt.tsv"
GT_01       = "PB_sample1_gt.tsv"

OUT_BASE    = os.path.join(
    os.path.dirname(__file__), "figures"
)
# ============================================================

def make_out_dir():
    """Create timestamped output directory under analysis/figures/."""
    ts = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
    out = os.path.join(OUT_BASE, f"sim_eda_{ts}")
    os.makedirs(out, exist_ok=True)
    print(f"[output] {out}")
    return out


# ---------------------------------------------------------------------------
# Q1: Reads & transcripts
# ---------------------------------------------------------------------------

def read_run_info(sample_dir):
    """
    Parse run_info.json for a kallisto output directory.

    Args:
        sample_dir (str): path to kallisto output folder.

    Returns:
        dict: selected fields from run_info.json.
    """
    path = os.path.join(sample_dir, "run_info.json")
    with open(path) as f:
        d = json.load(f)
    return {
        "n_targets":      d["n_targets"],
        "n_processed":    d["n_processed"],
        "n_pseudoaligned": d["n_pseudoaligned"],
        "p_pseudoaligned": d["p_pseudoaligned"],
        "n_unique":       d["n_unique"],
        "p_unique":       d["p_unique"],
    }


def count_covered_transcripts(ec_path):
    """
    Count how many unique transcript IDs appear across all ECs.

    Args:
        ec_path (str): path to matrix.ec file.

    Returns:
        int: number of unique transcripts covered.
    """
    tx_ids = set()
    with open(ec_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            for tid in parts[1].split(","):
                tx_ids.add(tid)
    return len(tx_ids)


def q1_reads_transcripts(out_dir, dirs):
    """
    Summarise read counts and transcript coverage for both samples.

    Args:
        out_dir (str): directory to write reads_transcripts.txt.
        dirs    (dict): {'21': path_to_ds21_dir, '01': path_to_ds01_dir}

    Returns:
        None
    """
    print("\n=== Q1: Reads & Transcripts ===")
    rows = []
    for label, d in dirs.items():
        info = read_run_info(d)
        ec_path = os.path.join(d, "matrix.ec")
        covered = count_covered_transcripts(ec_path)
        rows.append({
            "sample":            label,
            "reads_processed":   info["n_processed"],
            "reads_aligned":     info["n_pseudoaligned"],
            "pct_aligned":       info["p_pseudoaligned"],
            "reads_unique":      info["n_unique"],
            "pct_unique":        info["p_unique"],
            "transcripts_index": info["n_targets"],
            "transcripts_in_ec": covered,
        })

    df = pd.DataFrame(rows)
    print(df.to_string(index=False))

    out_path = os.path.join(out_dir, "reads_transcripts.txt")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  saved → {out_path}")


# ---------------------------------------------------------------------------
# Q2: EC breakdown
# ---------------------------------------------------------------------------

def parse_ec_sizes(ec_path):
    """
    Parse matrix.ec and return list of EC sizes (number of transcripts per EC).

    Args:
        ec_path (str): path to matrix.ec file.

    Returns:
        list[int]: one integer per EC, value = number of transcripts.
    """
    sizes = []
    with open(ec_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            n = len(parts[1].split(","))
            sizes.append(n)
    return sizes


def parse_ec_read_counts(mtx_path):
    """
    Parse count.mtx (MatrixMarket) and return per-EC read counts.

    Args:
        mtx_path (str): path to count.mtx file.

    Returns:
        dict[int, float]: ec_index (0-based) → read count.
    """
    counts = {}
    with open(mtx_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("%") or not line:
                continue
            parts = line.split()
            if len(parts) == 3 and parts[0] == "1" and parts[1] != parts[2]:
                # header line: "1 <n_ec> <n_ec>" — skip
                # real data lines have row=1 col=ec_index+1 value
                pass
            # MatrixMarket: row col value (1-indexed)
            if len(parts) == 3:
                try:
                    row, col, val = int(parts[0]), int(parts[1]), float(parts[2])
                    if row == 1:          # only one barcode (bulk)
                        counts[col - 1] = val   # convert to 0-based
                except ValueError:
                    continue
    return counts


def q2_ec_breakdown(out_dir, dirs):
    """
    Compute EC count, single/multi-tx split, size distribution,
    and total reads per EC. Plot and save.

    Args:
        out_dir (str): output directory.
        dirs    (dict): {'21': path, '01': path}

    Returns:
        dict: per-label size lists, for downstream use.
    """
    print("\n=== Q2: EC Breakdown ===")
    size_data = {}
    rows = []

    for label, d in dirs.items():
        ec_path  = os.path.join(d, "matrix.ec")
        mtx_path = os.path.join(d, "count.mtx")

        sizes  = parse_ec_sizes(ec_path)
        counts = parse_ec_read_counts(mtx_path)

        n_total  = len(sizes)
        n_single = sum(1 for s in sizes if s == 1)
        n_multi  = n_total - n_single
        total_reads = sum(counts.values())

        print(f"  ds_{label}: {n_total} ECs  |  single={n_single} ({100*n_single/n_total:.1f}%)  "
              f"multi={n_multi} ({100*n_multi/n_total:.1f}%)  |  total reads in ECs={total_reads:.0f}")

        rows.append({
            "sample":       label,
            "total_ecs":    n_total,
            "single_tx_ec": n_single,
            "multi_tx_ec":  n_multi,
            "pct_single":   round(100 * n_single / n_total, 1),
            "pct_multi":    round(100 * n_multi / n_total, 1),
            "total_reads_in_ecs": total_reads,
        })
        size_data[label] = sizes

    # Save text summary
    df = pd.DataFrame(rows)
    out_txt = os.path.join(out_dir, "ec_breakdown.txt")
    df.to_csv(out_txt, sep="\t", index=False)
    print(f"  saved → {out_txt}")

    # Plot EC-size distribution (cap at size 10 for readability)
    cap = 10
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=False)
    labels_plot = ["21 (Sample 2)", "01 (Sample 1)"]
    for ax, (label, sizes), title in zip(axes, size_data.items(), labels_plot):
        capped = [min(s, cap) for s in sizes]
        cnt = Counter(capped)
        xs = sorted(cnt.keys())
        ys = [cnt[x] for x in xs]
        xtick_labels = [str(x) if x < cap else f"≥{cap}" for x in xs]

        ax.bar(xs, ys, color="steelblue", edgecolor="white")
        ax.set_xticks(xs)
        ax.set_xticklabels(xtick_labels)
        ax.set_xlabel("EC size (# transcripts)")
        ax.set_ylabel("# ECs")
        ax.set_title(f"EC size distribution — ds_{title}")

        # annotate single vs multi counts
        n_s = cnt.get(1, 0)
        n_m = sum(v for k, v in cnt.items() if k > 1)
        ax.text(0.97, 0.97, f"single={n_s:,}\nmulti={n_m:,}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=9, bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    fig.suptitle("Equivalence Class Size Distribution", fontsize=13, fontweight="bold")
    fig.tight_layout()
    out_fig = os.path.join(out_dir, "ec_size_distribution.png")
    fig.savefig(out_fig, dpi=150)
    plt.close(fig)
    print(f"  saved → {out_fig}")

    return size_data


# ---------------------------------------------------------------------------
# Q3: GT correlation
# ---------------------------------------------------------------------------

def load_gt(gt_path):
    """
    Load a PB_sample*_gt.tsv file into a Series indexed by transcript_name.

    Args:
        gt_path (str): path to ground truth TSV (columns: index, transcript_name, tpm).

    Returns:
        pd.Series: TPM values indexed by transcript_name.
    """
    df = pd.read_csv(gt_path, index_col=0)
    return df.set_index("transcript_name")["tpm"]


def q3_gt_correlation(out_dir):
    """
    Compute and plot the Spearman/Pearson correlation between
    GT TPM vectors of the two samples.

    Args:
        out_dir (str): output directory.

    Returns:
        None
    """
    print("\n=== Q3: GT Correlation (ds_21 vs ds_01) ===")

    gt21 = load_gt(os.path.join(GT_DIR, GT_21))
    gt01 = load_gt(os.path.join(GT_DIR, GT_01))

    # Join on common transcripts
    merged = pd.DataFrame({"tpm_s2": gt21, "tpm_s1": gt01}).dropna()
    print(f"  Transcripts in GT21: {len(gt21):,}  |  GT01: {len(gt01):,}  |  shared: {len(merged):,}")

    # Log transform
    log_s2 = np.log10(merged["tpm_s2"] + 1)
    log_s1 = np.log10(merged["tpm_s1"] + 1)

    sp_raw, _ = spearmanr(merged["tpm_s2"], merged["tpm_s1"])
    pe_raw, _ = pearsonr(merged["tpm_s2"],  merged["tpm_s1"])
    sp_log, _ = spearmanr(log_s2, log_s1)
    pe_log, _ = pearsonr(log_s2,  log_s1)

    summary = (
        f"Transcripts in GT_21 (Sample 2): {len(gt21):,}\n"
        f"Transcripts in GT_01 (Sample 1): {len(gt01):,}\n"
        f"Shared transcripts:              {len(merged):,}\n\n"
        f"On raw TPM:\n"
        f"  Spearman rho = {sp_raw:.4f}\n"
        f"  Pearson  r   = {pe_raw:.4f}\n\n"
        f"On log10(TPM+1):\n"
        f"  Spearman rho = {sp_log:.4f}\n"
        f"  Pearson  r   = {pe_log:.4f}\n"
    )
    print(summary)
    out_txt = os.path.join(out_dir, "gt_correlation.txt")
    with open(out_txt, "w") as f:
        f.write(summary)
    print(f"  saved → {out_txt}")

    # Scatter plot (log scale)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.hexbin(log_s1, log_s2, gridsize=60, cmap="Blues", mincnt=1)
    ax.set_xlabel("log10(TPM+1)  — ds_01 (Sample 1)", fontsize=11)
    ax.set_ylabel("log10(TPM+1)  — ds_21 (Sample 2)", fontsize=11)
    ax.set_title("GT TPM Correlation", fontsize=13, fontweight="bold")

    # Diagonal reference line
    lo = min(log_s1.min(), log_s2.min())
    hi = max(log_s1.max(), log_s2.max())
    ax.plot([lo, hi], [lo, hi], "r--", lw=1, alpha=0.6, label="y=x")

    ax.text(0.05, 0.95,
            f"Spearman={sp_log:.4f}\nPearson={pe_log:.4f}\nn={len(merged):,}",
            transform=ax.transAxes, va="top", fontsize=10,
            bbox=dict(boxstyle="round", fc="white", alpha=0.8))
    ax.legend(fontsize=9)
    fig.tight_layout()
    out_fig = os.path.join(out_dir, "gt_correlation_scatter.png")
    fig.savefig(out_fig, dpi=150)
    plt.close(fig)
    print(f"  saved → {out_fig}")


# ---------------------------------------------------------------------------
# Q4: Theta (TPM) distribution
# ---------------------------------------------------------------------------

def q4_theta_distribution(out_dir):
    """
    Plot the GT TPM distribution for each sample as log10(TPM+1) histograms,
    restricted to transcripts with TPM > 0.

    Args:
        out_dir (str): output directory.

    Returns:
        None
    """
    print("\n=== Q4: Theta (TPM) Distribution ===")

    gt21 = load_gt(os.path.join(GT_DIR, GT_21))
    gt01 = load_gt(os.path.join(GT_DIR, GT_01))

    for label, gt in [("21 (Sample 2)", gt21), ("01 (Sample 1)", gt01)]:
        nonzero = gt[gt > 0]
        pct_nz  = 100 * len(nonzero) / len(gt)
        print(f"  ds_{label}: total={len(gt):,}  nonzero={len(nonzero):,} ({pct_nz:.1f}%)  "
              f"median_tpm={nonzero.median():.2f}  max_tpm={nonzero.max():.1f}")

    gt21_nz = gt21[gt21 > 0]
    gt01_nz = gt01[gt01 > 0]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    bin_edges = np.linspace(0, max(np.log10(gt21_nz + 1).max(),
                                   np.log10(gt01_nz + 1).max()) + 0.1, 60)

    for ax, (label, tpm), color in zip(
        axes,
        [("21 (Sample 2)", gt21_nz), ("01 (Sample 1)", gt01_nz)],
        ["steelblue", "darkorange"]
    ):
        log_vals = np.log10(tpm + 1)
        ax.hist(log_vals, bins=bin_edges, color=color, edgecolor="white", alpha=0.85)
        ax.set_xlabel("log10(TPM+1)", fontsize=11)
        ax.set_ylabel("# transcripts", fontsize=11)
        ax.set_title(f"TPM distribution — ds_{label}\n"
                     f"(nonzero: {len(tpm):,}, median={tpm.median():.1f})", fontsize=11)

    fig.suptitle("Ground-Truth TPM (Theta) Distribution", fontsize=13, fontweight="bold")
    fig.tight_layout()
    out_fig = os.path.join(out_dir, "theta_distribution.png")
    fig.savefig(out_fig, dpi=150)
    plt.close(fig)
    print(f"  saved → {out_fig}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """
    Run all four EDA sections and save results.

    Run:
        cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
        conda activate NanoCount_5
        python analysis/explore_sim_data.py
    """
    dirs = {
        "21": os.path.join(KAL_BASE, SAMPLE_21),
        "01": os.path.join(KAL_BASE, SAMPLE_01),
    }

    # Sanity-check paths
    for label, d in dirs.items():
        assert os.path.isdir(d), f"Missing directory: {d}"
    for gt_name in [GT_21, GT_01]:
        path = os.path.join(GT_DIR, gt_name)
        assert os.path.isfile(path), f"Missing GT file: {path}"

    out_dir = make_out_dir()

    q1_reads_transcripts(out_dir, dirs)
    q2_ec_breakdown(out_dir, dirs)
    q3_gt_correlation(out_dir)
    q4_theta_distribution(out_dir)

    print(f"\n[done] all outputs in {out_dir}")


if __name__ == "__main__":
    main()
