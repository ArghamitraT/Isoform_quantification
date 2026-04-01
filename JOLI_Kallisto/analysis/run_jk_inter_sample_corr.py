"""
run_jk_inter_sample_corr.py
============================
Compute and plot the inter-sample Spearman correlation between sim1 and sim2
theta estimates across EM rounds, using theta_snapshots.pkl from a JK single-
sample experiment folder.

Answers: how does the similarity between the two sample theta estimates evolve
across EM iterations?

Inputs:
  JK_EXP_DIR   : experiment folder containing sim1/ and sim2/ subfolders,
                 each with theta_snapshots.pkl
  SAMPLE_NAMES : names of the two sample subfolders to compare

Outputs (saved inside JK_EXP_DIR):
  jk_inter_sample_spearman.png  -- Spearman(theta_sim1, theta_sim2) vs EM round

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/run_jk_inter_sample_corr.py
"""

import os
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr

# ============================================================
# CONFIG — edit before running; do not edit below
# ============================================================
RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"
JK_EXP_DIR   = "exprmnt_2026_03_30__11_14_19"
SAMPLE_NAMES = ["sim1", "sim2"]   # must match subfolder names inside JK_EXP_DIR
# ============================================================
# END CONFIG
# ============================================================


def load_snapshots(exp_dir: str, sample_name: str) -> dict:
    """
    Load theta snapshots from a JK single-sample experiment subfolder.

    Args:
        exp_dir     : str -- absolute path to experiment folder.
        sample_name : str -- sample subfolder name (e.g. "sim1").

    Returns:
        dict with keys:
          "transcript_names" : list[str]
          "snapshots"        : list of (round_num, theta_array)
    """
    path = os.path.join(exp_dir, sample_name, "theta_snapshots.pkl")
    if not os.path.exists(path):
        raise FileNotFoundError(f"Snapshot not found: {path}")
    with open(path, "rb") as fh:
        data = pickle.load(fh)
    print(f"  Loaded {len(data['snapshots'])} snapshots for {sample_name}: {path}")
    return data


def align_transcripts(data1: dict, data2: dict) -> tuple:
    """
    Find shared transcripts between two snapshot dicts and return aligned index arrays.

    Args:
        data1 : dict -- snapshots dict for sample 1.
        data2 : dict -- snapshots dict for sample 2.

    Returns:
        (idx1, idx2) : (np.ndarray, np.ndarray) -- indices into each sample's theta
                       array for the shared transcript set.
    """
    names1 = {n: i for i, n in enumerate(data1["transcript_names"])}
    names2 = {n: i for i, n in enumerate(data2["transcript_names"])}
    shared = sorted(set(names1) & set(names2))
    print(f"  Shared transcripts: {len(shared)} "
          f"(of {len(names1)} / {len(names2)})")
    idx1 = np.array([names1[n] for n in shared])
    idx2 = np.array([names2[n] for n in shared])
    return idx1, idx2


def compute_correlations(data1: dict, data2: dict,
                         idx1: np.ndarray, idx2: np.ndarray) -> tuple:
    """
    Compute Spearman and Pearson correlation between the two sample thetas at each
    snapshot round.

    Args:
        data1 : dict      -- snapshots dict for sample 1.
        data2 : dict      -- snapshots dict for sample 2.
        idx1  : np.ndarray -- shared transcript indices for sample 1.
        idx2  : np.ndarray -- shared transcript indices for sample 2.

    Returns:
        (rounds, spearmans, pearsons) : three lists of equal length.
    """
    snaps1 = {r: theta for r, theta in data1["snapshots"]}
    snaps2 = {r: theta for r, theta in data2["snapshots"]}
    shared_rounds = sorted(set(snaps1) & set(snaps2))

    rounds, spearmans, pearsons = [], [], []
    for r in shared_rounds:
        t1 = snaps1[r][idx1]
        t2 = snaps2[r][idx2]
        sp = spearmanr(t1, t2)
        pe = pearsonr(t1, t2)
        spearmans.append(sp.statistic if hasattr(sp, "statistic") else sp[0])
        pearsons.append(pe.statistic if hasattr(pe, "statistic") else pe[0])
        rounds.append(r)

    return rounds, spearmans, pearsons


def plot_and_save(rounds: list, spearmans: list, pearsons: list,
                  sample_names: list, output_path: str) -> None:
    """
    Plot Spearman and Pearson inter-sample correlation vs EM round and save to PNG.

    Args:
        rounds       : list[int]   -- EM round numbers.
        spearmans    : list[float] -- Spearman values per round.
        pearsons     : list[float] -- Pearson values per round.
        sample_names : list[str]   -- the two sample names being compared.
        output_path  : str         -- path to save the PNG.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(f"JK inter-sample correlation: {sample_names[0]} vs {sample_names[1]}")

    for ax, vals, label in [
        (axes[0], spearmans, "Spearman"),
        (axes[1], pearsons,  "Pearson"),
    ]:
        ax.plot(rounds, vals, marker="o", markersize=3, linewidth=1.5)
        ax.set_xlabel("EM round")
        ax.set_ylabel(f"{label} correlation")
        ax.set_title(f"{label}: theta {sample_names[0]} vs theta {sample_names[1]}")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-0.05, 1.05)
        # Annotate final value
        ax.annotate(f"final: {vals[-1]:.4f}",
                    xy=(rounds[-1], vals[-1]),
                    xytext=(-60, 10), textcoords="offset points",
                    fontsize=9, arrowprops=dict(arrowstyle="->", lw=0.8))

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


def main() -> None:
    """
    Load snapshots, align transcripts, compute inter-sample correlation per round,
    plot and save.
    """
    exp_dir = os.path.join(RESULTS_BASE, JK_EXP_DIR)

    print("=" * 55)
    print("JK inter-sample theta correlation across EM rounds")
    print("=" * 55)
    print(f"  Experiment : {JK_EXP_DIR}")
    print(f"  Samples    : {SAMPLE_NAMES[0]} vs {SAMPLE_NAMES[1]}")
    print()

    data1 = load_snapshots(exp_dir, SAMPLE_NAMES[0])
    data2 = load_snapshots(exp_dir, SAMPLE_NAMES[1])

    idx1, idx2 = align_transcripts(data1, data2)
    rounds, spearmans, pearsons = compute_correlations(data1, data2, idx1, idx2)

    print(f"\n  Rounds computed : {len(rounds)}  (first={rounds[0]}, last={rounds[-1]})")
    print(f"  Final Spearman  : {spearmans[-1]:.6f}")
    print(f"  Final Pearson   : {pearsons[-1]:.6f}")

    figures_dir = os.path.join(exp_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)
    output_path = os.path.join(figures_dir, "jk_inter_sample_spearman.png")
    plot_and_save(rounds, spearmans, pearsons, SAMPLE_NAMES, output_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
