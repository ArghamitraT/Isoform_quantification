"""
plot_final_state.py
===================
Static version of the convergence animation: renders the 4-panel figure at the
FINAL training round only (no GIF).

One PNG per sample (sim1, sim2):

  ┌──────────────────────┬──────────────────────┐
  │  Panel 1: GT (static)│  Panel 2: JK theta   │
  ├──────────────────────┼──────────────────────┤
  │  Panel 3: JK MS alpha│  Panel 4: JK MS theta│
  └──────────────────────┴──────────────────────┘

Transcript universe per sample:
  - TP : GT TPM > 0  AND predicted non-zero (either method) → left portion
  - FP : GT TPM = 0  AND predicted non-zero (either method) → right of separator
  - FN are excluded.

Inputs:
  JK_SINGLE_DIR : experiment folder from run_joli_kallisto.sh (contains sim1/ sim2/)
  JK_MS_DIR     : experiment folder from run_multisample_joli.sh

Outputs (saved inside JK_MS_DIR):
  final_state_sim1.png
  final_state_sim2.png

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/plot_final_state.py
"""

import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# ============================================================
# CONFIG — edit before running; do not edit below
# ============================================================
RESULTS_BASE  = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# JK single-sample experiment folder (contains sim1/ and sim2/ subfolders)
JK_SINGLE_DIR = "exprmnt_2026_03_30__11_14_19"

# JK multi-sample experiment folder (contains snapshots.pkl + sim1/ sim2/)
JK_MS_DIR     = "exprmnt_2026_03_30__11_39_16"   # e.g. "exprmnt_2026_03_30__12_00_00"

SAMPLE_NAMES  = ["sim1", "sim2"]

GT_PATHS = {
    "sim1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv",
    "sim2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv",
}

COLORMAP     = "coolwarm"
LOG_SCALE_GT = True   # log y-axis for GT panel
FIG_SIZE     = (14, 8)
# ============================================================
# END CONFIG
# ============================================================


# ============================================================
# I/O helpers
# ============================================================

def load_gt(path: str) -> pd.DataFrame:
    """
    Load a ground truth file.

    Args:
        path : str -- path to GT TSV/CSV.

    Returns:
        pd.DataFrame with columns [transcript_id, tpm_gt].
    """
    df = pd.read_csv(path, index_col=0)
    id_col  = next((c for c in ["transcript_name", "transcript_id", "target_id"]
                    if c in df.columns), df.columns[0])
    val_col = next((c for c in ["tpm", "TPM"] if c in df.columns), df.columns[1])
    out = df[[id_col, val_col]].copy()
    out.columns = ["transcript_id", "tpm_gt"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["tpm_gt"]        = pd.to_numeric(out["tpm_gt"], errors="coerce").fillna(0.0)
    return out


def load_jk_single_final(exp_dir: str, sample_name: str) -> tuple:
    """
    Load the final theta snapshot from a JK single-sample run.

    Args:
        exp_dir     : str -- absolute path to experiment folder.
        sample_name : str -- sample subfolder name.

    Returns:
        (transcript_names: list[str], final_theta: np.ndarray)
    """
    path = os.path.join(exp_dir, sample_name, "theta_snapshots.pkl")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"JK single snapshot not found: {path}\n"
            "Re-run run_joli_kallisto.sh with SAVE_SNAPSHOTS=true."
        )
    with open(path, "rb") as fh:
        data = pickle.load(fh)
    _, final_theta = data["snapshots"][-1]
    print(f"  JK single {sample_name}: {len(data['snapshots'])} snapshots, "
          f"final round={data['snapshots'][-1][0]}")
    return data["transcript_names"], final_theta


def load_jk_ms_final(exp_dir: str) -> tuple:
    """
    Load the final alpha and per-sample theta from JK MS snapshots.

    Args:
        exp_dir : str -- absolute path to JK MS experiment folder.

    Returns:
        (sample_names: list[str], final_alpha: np.ndarray, final_thetas: list[np.ndarray])
    """
    path = os.path.join(exp_dir, "snapshots.pkl")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"JK MS snapshot not found: {path}\n"
            "Re-run run_multisample_joli.sh with SAVE_SNAPSHOTS=true."
        )
    with open(path, "rb") as fh:
        data = pickle.load(fh)
    final = data["snapshots"][-1]
    print(f"  JK MS: {len(data['snapshots'])} snapshots, final round={final['round']}")
    return data["sample_names"], final["alpha"], final["thetas"]


def load_transcript_names_from_abundance(exp_dir: str, sample_name: str) -> list:
    """
    Load transcript names from a sample's abundance.tsv inside an experiment folder.

    Args:
        exp_dir     : str -- absolute path to experiment folder.
        sample_name : str -- sample subfolder name.

    Returns:
        list[str] -- ordered transcript names.
    """
    path = os.path.join(exp_dir, sample_name, "abundance.tsv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"abundance.tsv not found: {path}")
    df = pd.read_csv(path, sep="\t")
    id_col = next((c for c in ["transcript_id", "target_id"] if c in df.columns),
                  df.columns[0])
    return list(df[id_col].astype(str))


# ============================================================
# Universe builder
# ============================================================

def build_universe(gt_df: pd.DataFrame, tx_names: list,
                   jk_theta: np.ndarray, jkms_theta: np.ndarray) -> tuple:
    """
    Build the sorted transcript universe (TP + FP) for one sample.

    TP = GT TPM > 0 AND predicted non-zero in at least one method.
    FP = GT TPM = 0 AND predicted non-zero in at least one method.
    FN are excluded.

    Args:
        gt_df      : pd.DataFrame -- [transcript_id, tpm_gt].
        tx_names   : list[str]    -- transcript names matching array indices.
        jk_theta   : np.ndarray   -- final JK single theta.
        jkms_theta : np.ndarray   -- final JK MS theta for this sample.

    Returns:
        (universe: pd.DataFrame, n_tp: int)
    """
    n = len(tx_names)
    df = pd.DataFrame({
        "transcript_id": tx_names,
        "jk_theta":      jk_theta[:n],
        "jkms_theta":    jkms_theta[:n],
    })
    df = df.merge(gt_df, on="transcript_id", how="left")
    df["tpm_gt"] = df["tpm_gt"].fillna(0.0)

    # GT-only universe: keep only transcripts where GT TPM > 0, sorted by GT TPM desc.
    df = df[df["tpm_gt"] > 0].copy()
    universe = df.sort_values("tpm_gt", ascending=False).reset_index(drop=True)
    universe["rank"] = np.arange(len(universe))
    n_tp = len(universe)
    print(f"  Universe: {n_tp} GT transcripts (FP excluded)")
    return universe, n_tp


# ============================================================
# Plotting
# ============================================================

def _panel_bounds(vals: np.ndarray) -> tuple:
    """
    Compute vmin/vmax using 1st–99th percentile of nonzero values.

    Args:
        vals : np.ndarray -- array of values.

    Returns:
        (vmin: float, vmax: float)
    """
    nz = vals[vals > 0]
    if len(nz) == 0:
        return 0.0, 1.0
    return float(np.percentile(nz, 1)), float(nz.max())


def plot_final_state(sample_name: str, universe: pd.DataFrame, n_tp: int,
                     jkms_alpha: np.ndarray, tx_index_ms: dict,
                     output_path: str) -> None:
    """
    Render and save the 4-panel static figure for one sample.

    Panels:
      [0,0] GT TPM   (static reference)
      [0,1] JK single theta (final round)
      [1,0] JK MS alpha    (final round, normalized)
      [1,1] JK MS theta    (final round)

    Args:
        sample_name  : str          -- e.g. "sim1".
        universe     : pd.DataFrame -- sorted universe with columns
                       [transcript_id, rank, tpm_gt, jk_theta, jkms_theta].
        n_tp         : int          -- number of TP transcripts.
        jkms_alpha   : np.ndarray   -- final JK MS alpha (T,).
        tx_index_ms  : dict         -- {transcript_id: array_index} for JK MS.
        output_path  : str          -- PNG save path.
    """
    ranks  = universe["rank"].values
    sep    = n_tp - 0.5

    gt_vals    = universe["tpm_gt"].values
    jk_vals    = universe["jk_theta"].values
    jkms_vals  = universe["jkms_theta"].values

    # Convert theta → pseudo-TPM (×1e6) so all panels share the same scale as GT.
    jk_tpm    = jk_vals   * 1e6
    jkms_tpm  = jkms_vals * 1e6

    # Alpha: show raw values — alpha is a Dirichlet concentration parameter,
    # not a probability. Normalizing by sum(alpha) gives misleadingly tiny values
    # because thousands of near-zero transcripts inflate the denominator.
    # Raw alpha: starts at ALPHA_INITIAL (e.g. 1.0), grows for expressed transcripts.
    alpha_tpm = np.array([jkms_alpha[tx_index_ms.get(tid, 0)]
                          for tid in universe["transcript_id"]])

    # Spearman vs GT (use TPM-scale values; rank order is identical to raw theta)
    gt_nz  = gt_vals > 0
    jk_sp  = spearmanr(gt_vals[gt_nz], jk_tpm[gt_nz]).statistic   \
              if gt_nz.sum() > 1 else float("nan")
    ms_sp  = spearmanr(gt_vals[gt_nz], jkms_tpm[gt_nz]).statistic  \
              if gt_nz.sum() > 1 else float("nan")

    # Theta panels share GT colorscale bounds so they are directly comparable.
    # Alpha uses its own bounds (raw Dirichlet concentration, independent scale).
    gt_vmin,  gt_vmax  = _panel_bounds(gt_vals)
    alp_vmin, alp_vmax = _panel_bounds(alpha_tpm)

    fig, axes = plt.subplots(2, 2, figsize=FIG_SIZE)
    fig.suptitle(
        f"{sample_name}  —  Final state\n"
        f"Theta panels share GT TPM scale  |  "
        f"Spearman vs GT:  JK={jk_sp:.4f}   JK MS={ms_sp:.4f}",
        fontsize=12, fontweight="bold"
    )

    panels = [
        (axes[0, 0], "GT  (reference TPM)",             gt_vals,   True,  gt_vmin,  gt_vmax),
        (axes[0, 1], "JK single-sample  θ × 1e6",       jk_tpm,    False, gt_vmin,  gt_vmax),
        (axes[1, 0], "JK MS  α  (raw Dirichlet conc.)", alpha_tpm, False, alp_vmin, alp_vmax),
        (axes[1, 1], "JK MS  θ × 1e6",                  jkms_tpm,  False, gt_vmin,  gt_vmax),
    ]

    for ax, title, vals, is_gt, vmin, vmax in panels:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        sc = ax.scatter(ranks, vals, c=vals, cmap=COLORMAP, norm=norm,
                        s=0.5, linewidths=0, rasterized=True)
        plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02)

        # TP / FP separator
        if n_tp < len(ranks):
            ax.axvline(x=sep, color="black", linewidth=1.0, linestyle="--", alpha=0.6)
            ax.text(sep + 1, vmax * 0.9, "FP →", fontsize=6, color="black", va="top")

        if is_gt and LOG_SCALE_GT:
            ax.set_yscale("log")

        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Transcript rank (GT order)", fontsize=7)
        ax.set_ylabel("Value", fontsize=7)
        ax.tick_params(labelsize=6)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path}")


# ============================================================
# Main
# ============================================================

def main() -> None:
    """
    Load final snapshots from JK single and JK MS, build universe per sample,
    and save one static 4-panel PNG per sample.
    """
    for name, val in [("JK_SINGLE_DIR", JK_SINGLE_DIR), ("JK_MS_DIR", JK_MS_DIR)]:
        if not val:
            print(f"ERROR: {name} is empty — set it in the CONFIG section.")
            sys.exit(1)

    jk_dir = os.path.join(RESULTS_BASE, JK_SINGLE_DIR)
    ms_dir = os.path.join(RESULTS_BASE, JK_MS_DIR)

    print("=" * 60)
    print("plot_final_state.py")
    print("=" * 60)
    print(f"  JK single dir : {JK_SINGLE_DIR}")
    print(f"  JK MS dir     : {JK_MS_DIR}")
    print()

    # Create figures/ subfolder inside JK MS experiment dir
    figures_dir = os.path.join(ms_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)
    print(f"  Output figures dir: {figures_dir}")

    # Load GT
    gt_dfs = {s: load_gt(GT_PATHS[s]) for s in SAMPLE_NAMES}

    # Load JK MS final state
    ms_sample_names, ms_alpha, ms_thetas = load_jk_ms_final(ms_dir)
    ms_sidx = {s: ms_sample_names.index(s) for s in SAMPLE_NAMES
               if s in ms_sample_names}

    # Transcript names + index from JK MS abundance.tsv
    tx_names_ms = load_transcript_names_from_abundance(ms_dir, ms_sample_names[0])
    tx_index_ms = {tid: i for i, tid in enumerate(tx_names_ms)}

    # Per-sample plots
    for sample in SAMPLE_NAMES:
        print(f"\n--- {sample} ---")

        # JK single final theta (aligned to JK single's own transcript order)
        jk_tx_names, jk_theta_raw = load_jk_single_final(jk_dir, sample)
        jk_tx_index = {tid: i for i, tid in enumerate(jk_tx_names)}

        # Re-index JK single theta to JK MS transcript order for a shared universe
        jk_theta_ms = np.array([
            jk_theta_raw[jk_tx_index[tid]] if tid in jk_tx_index else 0.0
            for tid in tx_names_ms
        ])

        s_idx = ms_sidx.get(sample, 0)
        jkms_theta = ms_thetas[s_idx]

        universe, n_tp = build_universe(
            gt_df      = gt_dfs[sample],
            tx_names   = tx_names_ms,
            jk_theta   = jk_theta_ms,
            jkms_theta = jkms_theta,
        )

        out_path = os.path.join(figures_dir, f"final_state_{sample}.png")
        plot_final_state(
            sample_name = sample,
            universe    = universe,
            n_tp        = n_tp,
            jkms_alpha  = ms_alpha,
            tx_index_ms = tx_index_ms,
            output_path = out_path,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
