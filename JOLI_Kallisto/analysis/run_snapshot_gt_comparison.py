"""
run_snapshot_gt_comparison.py
==============================
Compare per-iteration theta (and optionally alpha) snapshots from two experiments
against simulation ground truth.

Evaluates JK multi-sample (em_wrapper) and JK single-sample EM at specific
training iterations (e.g. 20, 25, 30) to track how quantification quality
improves across training. Useful for diagnosing convergence speed vs accuracy.

Snapshot formats expected
--------------------------
JK multi-sample  ({JK_MS_EXP}/snapshots.pkl):
    dict {
        "sample_names": list[str],
        "snapshots":    list of dicts {
            "round":  int,
            "alpha":  np.ndarray (T,)   -- shared Dirichlet hyperparameter
            "thetas": list[np.ndarray]  -- per-sample normalized theta (T,)
        }
    }

JK single-sample  ({JK_SS_EXP}/{sample}/theta_snapshots.pkl):
    dict {
        "transcript_names": list[str],
        "snapshots":        list of (round_num: int, theta: np.ndarray (T,))
    }

TPM conversion (no total_reads needed — cancels out in normalization)
----------------------------------------------------------------------
    rho[t]  = theta[t] / eff_len[t]
    tpm[t]  = rho[t] / sum(rho) * 1e6

eff_lens are read from the existing {exp}/{sample}/abundance.tsv (eff_length column),
which is already computed and saved by write_abundance().

Metrics (same as run_gt_comparison.py)
---------------------------------------
Spearman, Pearson, MAE, RMSE, Bray-Curtis across 3 universes:
  1. all        — outer join (includes TN pairs)
  2. active     — GT ∪ non-zero predicted (fair cross-method comparison)
  3. gt_nonzero — only transcripts present in GT

Plus contingency table: TP / FP / FN / TN.

Outputs (Rule 2 — multi-file → named subfolder, written to BOTH experiments)
-----------------------------------------------------------------------------
    {exp}/snapshot_gt_comparison_{timestamp}/
        metrics_summary.txt          -- compact table: iter × sample × metric
        detail_iter_{round}.txt      -- full formatted block per iteration
        figures/
            metrics_convergence_{timestamp}.png  -- samples × metrics grid
            alpha_metrics_{timestamp}.png        -- (JK MS only) alpha vs GT

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/run_snapshot_gt_comparison.py
"""

from __future__ import annotations

import pickle
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================

RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# Experiment folder names
JK_MS_EXP = "exprmnt_2026_03_30__22_37_55"   # JK multi-sample (em_wrapper)
JK_SS_EXP = "exprmnt_2026_03_30__22_41_19"   # JK single-sample

# Iterations to evaluate (must exist in snapshots)
TARGET_ITERS = [-1, 20, 25, 30]

# Also compare the shared Dirichlet alpha (JK MS only) against GT?
COMPARE_ALPHA = True

# Map sample subfolder name → absolute path to ground truth TSV.
SAMPLE_GT_MAP = {
    "sim1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv",
    "sim2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv",
}

# ============================================================
# END CONFIG
# ============================================================

TIMESTAMP = time.strftime("%Y_%m_%d__%H_%M_%S")
METRICS   = ["spearman", "pearson", "mae", "rmse", "bray_curtis"]
UNIVERSES = ["all", "active", "gt"]   # short keys; displayed as full names below
UNIVERSE_LABELS = {
    "all":    "All transcripts",
    "active": "Active universe (GT ∪ nonzero pred)",
    "gt":     "GT non-zero only",
}


# ============================================================
# Metric helpers (identical to run_gt_comparison.py)
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
    Compute Spearman, Pearson, MAE, RMSE, Bray-Curtis between two TPM arrays.

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
# I/O helpers
# ============================================================

def read_ground_truth(path: Path) -> pd.DataFrame:
    """
    Read a ground truth CSV (comma-separated, first column is unnamed index).

    Args:
        path : Path -- path to ground truth TSV/CSV.

    Returns:
        pd.DataFrame with columns [transcript_id, tpm_gt].
    """
    df = pd.read_csv(path, index_col=0)
    id_candidates  = ["transcript_name", "transcript_id", "target_id", "Name"]
    val_candidates = ["tpm", "TPM", "est_counts"]
    id_col  = next((c for c in id_candidates  if c in df.columns), df.columns[0])
    val_col = next((c for c in val_candidates if c in df.columns), df.columns[1])
    out = df[[id_col, val_col]].copy()
    out.columns = ["transcript_id", "tpm_gt"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["tpm_gt"] = pd.to_numeric(out["tpm_gt"], errors="coerce").fillna(0.0)
    return out


def read_eff_lens(abundance_path: Path) -> tuple[list[str], np.ndarray]:
    """
    Read transcript names and effective lengths from an existing abundance.tsv.

    Uses the already-computed eff_length column so we don't need to reload
    bustools data.

    Args:
        abundance_path : Path -- path to {exp}/{sample}/abundance.tsv.

    Returns:
        tuple(list[str], np.ndarray): transcript_names, eff_lens (float64, T).
    """
    df = pd.read_csv(abundance_path, sep="\t", comment="#")
    id_col  = "target_id"   if "target_id"  in df.columns else df.columns[0]
    len_col = "eff_length"  if "eff_length" in df.columns else df.columns[2]
    transcript_names = df[id_col].astype(str).tolist()
    eff_lens         = pd.to_numeric(df[len_col], errors="coerce").fillna(1.0).to_numpy(dtype=np.float64)
    return transcript_names, eff_lens


def theta_to_tpm(theta: np.ndarray, eff_lens: np.ndarray) -> np.ndarray:
    """
    Convert normalized theta to TPM using effective lengths.

    total_reads cancels out in the normalization, so only theta and eff_lens
    are needed:
        rho[t]  = theta[t] / eff_len[t]
        tpm[t]  = rho[t] / sum(rho) * 1e6

    Args:
        theta    : np.ndarray (float64, T) -- normalized abundances (sums to ~1).
        eff_lens : np.ndarray (float64, T) -- effective transcript lengths.

    Returns:
        np.ndarray (float64, T) -- TPM values (sums to ~1e6).
    """
    safe_eff = np.where(eff_lens > 0, eff_lens, 1.0)
    rho      = theta / safe_eff
    rho_sum  = rho.sum()
    return (rho / rho_sum * 1e6) if rho_sum > 0 else np.zeros_like(rho)


def alpha_to_tpm(alpha: np.ndarray, eff_lens: np.ndarray) -> np.ndarray:
    """
    Normalize Dirichlet alpha to TPM using effective lengths.

    alpha is first normalized to a probability vector (same scale as theta),
    then converted to TPM via eff_lens.

    Args:
        alpha    : np.ndarray (float64, T) -- raw Dirichlet concentration values.
        eff_lens : np.ndarray (float64, T) -- effective transcript lengths.

    Returns:
        np.ndarray (float64, T) -- TPM values (sums to ~1e6).
    """
    alpha_sum = alpha.sum()
    theta_alpha = alpha / alpha_sum if alpha_sum > 0 else alpha.copy()
    return theta_to_tpm(theta_alpha, eff_lens)


# ============================================================
# GT comparison (one theta/alpha array vs one ground truth)
# ============================================================

def compare_to_gt(
    tpm_pred:         np.ndarray,
    transcript_names: list[str],
    gt_df:            pd.DataFrame,
    label:            str,
) -> dict:
    """
    Compare a TPM prediction array against ground truth across 3 universes.

    Universes:
      1. all        — outer join, includes TN (0,0) pairs
      2. active     — GT ∪ non-zero predicted, excludes TN pairs
      3. gt_nonzero — only transcripts with GT > 0

    Args:
        tpm_pred         : np.ndarray    -- predicted TPM (T,), ordered by transcript index.
        transcript_names : list[str]     -- transcript name at each index.
        gt_df            : pd.DataFrame  -- columns [transcript_id, tpm_gt].
        label            : str           -- label for printing (e.g. "sim1 iter=20").

    Returns:
        dict with contingency counts, universe sizes, and prefixed metric keys.
    """
    pred_df = pd.DataFrame({
        "transcript_id": transcript_names,
        "tpm_pred":      tpm_pred,
    })

    # Universe 1: outer join (all transcripts, including TN)
    merged_all = pred_df.merge(gt_df, on="transcript_id", how="outer")
    merged_all["tpm_pred"] = merged_all["tpm_pred"].fillna(0.0)
    merged_all["tpm_gt"]   = merged_all["tpm_gt"].fillna(0.0)

    pred_nz = merged_all["tpm_pred"] > 0
    gt_nz   = merged_all["tpm_gt"]   > 0
    n_tp = int(( pred_nz &  gt_nz).sum())
    n_fp = int(( pred_nz & ~gt_nz).sum())
    n_fn = int((~pred_nz &  gt_nz).sum())
    n_tn = int((~pred_nz & ~gt_nz).sum())
    n_total_all = len(merged_all)

    x_all = merged_all["tpm_pred"].to_numpy()
    y_all = merged_all["tpm_gt"].to_numpy()

    # Universe 2: active (GT ∪ nonzero pred), excludes TN
    pred_nz_df    = pred_df[pred_df["tpm_pred"] > 0]
    merged_active = pred_nz_df.merge(gt_df, on="transcript_id", how="outer")
    merged_active["tpm_pred"] = merged_active["tpm_pred"].fillna(0.0)
    merged_active["tpm_gt"]   = merged_active["tpm_gt"].fillna(0.0)
    x_active = merged_active["tpm_pred"].to_numpy()
    y_active = merged_active["tpm_gt"].to_numpy()
    n_active = len(merged_active)

    # Universe 3: GT non-zero only
    gt_nonzero_mask = merged_all["tpm_gt"] > 0
    x_gt = merged_all.loc[gt_nonzero_mask, "tpm_pred"].to_numpy()
    y_gt = merged_all.loc[gt_nonzero_mask, "tpm_gt"].to_numpy()
    n_gt_nonzero = int(gt_nonzero_mask.sum())

    results = {
        "label":        label,
        "n_total_all":  n_total_all,
        "n_tp": n_tp, "n_fp": n_fp, "n_fn": n_fn, "n_tn": n_tn,
        "n_active":     n_active,
        "n_gt_nonzero": n_gt_nonzero,
        **{f"all_{k}":    v for k, v in compute_metrics(x_all,    y_all).items()},
        **{f"active_{k}": v for k, v in compute_metrics(x_active, y_active).items()},
        **{f"gt_{k}":     v for k, v in compute_metrics(x_gt,     y_gt).items()},
    }
    return results


# ============================================================
# Formatting
# ============================================================

def _metric_block(prefix: str, results: dict) -> list[str]:
    """
    Build formatted metric lines for one universe prefix (all/active/gt).

    Args:
        prefix  : str  -- universe prefix key.
        results : dict -- output of compare_to_gt().

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


def format_result_block(r: dict) -> str:
    """
    Format one result dict (from compare_to_gt) into a human-readable text block.

    Args:
        r : dict -- output of compare_to_gt().

    Returns:
        str -- formatted multi-line block.
    """
    lines = [
        f"  Label : {r['label']}",
        "",
        "  Contingency table",
        "  " + "-" * 50,
        f"  True  positives  (pred>0, GT>0) : {r['n_tp']:>10,}",
        f"  False positives  (pred>0, GT=0) : {r['n_fp']:>10,}",
        f"  False negatives  (pred=0, GT>0) : {r['n_fn']:>10,}",
        f"  True  negatives  (pred=0, GT=0) : {r['n_tn']:>10,}",
        f"  Total (outer join universe)     : {r['n_total_all']:>10,}",
        "",
        f"  --- Metrics (all transcripts, N={r['n_total_all']:,}) ---",
        *_metric_block("all", r),
        "",
        f"  --- Metrics (active universe: GT ∪ non-zero pred, N={r['n_active']:,}) ---",
        "  (excludes true-negative (0,0) pairs — fair cross-method comparison)",
        *_metric_block("active", r),
        "",
        f"  --- Metrics (GT non-zero transcripts only, N={r['n_gt_nonzero']:,}) ---",
        "  (only transcripts present in the simulation ground truth)",
        *_metric_block("gt", r),
    ]
    return "\n".join(lines)


# ============================================================
# Snapshot loading
# ============================================================

def load_ms_snapshots(exp_dir: Path) -> dict:
    """
    Load JK multi-sample snapshots.pkl from an experiment folder.

    Expected format: dict { sample_names, snapshots: [{round, alpha, thetas}] }

    Args:
        exp_dir : Path -- experiment root directory.

    Returns:
        dict { "sample_names": list[str], "snapshots": list[dict] }

    Raises:
        FileNotFoundError if snapshots.pkl is not present.
    """
    snap_path = exp_dir / "snapshots.pkl"
    if not snap_path.exists():
        raise FileNotFoundError(f"snapshots.pkl not found: {snap_path}")
    with open(snap_path, "rb") as fh:
        data = pickle.load(fh)
    rounds = [s["round"] for s in data["snapshots"]]
    print(f"[load_ms_snapshots] Loaded {len(rounds)} snapshots, rounds: {rounds}")
    return data


def load_ss_snapshots(exp_dir: Path, sample: str) -> dict:
    """
    Load JK single-sample theta_snapshots.pkl for one sample.

    Expected format: dict { transcript_names, snapshots: [(round_num, theta)] }

    Args:
        exp_dir : Path -- experiment root directory.
        sample  : str  -- sample subfolder name (e.g. "sim1").

    Returns:
        dict { "transcript_names": list[str], "snapshots": list[tuple] }

    Raises:
        FileNotFoundError if theta_snapshots.pkl is not present.
    """
    snap_path = exp_dir / sample / "theta_snapshots.pkl"
    if not snap_path.exists():
        raise FileNotFoundError(f"theta_snapshots.pkl not found: {snap_path}")
    with open(snap_path, "rb") as fh:
        data = pickle.load(fh)
    rounds = [s[0] for s in data["snapshots"]]
    print(f"[load_ss_snapshots] {sample}: {len(rounds)} snapshots, rounds: {rounds}")
    return data


def get_ms_snapshot(ms_data: dict, target_iter: int) -> dict | None:
    """
    Retrieve the snapshot dict for a given iteration from JK MS snapshots.

    Args:
        ms_data     : dict -- output of load_ms_snapshots().
        target_iter : int  -- target round number.

    Returns:
        Matching snapshot dict or None if not found.
    """
    for snap in ms_data["snapshots"]:
        if snap["round"] == target_iter:
            return snap
    return None


def get_ss_snapshot(ss_data: dict, target_iter: int) -> np.ndarray | None:
    """
    Retrieve the theta array for a given iteration from JK SS snapshots.

    Args:
        ss_data     : dict -- output of load_ss_snapshots().
        target_iter : int  -- target round number.

    Returns:
        theta np.ndarray or None if not found.
    """
    for (round_num, theta) in ss_data["snapshots"]:
        if round_num == target_iter:
            return theta
    return None


# ============================================================
# Summary table builder
# ============================================================

def build_summary_table(all_records: list[dict]) -> str:
    """
    Build a compact plain-text table from a list of result dicts.

    Columns: label | universe | spearman | pearson | mae | rmse | bray_curtis

    Args:
        all_records : list[dict] -- list of compare_to_gt() outputs.

    Returns:
        str -- tab-formatted table.
    """
    header = f"{'label':<40} {'universe':<10} {'spearman':>10} {'pearson':>10} {'mae':>12} {'rmse':>12} {'bray_curtis':>12}"
    sep    = "-" * len(header)
    rows   = [header, sep]
    for r in all_records:
        for u in UNIVERSES:
            row = (
                f"{r['label']:<40} "
                f"{u:<10} "
                f"{fmt(r[f'{u}_spearman']):>10} "
                f"{fmt(r[f'{u}_pearson']):>10} "
                f"{fmt(r[f'{u}_mae']):>12} "
                f"{fmt(r[f'{u}_rmse']):>12} "
                f"{fmt(r[f'{u}_bray_curtis']):>12}"
            )
            rows.append(row)
        rows.append("")   # blank line between records
    return "\n".join(rows)


# ============================================================
# Figure generation
# ============================================================

def _plot_sample_bar_figure(
    sample_name:     str,
    iter_dict:       dict[int, dict],
    target_iters:    list[int],
    title:           str,
    out_path:        Path,
) -> None:
    """
    Bar plot of all metrics vs iteration for ONE sample.

    Layout: 1 row × 5 subplots (one per metric).
    Each subplot: grouped bars where
      - x-groups = iterations (e.g. -1, 20, 25, 30)
      - bars within group = universes (all / active / gt-nonzero), 3 bars per group.

    Per CLAUDE.md figure rules:
      - No overlapping plots — separate subplots per metric.
      - Group related figures onto one page using subplots.

    Args:
        sample_name  : str              -- sample label (used in y-axis / title).
        iter_dict    : dict[int, dict]  -- {iter: compare_to_gt() result dict}.
        target_iters : list[int]        -- ordered iteration numbers (x-axis groups).
        title        : str              -- figure suptitle.
        out_path     : Path             -- where to save the PNG.
    """
    n_metrics  = len(METRICS)
    n_univs    = len(UNIVERSES)
    bar_width  = 0.22
    x          = np.arange(len(target_iters))

    # One distinct color per universe — no overlapping bars, no transparency layering
    univ_colors = {
        "all":    "#4C72B0",   # muted blue
        "active": "#DD8452",   # muted orange
        "gt":     "#55A868",   # muted green
    }
    univ_display = {
        "all":    "All transcripts",
        "active": "Active (GT ∪ nonzero pred)",
        "gt":     "GT non-zero only",
    }

    fig, axes = plt.subplots(1, n_metrics,
                             figsize=(4.5 * n_metrics, 4.5),
                             squeeze=False)

    for m_idx, metric in enumerate(METRICS):
        ax = axes[0][m_idx]

        for u_idx, univ in enumerate(UNIVERSES):
            key    = f"{univ}_{metric}"
            offset = (u_idx - n_univs / 2 + 0.5) * bar_width
            vals   = []
            for it in target_iters:
                r = iter_dict.get(it)
                v = r.get(key, float("nan")) if r is not None else float("nan")
                vals.append(v)

            bars = ax.bar(
                x + offset, vals, bar_width,
                label  = univ_display[univ],
                color  = univ_colors[univ],
                edgecolor = "white",
                linewidth = 0.5,
            )

            # Annotate bar tops with values (skip NaN)
            for bar, v in zip(bars, vals):
                if not np.isnan(v):
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height() + 0.002 * (ax.get_ylim()[1] - ax.get_ylim()[0] + 1e-9),
                        f"{v:.3f}",
                        ha="center", va="bottom",
                        fontsize=5.5, rotation=90,
                    )

        ax.set_title(metric.replace("_", " ").title(), fontsize=10, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels([str(it) for it in target_iters], fontsize=9)
        ax.set_xlabel("Iteration", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(True, axis="y", linestyle="--", alpha=0.35)
        ax.set_axisbelow(True)

        # Legend only on first subplot to avoid repetition
        if m_idx == 0:
            ax.legend(fontsize=7, loc="best")

    fig.suptitle(title, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] Saved: {out_path}")


# ============================================================
# Per-experiment analysis
# ============================================================

def run_ms_experiment(
    exp_dir:    Path,
    ms_data:    dict,
    target_iters: list[int],
    sample_gt_map: dict,
) -> tuple[list[dict], dict[str, dict[str, list[dict]]]]:
    """
    Run GT comparison for all target iterations for the JK multi-sample experiment.

    For each (iteration, sample): load theta from snapshot, convert to TPM,
    compare against GT. Optionally also compare alpha against GT.

    Args:
        exp_dir       : Path -- experiment root directory.
        ms_data       : dict -- output of load_ms_snapshots().
        target_iters  : list[int]  -- iterations to evaluate.
        sample_gt_map : dict       -- {sample_name: gt_path}.

    Returns:
        tuple:
          all_records          : list[dict] -- flat list of compare_to_gt() results.
          records_by_series    : dict[series_label, dict[iter → result_dict]]
                                 (for plotting: series_label = sample name or "alpha")
    """
    sample_names = ms_data["sample_names"]
    all_records  = []

    # {series_label: {iter: result_dict}}  — used for figure plotting
    series_dict: dict[str, dict[int, dict]] = {}
    for sname in sample_names:
        series_dict[f"JK-MS {sname}"] = {}
    if COMPARE_ALPHA:
        series_dict["JK-MS alpha"] = {}

    for it in target_iters:
        snap = get_ms_snapshot(ms_data, it)
        if snap is None:
            print(f"[JK MS] iter={it} not found in snapshots, skipping.")
            continue

        for s_idx, sname in enumerate(sample_names):
            if sname not in sample_gt_map:
                continue
            gt_df = read_ground_truth(Path(sample_gt_map[sname]))

            # Load eff_lens from existing abundance.tsv
            abund_path = exp_dir / sname / "abundance.tsv"
            if not abund_path.exists():
                print(f"[JK MS] abundance.tsv not found for {sname}, skipping.")
                continue
            transcript_names, eff_lens = read_eff_lens(abund_path)

            # Convert theta to TPM
            theta   = snap["thetas"][s_idx]
            tpm_pred = theta_to_tpm(theta, eff_lens)

            label  = f"JK-MS {sname} iter={it}"
            result = compare_to_gt(tpm_pred, transcript_names, gt_df, label)
            all_records.append(result)
            series_dict[f"JK-MS {sname}"][it] = result

            print(f"  [JK MS] {sname} iter={it}: "
                  f"Spearman(active)={fmt(result['active_spearman'])}, "
                  f"TP={result['n_tp']}, FP={result['n_fp']}, FN={result['n_fn']}")

        # Alpha comparison (uses sim1 eff_lens as reference — same transcriptome)
        if COMPARE_ALPHA:
            ref_sname = next((s for s in sample_names if s in sample_gt_map), None)
            if ref_sname is not None:
                abund_path = exp_dir / ref_sname / "abundance.tsv"
                if abund_path.exists():
                    transcript_names, eff_lens = read_eff_lens(abund_path)
                    tpm_alpha = alpha_to_tpm(snap["alpha"], eff_lens)
                    # Compare alpha against each sample's GT and store separately
                    for sname in sample_names:
                        if sname not in sample_gt_map:
                            continue
                        gt_df  = read_ground_truth(Path(sample_gt_map[sname]))
                        label  = f"JK-MS alpha vs {sname} iter={it}"
                        result = compare_to_gt(tpm_alpha, transcript_names, gt_df, label)
                        all_records.append(result)
                        # For plotting: one alpha series per GT sample
                        key = f"JK-MS alpha vs {sname}"
                        if key not in series_dict:
                            series_dict[key] = {}
                        series_dict[key][it] = result

    return all_records, series_dict


def run_ss_experiment(
    exp_dir:    Path,
    target_iters: list[int],
    sample_gt_map: dict,
) -> tuple[list[dict], dict[str, dict[str, list[dict]]]]:
    """
    Run GT comparison for all target iterations for the JK single-sample experiment.

    For each (sample, iteration): load theta from theta_snapshots.pkl,
    convert to TPM, compare against GT.

    Args:
        exp_dir       : Path -- experiment root directory.
        target_iters  : list[int]  -- iterations to evaluate.
        sample_gt_map : dict       -- {sample_name: gt_path}.

    Returns:
        tuple:
          all_records       : list[dict] -- flat list of compare_to_gt() results.
          series_dict       : dict[series_label, dict[iter → result_dict]]
    """
    all_records = []
    series_dict: dict[str, dict[int, dict]] = {}

    for sname, gt_path in sample_gt_map.items():
        snap_file = exp_dir / sname / "theta_snapshots.pkl"
        if not snap_file.exists():
            print(f"[JK SS] theta_snapshots.pkl not found for {sname}, skipping.")
            continue

        ss_data = load_ss_snapshots(exp_dir, sname)
        gt_df   = read_ground_truth(Path(gt_path))

        abund_path = exp_dir / sname / "abundance.tsv"
        if not abund_path.exists():
            print(f"[JK SS] abundance.tsv not found for {sname}, skipping.")
            continue
        transcript_names, eff_lens = read_eff_lens(abund_path)

        series_key = f"JK-SS {sname}"
        series_dict[series_key] = {}

        for it in target_iters:
            theta = get_ss_snapshot(ss_data, it)
            if theta is None:
                print(f"[JK SS] {sname} iter={it} not found in snapshots, skipping.")
                continue

            tpm_pred = theta_to_tpm(theta, eff_lens)
            label    = f"JK-SS {sname} iter={it}"
            result   = compare_to_gt(tpm_pred, transcript_names, gt_df, label)
            all_records.append(result)
            series_dict[series_key][it] = result

            print(f"  [JK SS] {sname} iter={it}: "
                  f"Spearman(active)={fmt(result['active_spearman'])}, "
                  f"TP={result['n_tp']}, FP={result['n_fp']}, FN={result['n_fn']}")

    return all_records, series_dict


# ============================================================
# Save outputs to one experiment folder
# ============================================================

def save_outputs(
    out_dir:       Path,
    all_records:   list[dict],
    series_dict:   dict[str, dict[int, dict]],
    target_iters:  list[int],
    exp_label:     str,
) -> None:
    """
    Save text reports and figures for one experiment folder.

    Outputs (Rule 2 — multi-file → named subfolder with timestamp):
      {out_dir}/
          metrics_summary.txt
          detail_iter_{round}.txt   (one per target iteration)
          figures/
              metrics_convergence_{timestamp}.png
              alpha_metrics_{timestamp}.png    (if alpha series present)

    Args:
        out_dir      : Path           -- subfolder path (already named + timestamped).
        all_records  : list[dict]     -- flat list of compare_to_gt() results.
        series_dict  : dict           -- {series_label: {iter: result_dict}}.
        target_iters : list[int]      -- iterations in order.
        exp_label    : str            -- experiment short label for titles.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # --- metrics_summary.txt ---
    summary_path = out_dir / "metrics_summary.txt"
    summary_text = build_summary_table(all_records)
    with open(summary_path, "w") as fh:
        fh.write(f"Snapshot GT comparison — {exp_label}\n")
        fh.write(f"Timestamp : {TIMESTAMP}\n")
        fh.write(f"Iterations: {target_iters}\n\n")
        fh.write(summary_text)
    print(f"[save] {summary_path}")

    # --- detail_iter_{round}.txt — one per iteration ---
    # Collect all records for each iteration
    iter_records: dict[int, list[dict]] = {it: [] for it in target_iters}
    for r in all_records:
        # Extract iteration from label (format: "... iter={N}")
        label = r["label"]
        try:
            it_str = label.split("iter=")[-1]
            it_val = int(it_str.split()[0])
        except (ValueError, IndexError):
            continue
        if it_val in iter_records:
            iter_records[it_val].append(r)

    for it, recs in iter_records.items():
        detail_path = out_dir / f"detail_iter_{it}.txt"
        with open(detail_path, "w") as fh:
            fh.write(f"Snapshot GT comparison — {exp_label} — iter={it}\n")
            fh.write(f"Timestamp : {TIMESTAMP}\n\n")
            for r in recs:
                fh.write("=" * 60 + "\n")
                fh.write(format_result_block(r) + "\n\n")
        print(f"[save] {detail_path}")

    # --- Figures: one bar-plot figure per sample (theta), one per alpha series ---
    # Separate theta series (one per sample) from alpha series
    theta_series = {k: v for k, v in series_dict.items() if "alpha" not in k}
    alpha_series = {k: v for k, v in series_dict.items() if "alpha" in k}

    # One figure per sample (theta)
    for series_label, iter_dict in theta_series.items():
        # Derive a safe filename slug from the series label (e.g. "JK-MS sim1" → "sim1")
        slug     = series_label.replace(" ", "_").replace("-", "_").lower()
        fig_path = fig_dir / f"{slug}_{TIMESTAMP}.png"
        _plot_sample_bar_figure(
            sample_name  = series_label,
            iter_dict    = iter_dict,
            target_iters = target_iters,
            title        = f"{exp_label} — {series_label} — metrics vs GT",
            out_path     = fig_path,
        )

    # One figure per alpha series (JK MS only, one per GT sample compared)
    for series_label, iter_dict in alpha_series.items():
        slug          = series_label.replace(" ", "_").replace("-", "_").lower()
        fig_path_alpha = fig_dir / f"{slug}_{TIMESTAMP}.png"
        _plot_sample_bar_figure(
            sample_name  = series_label,
            iter_dict    = iter_dict,
            target_iters = target_iters,
            title        = f"{exp_label} — {series_label} — metrics vs GT",
            out_path     = fig_path_alpha,
        )


# ============================================================
# Main
# ============================================================

def main() -> None:
    """
    Run snapshot GT comparison for JK MS and JK SS experiments.

    Saves text reports and figures to a timestamped subfolder inside
    each experiment folder (Rule 1: both experiments get their own copy).

    Run:
        cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
        conda activate NanoCount_5
        python analysis/run_snapshot_gt_comparison.py
    """
    results_base = Path(RESULTS_BASE)
    ms_exp_dir   = results_base / JK_MS_EXP
    ss_exp_dir   = results_base / JK_SS_EXP

    out_subdir = f"snapshot_gt_comparison_{TIMESTAMP}"

    # ---- JK Multi-sample ----
    print("\n" + "=" * 70)
    print(f"JK Multi-sample: {JK_MS_EXP}")
    print("=" * 70)

    ms_data = load_ms_snapshots(ms_exp_dir)

    ms_records, ms_series = run_ms_experiment(
        exp_dir       = ms_exp_dir,
        ms_data       = ms_data,
        target_iters  = TARGET_ITERS,
        sample_gt_map = SAMPLE_GT_MAP,
    )

    save_outputs(
        out_dir      = ms_exp_dir / out_subdir,
        all_records  = ms_records,
        series_dict  = ms_series,
        target_iters = TARGET_ITERS,
        exp_label    = f"JK-MS ({JK_MS_EXP})",
    )

    # ---- JK Single-sample ----
    print("\n" + "=" * 70)
    print(f"JK Single-sample: {JK_SS_EXP}")
    print("=" * 70)

    ss_records, ss_series = run_ss_experiment(
        exp_dir       = ss_exp_dir,
        target_iters  = TARGET_ITERS,
        sample_gt_map = SAMPLE_GT_MAP,
    )

    save_outputs(
        out_dir      = ss_exp_dir / out_subdir,
        all_records  = ss_records,
        series_dict  = ss_series,
        target_iters = TARGET_ITERS,
        exp_label    = f"JK-SS ({JK_SS_EXP})",
    )

    print("\n[done] Snapshot GT comparison complete.")
    print(f"  JK MS output: {ms_exp_dir / out_subdir}")
    print(f"  JK SS output: {ss_exp_dir / out_subdir}")


if __name__ == "__main__":
    main()
