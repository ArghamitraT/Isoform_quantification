"""
plot_interactive_convergence.py
================================
Interactive convergence explorer: drag a slider to scrub through training rounds
and see how JK single-sample theta, JK MS alpha, and JK MS theta evolve vs GT.

One standalone HTML file per sample (sim1, sim2) — open in any browser, no server
needed. The slider at the bottom maps to snapshot rounds.

Layout (2×2 subplots per page):
  ┌──────────────────────┬──────────────────────┐
  │  Panel 1: GT (static)│  Panel 2: JK θ       │
  ├──────────────────────┼──────────────────────┤
  │  Panel 3: JK MS α    │  Panel 4: JK MS θ    │
  └──────────────────────┴──────────────────────┘

Panels: theta × 1e6 (TPM scale); alpha shown as raw values (not normalised).
Colorscale: RdBu_r (blue=low, red=high — Plotly equivalent of coolwarm).

Inputs:
  JK_SINGLE_DIR : experiment folder from run_joli_kallisto.sh (sim1/ and sim2/ inside)
  JK_MS_DIR     : experiment folder from run_multisample_joli.sh (snapshots.pkl inside)

Outputs (saved inside JK_MS_DIR/figures/):
  interactive_convergence_sim1.html
  interactive_convergence_sim2.html

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/plot_interactive_convergence.py
"""

import os
import pickle
import sys
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ============================================================
# CONFIG — edit before running; do not edit below
# ============================================================
RESULTS_BASE  = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# JK single-sample experiment folder (contains sim1/ and sim2/ subfolders)
JK_SINGLE_DIR = "exprmnt_2026_03_30__22_41_19"

# JK multi-sample experiment folder (contains snapshots.pkl + sim1/ sim2/)
JK_MS_DIR     = "exprmnt_2026_03_30__22_37_55"   # e.g. "exprmnt_2026_03_30__12_00_00"

SAMPLE_NAMES  = ["sim1", "sim2"]

GT_PATHS = {
    "sim1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv",
    "sim2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv",
}

COLORSCALE = "RdBu_r"   # Plotly equivalent of matplotlib coolwarm

# Restrict plot to the top-N GT transcripts by TPM (None = use all)
TOP_N = 50
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


def load_jk_single_snapshots(exp_dir: str, sample_name: str) -> dict:
    """
    Load theta snapshots from a JK single-sample run.

    Args:
        exp_dir     : str -- absolute path to experiment folder.
        sample_name : str -- sample subfolder name.

    Returns:
        dict with keys "transcript_names" and "snapshots" (list of (round, theta)).
    """
    path = os.path.join(exp_dir, sample_name, "theta_snapshots.pkl")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"JK single snapshot not found: {path}\n"
            "Re-run run_joli_kallisto.sh with SAVE_SNAPSHOTS=true."
        )
    with open(path, "rb") as fh:
        data = pickle.load(fh)
    print(f"  JK single {sample_name}: {len(data['snapshots'])} snapshots")
    return data


def load_jk_ms_snapshots(exp_dir: str) -> dict:
    """
    Load alpha + theta snapshots from a JK MS run.

    Args:
        exp_dir : str -- absolute path to JK MS experiment folder.

    Returns:
        dict with keys "sample_names" and "snapshots"
        (list of {"round", "alpha", "thetas"}).
    """
    path = os.path.join(exp_dir, "snapshots.pkl")
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"JK MS snapshot not found: {path}\n"
            "Re-run run_multisample_joli.sh with SAVE_SNAPSHOTS=true."
        )
    with open(path, "rb") as fh:
        data = pickle.load(fh)
    print(f"  JK MS: {len(data['snapshots'])} snapshots")
    return data


def load_transcript_names(exp_dir: str, sample_name: str) -> list:
    """
    Load transcript names from a sample's abundance.tsv.

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

def build_universe(gt_df: pd.DataFrame,
                   top_n: int | None = None) -> tuple:
    """
    Build sorted transcript universe purely from ground truth.

    Universe = all GT non-zero transcripts, sorted by GT TPM descending.
    The x-axis rank is determined solely by GT TPM — methods that did not
    predict a transcript get 0; they do not exclude transcripts from the axis.

    Args:
        gt_df : pd.DataFrame -- [transcript_id, tpm_gt].
        top_n : int | None   -- if set, restrict to top-N GT transcripts by TPM.

    Returns:
        (universe: pd.DataFrame with columns [transcript_id, tpm_gt, rank], n_tp: int)
    """
    universe = (gt_df[gt_df["tpm_gt"] > 0]
                .sort_values("tpm_gt", ascending=False)
                .reset_index(drop=True))

    if top_n is not None:
        universe = universe.head(top_n).reset_index(drop=True)

    universe["rank"] = np.arange(len(universe))
    n_tp = len(universe)
    print(f"  Universe: {n_tp} GT transcripts{f' (top {top_n} by GT TPM)' if top_n else ''}")
    return universe, n_tp


# ============================================================
# Per-panel colorscale bounds (stable across all frames)
# ============================================================

def panel_bounds(all_vals: list) -> tuple:
    """
    Compute stable vmin/vmax across all frames using 1st–99th percentile of nonzeros.

    Args:
        all_vals : list[np.ndarray] -- one array per frame.

    Returns:
        (vmin: float, vmax: float)
    """
    stacked = np.concatenate([v[v > 0] for v in all_vals if (v > 0).any()])
    if len(stacked) == 0:
        return 0.0, 1.0
    # Use actual max (not 99th percentile) so the full alpha range is visible.
    return float(np.percentile(stacked, 1)), float(stacked.max())


# ============================================================
# HTML builder
# ============================================================

def build_interactive_html(
    sample_name:    str,
    universe:       pd.DataFrame,
    n_tp:           int,
    gt_vals:        np.ndarray,
    gt_tpm_vmin:    float,
    gt_tpm_vmax:    float,
    frames_data:    list,   # list of {"round": int, "jk": arr, "alpha": arr, "jkms": arr}
    output_path:    str,
) -> None:
    """
    Build and save a Plotly interactive HTML with a round slider.

    Args:
        sample_name   : str           -- e.g. "sim1".
        universe      : pd.DataFrame  -- sorted universe with [transcript_id, rank, tpm_gt].
        n_tp          : int           -- number of TP transcripts.
        gt_vals       : np.ndarray    -- GT TPM per universe position.
        gt_tpm_vmin   : float         -- GT colorscale min.
        gt_tpm_vmax   : float         -- GT colorscale max.
        frames_data   : list[dict]    -- per-round values for panels 2/3/4.
        output_path   : str           -- path to save HTML.
    """
    ranks = universe["rank"].values
    sep   = n_tp - 0.5   # TP/FP boundary position

    # Theta panels (JK θ and JK MS θ) share the same colorscale as GT so they are
    # directly comparable. GT bounds are passed in as gt_tpm_vmin/vmax.
    jk_vmin,   jk_vmax   = gt_tpm_vmin, gt_tpm_vmax
    jkms_vmin, jkms_vmax = gt_tpm_vmin, gt_tpm_vmax

    # Alpha uses its own scale: raw Dirichlet concentration, which can be large.
    # Use the global max across all frames so the colorscale is stable.
    alpha_vmin, alpha_vmax = panel_bounds([f["alpha"] for f in frames_data])

    # ---- Helper: build one scatter trace ----
    def make_scatter(x, y, vals, vmin, vmax, row, col, name, showscale=True,
                     colorbar_title="TPM"):
        return go.Scatter(
            x=x, y=y,
            mode="markers",
            marker=dict(
                color=vals,
                colorscale=COLORSCALE,
                cmin=vmin,
                cmax=vmax,
                size=3,
                showscale=showscale,
                colorbar=dict(
                    len=0.45,
                    thickness=12,
                    x=0.48 if col == 1 else 1.01,
                    y=0.78 if row == 1 else 0.22,
                    title=dict(text=colorbar_title, font=dict(size=9)),
                    tickfont=dict(size=8),
                ),
            ),
            name=name,
            showlegend=False,
            xaxis=f"x{(row-1)*2 + col}",
            yaxis=f"y{(row-1)*2 + col}",
        )

    # ---- Static GT trace (panel 1, row=1 col=1) ----
    gt_trace = make_scatter(
        ranks, gt_vals, gt_vals, gt_tpm_vmin, gt_tpm_vmax,
        row=1, col=1, name="GT"
    )

    # ---- Initial frame traces (round 0) ----
    f0 = frames_data[0]
    jk_trace    = make_scatter(ranks, f0["jk"],    f0["jk"],    jk_vmin,    jk_vmax,    row=1, col=2, name="JK θ")
    alpha_trace = make_scatter(ranks, f0["alpha"], f0["alpha"], alpha_vmin, alpha_vmax, row=2, col=1, name="JK MS α",
                               colorbar_title="α value")
    jkms_trace  = make_scatter(ranks, f0["jkms"],  f0["jkms"],  jkms_vmin,  jkms_vmax,  row=2, col=2, name="JK MS θ")

    # Helper: compute linear y-axis range [0, max * 1.05] for a data array.
    # Used to update per-panel y-axis in each frame so the axis follows the data.
    def yrange(vals, pad=0.05):
        nz = vals[vals > 0]
        if len(nz) == 0:
            return [0, 1]
        return [0, float(vals.max()) * (1 + pad)]

    # ---- Plotly frames (one per snapshot round) ----
    plotly_frames = []
    for fd in frames_data:
        rnd = fd["round"]
        # All universe entries have GT > 0 (universe is purely GT-driven), so compute
        # Spearman directly over the full universe — GT values vs predicted values.
        # Predicted values are 0 for transcripts not in the method's index.
        jk_sp = float(spearmanr(gt_vals, fd["jk"]).statistic
                      if len(gt_vals) > 1 else float("nan"))
        ms_sp = float(spearmanr(gt_vals, fd["jkms"]).statistic
                      if len(gt_vals) > 1 else float("nan"))

        # Per-frame colorscale bounds: let each panel's color track its own data range
        jk_cmax    = float(fd["jk"].max())    if fd["jk"].max()    > 0 else 1.0
        alpha_cmax = float(fd["alpha"].max()) if fd["alpha"].max() > 0 else 1.0
        jkms_cmax  = float(fd["jkms"].max())  if fd["jkms"].max()  > 0 else 1.0

        frame = go.Frame(
            name=str(rnd),
            data=[
                # GT trace stays the same — reuse it unchanged (trace index 0)
                go.Scatter(x=ranks, y=gt_vals,
                           marker=dict(color=gt_vals, colorscale=COLORSCALE,
                                       cmin=gt_tpm_vmin, cmax=gt_tpm_vmax, size=3,
                                       showscale=True,
                                       colorbar=dict(len=0.45, thickness=12,
                                                     x=0.48, y=0.78))),
                go.Scatter(x=ranks, y=fd["jk"],
                           marker=dict(color=fd["jk"], colorscale=COLORSCALE,
                                       cmin=gt_tpm_vmin, cmax=gt_tpm_vmax, size=3,
                                       showscale=True,
                                       colorbar=dict(len=0.45, thickness=12,
                                                     x=1.01, y=0.78,
                                                     title=dict(text="TPM", font=dict(size=9)),
                                                     tickfont=dict(size=8)))),
                go.Scatter(x=ranks, y=fd["alpha"],
                           marker=dict(color=fd["alpha"], colorscale=COLORSCALE,
                                       cmin=0, cmax=alpha_cmax, size=3,
                                       showscale=True,
                                       colorbar=dict(len=0.45, thickness=12,
                                                     x=0.48, y=0.22,
                                                     title=dict(text="α value", font=dict(size=9)),
                                                     tickfont=dict(size=8)))),
                go.Scatter(x=ranks, y=fd["jkms"],
                           marker=dict(color=fd["jkms"], colorscale=COLORSCALE,
                                       cmin=gt_tpm_vmin, cmax=gt_tpm_vmax, size=3,
                                       showscale=True,
                                       colorbar=dict(len=0.45, thickness=12,
                                                     x=1.01, y=0.22,
                                                     title=dict(text="TPM", font=dict(size=9)),
                                                     tickfont=dict(size=8)))),
            ],
            traces=[0, 1, 2, 3],
            layout=go.Layout(
                title=dict(
                    text=(f"{sample_name}  —  Round {rnd}  |  "
                          f"Spearman vs GT:  JK={jk_sp:.4f}   JK MS={ms_sp:.4f}"),
                ),
                # Update y-axis ranges to follow per-round data values
                yaxis2=dict(range=yrange(fd["jk"])),
                yaxis3=dict(range=yrange(fd["alpha"])),
                yaxis4=dict(range=yrange(fd["jkms"])),
            ),
        )
        plotly_frames.append(frame)

    # ---- Slider steps ----
    steps = []
    for fd in frames_data:
        rnd = fd["round"]
        steps.append(dict(
            method="animate",
            label=str(rnd),
            args=[[str(rnd)], dict(
                mode="immediate",
                frame=dict(duration=0, redraw=True),
                transition=dict(duration=0),
            )],
        ))

    slider = dict(
        active=0,
        steps=steps,
        currentvalue=dict(prefix="Round: ", font=dict(size=13)),
        pad=dict(t=50),
        len=0.9,
        x=0.05,
    )

    # ---- TP/FP separator lines (added as shapes) ----
    shapes = []
    if n_tp < len(ranks):
        for xref, yref in [("x1", "y1"), ("x2", "y2"), ("x3", "y3"), ("x4", "y4")]:
            shapes.append(dict(
                type="line",
                x0=sep, x1=sep, y0=0, y1=1,
                xref=xref, yref="paper",
                line=dict(color="black", width=1, dash="dash"),
            ))

    # ---- Subplot layout ----
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            "GT  (reference TPM)",
            "JK single-sample  θ × 1e6",
            "JK MS  α (raw)",
            "JK MS  θ × 1e6",
        ],
        horizontal_spacing=0.12,
        vertical_spacing=0.15,
    )

    for trace, row, col in [
        (gt_trace,    1, 1),
        (jk_trace,    1, 2),
        (alpha_trace, 2, 1),
        (jkms_trace,  2, 2),
    ]:
        fig.add_trace(trace, row=row, col=col)

    # Add TP/FP separator annotation text
    if n_tp < len(ranks):
        for row, col in [(1,1),(1,2),(2,1),(2,2)]:
            fig.add_annotation(
                x=sep, y=1.02, xref=f"x{(row-1)*2+col}", yref="paper",
                text="← TP | FP →", showarrow=False,
                font=dict(size=8, color="black"),
            )

    # Initial y-axis ranges — set from first frame so the axes start correct
    f0 = frames_data[0]
    fig.update_layout(
        title=dict(
            text=(f"{sample_name}  —  Round {f0['round']}  |  "
                  f"(drag slider below to scrub through rounds)"),
            font=dict(size=13),
        ),
        height=750,
        sliders=[slider],
        shapes=shapes,
        margin=dict(t=80, b=100),
        yaxis2=dict(range=yrange(f0["jk"])),
        yaxis3=dict(range=yrange(f0["alpha"])),
        yaxis4=dict(range=yrange(f0["jkms"])),
    )

    # Log-scale y-axis for GT panel (panel 1) — static reference, never changes
    gt_nz_vals = gt_vals[gt_vals > 0]
    fig.update_yaxes(
        type="log",
        range=[np.log10(gt_nz_vals.min()) - 0.2, np.log10(gt_nz_vals.max()) + 0.2],
        row=1, col=1,
    )

    fig.frames = plotly_frames

    fig.write_html(output_path, include_plotlyjs="cdn")
    print(f"  Saved interactive HTML: {output_path}")


# ============================================================
# Main
# ============================================================

def main() -> None:
    """
    Load snapshots, build universe per sample, render interactive HTML per sample.
    """
    for name, val in [("JK_SINGLE_DIR", JK_SINGLE_DIR), ("JK_MS_DIR", JK_MS_DIR)]:
        if not val:
            print(f"ERROR: {name} is empty — set it in the CONFIG section.")
            sys.exit(1)

    jk_dir = os.path.join(RESULTS_BASE, JK_SINGLE_DIR)
    ms_dir = os.path.join(RESULTS_BASE, JK_MS_DIR)

    figures_dir = os.path.join(ms_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)

    print("=" * 60)
    print("plot_interactive_convergence.py")
    print("=" * 60)
    print(f"  JK single dir    : {JK_SINGLE_DIR}")
    print(f"  JK MS dir        : {JK_MS_DIR}")
    print(f"  Output figures   : {figures_dir}")
    print()

    # Load GT
    gt_dfs = {s: load_gt(GT_PATHS[s]) for s in SAMPLE_NAMES}

    # Load snapshots
    jk_snaps   = {s: load_jk_single_snapshots(jk_dir, s) for s in SAMPLE_NAMES}
    ms_snaps   = load_jk_ms_snapshots(ms_dir)
    ms_snames  = ms_snaps["sample_names"]
    ms_sidx    = {s: ms_snames.index(s) for s in SAMPLE_NAMES if s in ms_snames}

    # Transcript index from JK MS (shared reference)
    tx_names_ms = load_transcript_names(ms_dir, ms_snames[0])
    tx_index_ms = {tid: i for i, tid in enumerate(tx_names_ms)}

    # JK single transcript index (may differ per sample)
    tx_index_jk = {
        s: {tid: i for i, tid in enumerate(jk_snaps[s]["transcript_names"])}
        for s in SAMPLE_NAMES
    }

    # Round index lookups
    jk_ridx = {
        s: {r: i for i, (r, _) in enumerate(jk_snaps[s]["snapshots"])}
        for s in SAMPLE_NAMES
    }
    ms_ridx = {snap["round"]: i for i, snap in enumerate(ms_snaps["snapshots"])}

    # Use JK MS rounds as the reference frame set
    common_rounds = sorted(ms_ridx.keys())

    for sample in SAMPLE_NAMES:
        print(f"\n--- {sample} ---")
        s_idx  = ms_sidx.get(sample, 0)
        gt_df  = gt_dfs[sample]

        # Universe is built purely from GT — no final-round theta needed here.
        universe, n_tp = build_universe(gt_df, top_n=TOP_N)

        tx_ids  = universe["transcript_id"].tolist()
        gt_vals = universe["tpm_gt"].values

        # GT colorscale bounds
        gt_nz      = gt_vals[gt_vals > 0]
        gt_tpm_vmin = float(np.percentile(gt_nz, 1))  if len(gt_nz) > 0 else 0.0
        gt_tpm_vmax = float(np.percentile(gt_nz, 99)) if len(gt_nz) > 0 else 1.0

        # Build per-round data for panels 2/3/4
        # Filter to rounds present in both JK single and JK MS
        avail_rounds = [r for r in common_rounds
                        if r in jk_ridx[sample] and r in ms_ridx]

        frames_data = []
        for rnd in avail_rounds:
            _, jk_theta = jk_snaps[sample]["snapshots"][jk_ridx[sample][rnd]]
            ms_snap     = ms_snaps["snapshots"][ms_ridx[rnd]]
            alpha_raw   = ms_snap["alpha"]
            jkms_theta  = ms_snap["thetas"][s_idx]

            # Re-index JK single theta to universe order
            jk_vals = np.array([
                jk_theta[tx_index_jk[sample].get(tid, 0)]
                if tid in tx_index_jk[sample] else 0.0
                for tid in tx_ids
            ]) * 1e6   # → pseudo-TPM

            # Alpha: raw values (no normalization).
            # Use explicit membership guard: .get(tid, 0) would default to INDEX 0
            # (the first transcript's value), not the value 0.0.
            alpha_u   = np.array([
                alpha_raw[tx_index_ms[tid]] if tid in tx_index_ms else 0.0
                for tid in tx_ids
            ])
            alpha_tpm = alpha_u

            # JK MS theta ×1e6 — same guard as above
            jkms_vals = np.array([
                jkms_theta[tx_index_ms[tid]] if tid in tx_index_ms else 0.0
                for tid in tx_ids
            ]) * 1e6

            frames_data.append({
                "round": rnd,
                "jk":    jk_vals,
                "alpha": alpha_tpm,
                "jkms":  jkms_vals,
            })

        # --- Diagnostic: show top-5 GT transcripts with final-round values ---
        # Helps explain any 0 values in the MS panel (e.g. method didn't index
        # the transcript, or EM genuinely assigned 0 at that snapshot round).
        if frames_data:
            fd_last = frames_data[-1]
            n_diag = min(5, n_tp)
            n_missing_jk   = sum(1 for tid in tx_ids[:n_tp] if tid not in tx_index_jk[sample])
            n_missing_jkms = sum(1 for tid in tx_ids[:n_tp] if tid not in tx_index_ms)
            print(f"\n  Diagnostic — top {n_diag} GT transcripts (round {fd_last['round']}):")
            print(f"  {'Rank':<6} {'transcript_id':<35} {'GT_TPM':>10} {'JK_theta*1e6':>14} {'JKMS_theta*1e6':>16}")
            for i in range(n_diag):
                print(f"  {i:<6} {tx_ids[i]:<35} {gt_vals[i]:>10.2f} "
                      f"{fd_last['jk'][i]:>14.2f} {fd_last['jkms'][i]:>16.2f}")
            print(f"  Transcripts in universe not in JK index  : {n_missing_jk}")
            print(f"  Transcripts in universe not in JKMS index: {n_missing_jkms}")

        print(f"  Building HTML with {len(frames_data)} frames...")
        timestamp  = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")
        top_suffix = f"_top{TOP_N}" if TOP_N is not None else ""
        out_path = os.path.join(figures_dir,
                                f"interactive_convergence{top_suffix}_{sample}_{timestamp}.html")
        build_interactive_html(
            sample_name  = sample,
            universe     = universe,
            n_tp         = n_tp,
            gt_vals      = gt_vals,
            gt_tpm_vmin  = gt_tpm_vmin,
            gt_tpm_vmax  = gt_tpm_vmax,
            frames_data  = frames_data,
            output_path  = out_path,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
