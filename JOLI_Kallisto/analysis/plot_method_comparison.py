"""
plot_method_comparison.py
=========================
Bar-chart comparison of LK, JK_SS, and JK_MS across all simulation samples.

Produces 3 figures:
  1. contingency_{timestamp}.png  — TP / FP / FN / TN
  2. active_{timestamp}.png       — Active-universe metrics
  3. gt_nonzero_{timestamp}.png   — GT-non-zero metrics

Layout of each figure:
  rows : samples   (sim1/LR, sim2/LR, sim_sr_s1/SR, sim_sr_s2/SR)
  cols : metrics   (varies per figure)
  bars : methods   (LK=blue, JK_SS=orange, JK_MS=green)

Outputs saved to figures/ inside every experiment folder involved
(following analysis output convention Rule 3 — 1-2 files → figures/ directly).

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/plot_method_comparison.py
"""

import os
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# CONFIG — edit experiment folder lists if you rerun experiments
# ============================================================

RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

# All folders that will receive a copy of the figures
ALL_FOLDERS = [
    "exprmnt_2026_07_18__23_43_13",   # LK LR
    "exprmnt_2026_07_18__23_48_54",   # LK SR
    "exprmnt_2026_07_18__23_45_57",   # JK_SS LR
    "exprmnt_2026_07_19__00_08_14",   # JK_SS SR
    "exprmnt_2026_07_19__00_21_20",   # JK_MS run1
    "exprmnt_2026_07_19__00_26_49",   # JK_MS run2
]

# ============================================================
# DATA — extracted from gt_comparison files
# Structure: DATA[sample][method] = {contingency, active, gt_nz}
# ============================================================

DATA = {
    "sim1\n(LR)": {
        "LK": {
            "contingency": {"TP": 62199, "FP": 6092,  "FN": 8956,  "TN": 0},
            "active":  {"Spearman": 0.770109, "Pearson": 0.966390, "MAE": 3.658319, "RMSE": 24.248716, "Bray-Curtis": 0.141297},
            "gt_nz":   {"Spearman": 0.823903, "Pearson": 0.970511, "MAE": 3.587216, "RMSE": 23.720277, "Bray-Curtis": 0.129393},
        },
        "JK_SS": {
            "contingency": {"TP": 62186, "FP": 6084,  "FN": 8969,  "TN": 141016},
            "active":  {"Spearman": 0.770096, "Pearson": 0.966390, "MAE": 3.658698, "RMSE": 24.249970, "Bray-Curtis": 0.141297},
            "gt_nz":   {"Spearman": 0.823906, "Pearson": 0.970511, "MAE": 3.587217, "RMSE": 23.720275, "Bray-Curtis": 0.129393},
        },
        "JK_MS": {
            "contingency": {"TP": 69285, "FP": 16066, "FN": 1870,  "TN": 131034},
            "active":  {"Spearman": 0.846255, "Pearson": 0.961182, "MAE": 3.244403, "RMSE": 24.875246, "Bray-Curtis": 0.141490},
            "gt_nz":   {"Spearman": 0.839724, "Pearson": 0.961955, "MAE": 3.738431, "RMSE": 27.205607, "Bray-Curtis": 0.134142},
        },
    },
    "sim2\n(LR)": {
        "LK": {
            "contingency": {"TP": 69144, "FP": 5854,  "FN": 9833,  "TN": 0},
            "active":  {"Spearman": 0.775123, "Pearson": 0.944948, "MAE": 3.483899, "RMSE": 27.848823, "Bray-Curtis": 0.147771},
            "gt_nz":   {"Spearman": 0.824680, "Pearson": 0.949769, "MAE": 3.402283, "RMSE": 27.617673, "Bray-Curtis": 0.136179},
        },
        "JK_SS": {
            "contingency": {"TP": 69139, "FP": 5850,  "FN": 9838,  "TN": 133428},
            "active":  {"Spearman": 0.775128, "Pearson": 0.944948, "MAE": 3.484063, "RMSE": 27.849474, "Bray-Curtis": 0.147771},
            "gt_nz":   {"Spearman": 0.824693, "Pearson": 0.949769, "MAE": 3.402283, "RMSE": 27.617666, "Bray-Curtis": 0.136179},
        },
        "JK_MS": {
            "contingency": {"TP": 76878, "FP": 14939, "FN": 2099,  "TN": 124339},
            "active":  {"Spearman": 0.844060, "Pearson": 0.956019, "MAE": 3.072466, "RMSE": 24.065583, "Bray-Curtis": 0.144277},
            "gt_nz":   {"Spearman": 0.839697, "Pearson": 0.956742, "MAE": 3.464301, "RMSE": 25.976593, "Bray-Curtis": 0.137831},
        },
    },
    "sim_sr_s1\n(SR)": {
        "LK": {
            "contingency": {"TP": 81096,  "FP": 10085, "FN": 34238, "TN": 92836},
            "active":  {"Spearman": 0.644334, "Pearson": 0.983264, "MAE": 1.692896, "RMSE": 13.778295, "Bray-Curtis": 0.106161},
            "gt_nz":   {"Spearman": 0.743135, "Pearson": 0.987396, "MAE": 1.597570, "RMSE": 12.736116, "Bray-Curtis": 0.093438},
        },
        "JK_SS": {
            "contingency": {"TP": 78361,  "FP": 9419,  "FN": 36973, "TN": 93502},
            "active":  {"Spearman": 0.645146, "Pearson": 0.983263, "MAE": 1.701944, "RMSE": 13.815013, "Bray-Curtis": 0.106161},
            "gt_nz":   {"Spearman": 0.742936, "Pearson": 0.987396, "MAE": 1.597577, "RMSE": 12.736090, "Bray-Curtis": 0.093439},
        },
        "JK_MS": {
            "contingency": {"TP": 115021, "FP": 44934, "FN": 313,   "TN": 57987},
            "active":  {"Spearman": 0.838207, "Pearson": 0.983130, "MAE": 1.140292, "RMSE": 12.319430, "Bray-Curtis": 0.091376},
            "gt_nz":   {"Spearman": 0.899464, "Pearson": 0.987011, "MAE": 1.317617, "RMSE": 13.010904, "Bray-Curtis": 0.077171},
        },
    },
    "sim_sr_s2\n(SR)": {
        "LK": {
            "contingency": {"TP": 87987,  "FP": 9641,  "FN": 32425, "TN": 88202},
            "active":  {"Spearman": 0.663420, "Pearson": 0.980212, "MAE": 1.720256, "RMSE": 13.706350, "Bray-Curtis": 0.111862},
            "gt_nz":   {"Spearman": 0.750147, "Pearson": 0.985182, "MAE": 1.625917, "RMSE": 12.543100, "Bray-Curtis": 0.099277},
        },
        "JK_SS": {
            "contingency": {"TP": 84686,  "FP": 8951,  "FN": 35726, "TN": 88892},
            "active":  {"Spearman": 0.663388, "Pearson": 0.980210, "MAE": 1.729430, "RMSE": 13.742835, "Bray-Curtis": 0.111862},
            "gt_nz":   {"Spearman": 0.749373, "Pearson": 0.985182, "MAE": 1.625917, "RMSE": 12.543070, "Bray-Curtis": 0.099277},
        },
        "JK_MS": {
            "contingency": {"TP": 120081, "FP": 42433, "FN": 331,   "TN": 55410},
            "active":  {"Spearman": 0.841523, "Pearson": 0.980515, "MAE": 1.139878, "RMSE": 12.200808, "Bray-Curtis": 0.092812},
            "gt_nz":   {"Spearman": 0.895516, "Pearson": 0.985253, "MAE": 1.282341, "RMSE": 12.554666, "Bray-Curtis": 0.078429},
        },
    },
}

# ============================================================
# Plot config
# ============================================================

METHODS       = ["LK", "JK_SS", "JK_MS"]
COLORS        = {"LK": "#cccccc", "JK_SS": "#888888", "JK_MS": "#333333"}   # light→dark grey
SAMPLES       = list(DATA.keys())

# Samples that were fed together as one JK_MS joint run get the same highlight color
# on their y-axis label, so the actual MS pairing is visible at a glance:
#   Run1 = sim1(LR)      + sim_sr_s2(SR)  -> pink
#   Run2 = sim2(LR)      + sim_sr_s1(SR)  -> blue
PAIR_COLORS = {
    "sim1\n(LR)":      "#ff69b4",  # pink  — JK_MS Run1
    "sim_sr_s2\n(SR)": "#ff69b4",  # pink  — JK_MS Run1
    "sim2\n(LR)":      "#5b9bd5",  # blue  — JK_MS Run2
    "sim_sr_s1\n(SR)": "#5b9bd5",  # blue  — JK_MS Run2
}


def _highlight_ylabel(ax, sample: str, fontsize: int = 8, labelpad: int = 6) -> None:
    """
    Set the row's y-axis label with a colored highlight box matching its
    JK_MS input pairing (see PAIR_COLORS).

    Args:
        ax       : matplotlib Axes
        sample   : str  -- sample key, e.g. "sim1\\n(LR)"
        fontsize : int
        labelpad : int

    Returns:
        None
    """
    color = PAIR_COLORS.get(sample)
    kwargs = dict(fontsize=fontsize, fontweight="bold", labelpad=labelpad)
    if color:
        kwargs["bbox"] = dict(facecolor=color, alpha=0.35, edgecolor=color,
                               boxstyle="round,pad=0.35")
    ax.set_ylabel(sample.replace("\n", " "), **kwargs)

# FP TPM% = sum(predicted TPM for FP transcripts) / sum(GT TPM) * 100
# Computed from abundance.tsv vs ground truth; keyed to same sample labels as DATA
FP_TPM_PCT = {
    "sim1\n(LR)":      {"LK": 2.73, "JK_SS": 2.73, "JK_MS": 1.70},
    "sim2\n(LR)":      {"LK": 2.68, "JK_SS": 2.68, "JK_MS": 1.50},
    "sim_sr_s1\n(SR)": {"LK": 2.81, "JK_SS": 2.81, "JK_MS": 3.08},
    "sim_sr_s2\n(SR)": {"LK": 2.79, "JK_SS": 2.79, "JK_MS": 3.12},
}

CONTINGENCY_METRICS = ["TP", "FP", "FN", "TN"]
PERF_METRICS        = ["Spearman", "Pearson", "MAE", "RMSE", "Bray-Curtis"]

# Metrics where lower = better (winner gets a downward arrow annotation)
LOWER_BETTER = {"MAE", "RMSE", "Bray-Curtis", "FP", "FN"}


# ============================================================
# Helpers
# ============================================================

def _bar_group(ax, values: dict, title: str, ylabel: str = "") -> None:
    """
    Draw a grouped bar chart with one group per method.

    Args:
        ax     : matplotlib Axes
        values : dict method → float
        title  : str -- subplot title
        ylabel : str -- y-axis label
    """
    x     = np.arange(len(METHODS))
    width = 0.55
    bars  = [values[m] for m in METHODS]
    cols  = [COLORS[m] for m in METHODS]

    rects = ax.bar(x, bars, width, color=cols, edgecolor="white", linewidth=0.5)

    # Mark winner
    lower = title in LOWER_BETTER
    wi    = int(np.argmin(bars) if lower else np.argmax(bars))
    ax.bar(x[wi:wi+1], bars[wi:wi+1], width,
           color=cols[wi], edgecolor="red", linewidth=2.0)

    ax.set_title(title, fontsize=8, pad=3)
    ax.set_xticks(x)
    ax.set_xticklabels(METHODS, fontsize=7, rotation=15, ha="right")
    ax.tick_params(axis="y", labelsize=7)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def make_figure(set_key: str, metric_keys: list, title: str, fname_prefix: str,
                timestamp: str) -> plt.Figure:
    """
    Build a figure: rows=samples, cols=metrics, bars=methods.

    Args:
        set_key      : str  -- data key: "contingency" | "active" | "gt_nz"
        metric_keys  : list -- metric names to plot (columns)
        title        : str  -- figure suptitle
        fname_prefix : str  -- filename prefix before timestamp
        timestamp    : str  -- datetime string

    Returns:
        plt.Figure
    """
    n_rows = len(SAMPLES)
    n_cols = len(metric_keys)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(2.8 * n_cols, 3.0 * n_rows),
        squeeze=False,
    )
    fig.suptitle(title, fontsize=12, fontweight="bold", y=1.01)

    for r, sample in enumerate(SAMPLES):
        for c, metric in enumerate(metric_keys):
            ax = axes[r][c]
            values = {m: DATA[sample][m][set_key][metric] for m in METHODS}

            # Row label on leftmost column
            ylabel = sample.replace("\n", " ") if c == 0 else ""
            _bar_group(ax, values, metric, ylabel=ylabel)

            # Sample label on the left axis, highlighted by JK_MS input pairing
            if c == 0:
                _highlight_ylabel(ax, sample)

    # Legends: methods + JK_MS pairing
    method_handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[m], label=m) for m in METHODS]
    fig.legend(handles=method_handles, loc="upper right", fontsize=9,
               framealpha=0.9, title="Method", title_fontsize=9)

    pair_handles = [
        plt.Rectangle((0, 0), 1, 1, color="#ff69b4", alpha=0.35, label="Run1: sim1(LR)+sim_sr_s2(SR)"),
        plt.Rectangle((0, 0), 1, 1, color="#5b9bd5", alpha=0.35, label="Run2: sim2(LR)+sim_sr_s1(SR)"),
    ]
    fig.legend(handles=pair_handles, loc="upper left", fontsize=8,
               framealpha=0.9, title="JK_MS pairing", title_fontsize=8)

    fig.tight_layout()
    return fig


# ============================================================
# FP-rate figure
# ============================================================

def make_fp_figure(timestamp: str) -> plt.Figure:
    """
    Build a 4-row × 3-col figure showing, per sample:
      Col 1 — Total predicted non-zero transcripts (TP+FP), stacked TP/FP bars
      Col 2 — FP% = FP / (TP+FP) * 100  (biased: larger predicted set → smaller %)
      Col 3 — θ_pred_FP / θ_GT_total * 100  (fair: expression-weighted FP vs GT budget)

    Args:
        timestamp : str -- datetime string for the title.

    Returns:
        plt.Figure
    """
    n_rows = len(SAMPLES)
    fig, axes = plt.subplots(n_rows, 3, figsize=(10.5, 3.0 * n_rows), squeeze=False)
    fig.suptitle(
        "False Positive Analysis — LK / JK_SS / JK_MS",
        fontsize=12, fontweight="bold", y=1.01,
    )

    x        = np.arange(len(METHODS))
    width    = 0.5
    tp_color = "#aaaaaa"   # grey for TP portion (matches greyscale theme)
    fp_color = "#e05c5c"   # red for FP portion
    bar_cols = [COLORS[m] for m in METHODS]

    for r, sample in enumerate(SAMPLES):
        ax_stack = axes[r][0]
        ax_pct   = axes[r][1]
        ax_tpm   = axes[r][2]

        tps  = [DATA[sample][m]["contingency"]["TP"] for m in METHODS]
        fps  = [DATA[sample][m]["contingency"]["FP"] for m in METHODS]
        tots = [tp + fp for tp, fp in zip(tps, fps)]
        pcts = [100 * fp / tot if tot > 0 else 0 for fp, tot in zip(fps, tots)]
        tpm_pcts = [FP_TPM_PCT[sample][m] for m in METHODS]

        # --- Col 1: stacked bar (TP bottom, FP on top) ---
        ax_stack.bar(x, tps, width, color=tp_color, label="TP")
        ax_stack.bar(x, fps, width, bottom=tps, color=fp_color, label="FP")
        ax_stack.set_xticks(x)
        ax_stack.set_xticklabels(METHODS, fontsize=8, rotation=15, ha="right")
        _highlight_ylabel(ax_stack, sample, fontsize=9, labelpad=4)
        ax_stack.tick_params(axis="y", labelsize=7)
        ax_stack.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v/1000)}k")
        )
        ax_stack.spines["top"].set_visible(False)
        ax_stack.spines["right"].set_visible(False)
        if r == 0:
            ax_stack.set_title("Total predicted\n(TP + FP)", fontsize=9)
            ax_stack.legend(fontsize=7, loc="upper left")

        # --- Col 2: FP% = FP / (TP+FP) bar ---
        wi = int(np.argmin(pcts))
        rects2 = ax_pct.bar(x, pcts, width, color=bar_cols, edgecolor="white", linewidth=0.5)
        ax_pct.bar(x[wi:wi+1], pcts[wi:wi+1], width,
                   color=bar_cols[wi], edgecolor="red", linewidth=2.0)
        for rect, pct in zip(rects2, pcts):
            ax_pct.text(rect.get_x() + rect.get_width() / 2,
                        rect.get_height() + 0.3,
                        f"{pct:.1f}%", ha="center", va="bottom", fontsize=8, fontweight="bold")
        ax_pct.set_xticks(x)
        ax_pct.set_xticklabels(METHODS, fontsize=8, rotation=15, ha="right")
        ax_pct.set_ylabel("FP%", fontsize=8)
        ax_pct.tick_params(axis="y", labelsize=7)
        ax_pct.set_ylim(0, max(pcts) * 1.3)
        ax_pct.spines["top"].set_visible(False)
        ax_pct.spines["right"].set_visible(False)
        if r == 0:
            ax_pct.set_title("FP / (TP+FP) × 100\n(biased by predicted set size)", fontsize=8)

        # --- Col 3: θ_pred_FP / θ_GT_total bar ---
        wi3 = int(np.argmin(tpm_pcts))
        rects3 = ax_tpm.bar(x, tpm_pcts, width, color=bar_cols, edgecolor="white", linewidth=0.5)
        ax_tpm.bar(x[wi3:wi3+1], tpm_pcts[wi3:wi3+1], width,
                   color=bar_cols[wi3], edgecolor="red", linewidth=2.0)
        for rect, pct in zip(rects3, tpm_pcts):
            ax_tpm.text(rect.get_x() + rect.get_width() / 2,
                        rect.get_height() + 0.02,
                        f"{pct:.2f}%", ha="center", va="bottom", fontsize=8, fontweight="bold")
        ax_tpm.set_xticks(x)
        ax_tpm.set_xticklabels(METHODS, fontsize=8, rotation=15, ha="right")
        ax_tpm.set_ylabel("θ_pred_FP / θ_GT_total (%)", fontsize=8)
        ax_tpm.tick_params(axis="y", labelsize=7)
        ax_tpm.set_ylim(0, max(tpm_pcts) * 1.3)
        ax_tpm.spines["top"].set_visible(False)
        ax_tpm.spines["right"].set_visible(False)
        if r == 0:
            ax_tpm.set_title(
                "θ_pred_FP / θ_GT_total × 100\n(fair: expression-weighted FP)", fontsize=8
            )

    # Legends: methods + JK_MS pairing
    method_handles = [
        plt.Rectangle((0, 0), 1, 1, color=COLORS[m], label=m) for m in METHODS
    ]
    fig.legend(handles=method_handles, loc="upper right", fontsize=8,
               title="Method", title_fontsize=8, framealpha=0.9)

    pair_handles = [
        plt.Rectangle((0, 0), 1, 1, color="#ff69b4", alpha=0.35, label="Run1: sim1(LR)+sim_sr_s2(SR)"),
        plt.Rectangle((0, 0), 1, 1, color="#5b9bd5", alpha=0.35, label="Run2: sim2(LR)+sim_sr_s1(SR)"),
    ]
    fig.legend(handles=pair_handles, loc="upper left", fontsize=8,
               framealpha=0.9, title="JK_MS pairing", title_fontsize=8)

    fig.tight_layout()
    return fig


# ============================================================
# Main
# ============================================================

def main() -> None:
    """
    Generate 3 comparison figures and save to figures/ in each experiment folder.

    Run:
        cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
        conda activate NanoCount_5
        python analysis/plot_method_comparison.py
    """
    timestamp = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")

    figures = [
        ("contingency", CONTINGENCY_METRICS,
         "Contingency Table — LK vs JK_SS vs JK_MS",
         f"contingency_{timestamp}"),
        ("active", PERF_METRICS,
         "Active Universe Metrics — LK vs JK_SS vs JK_MS",
         f"active_{timestamp}"),
        ("gt_nz", PERF_METRICS,
         "GT Non-zero Metrics — LK vs JK_SS vs JK_MS",
         f"gt_nonzero_{timestamp}"),
    ]

    figs = []
    for set_key, metrics, title, fname in figures:
        print(f"[plot_method_comparison] Building: {fname}")
        fig = make_figure(set_key, metrics, title, fname, timestamp)
        figs.append((fname, fig))

    # Save to figures/ in every experiment folder
    for folder in ALL_FOLDERS:
        fig_dir = Path(RESULTS_BASE) / folder / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        for fname, fig in figs:
            out = fig_dir / f"{fname}.png"
            fig.savefig(out, dpi=150, bbox_inches="tight")
            print(f"  Saved: {out}")

    for _, fig in figs:
        plt.close(fig)

    # --- FP% figure ---
    fp_fig = make_fp_figure(timestamp)
    for folder in ALL_FOLDERS:
        fig_dir = Path(RESULTS_BASE) / folder / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        out = fig_dir / f"fp_rate_{timestamp}.png"
        fp_fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"  Saved: {out}")
    plt.close(fp_fig)

    print(f"[plot_method_comparison] Done.")


if __name__ == "__main__":
    main()
