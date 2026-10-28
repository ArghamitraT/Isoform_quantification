"""
plot_method_comparison_realdataReplicate.py
=============================================
Bar-chart comparison of LK, JK_SS, and JK_MS across the real-data replicate
pairs (analogous to plot_method_comparison.py, but for real data +
replicate comparisons instead of simulation + ground-truth comparisons).

Parses replicate_comparison_*.txt reports (written by
analysis/run_replicate_comparison.py) for all 3 methods x 4 replicate
groups and produces 4 figures:
  1. all_transcripts_{timestamp}.png     — Spearman/Pearson/MAE/RMSE/Bray-Curtis,
                                            "all transcripts" universe
  2. active_universe_{timestamp}.png     — same 5 metrics, "active universe"
  3. concordant_positive_{timestamp}.png — same 5 metrics, "concordant positive" universe
  4. concordance_counts_{timestamp}.png  — Concordant Positive / Rep1-only /
                                            Rep2-only / Concordant Negative counts
                                            (no precision/recall/F1)

Layout of each figure (same structure as plot_method_comparison.py):
  rows : replicate groups (flnc_31v32 "LR", flnc_51v52 "LR", ill_31v32 "SR", ill_51v52 "SR")
  cols : metrics (varies per figure)
  bars : methods (LK=light grey, JK_SS=mid grey, JK_MS=dark grey)

Outputs saved to figures/ inside every experiment folder involved (LK, JK_SS,
and JK_MS folders alike), following analysis output convention Rule 3 and
matching plot_method_comparison.py's precedent of copying into every
contributing folder.

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/plot_method_comparison_realdataReplicate.py
"""

import re
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# CONFIG — edit report paths / folders if you rerun the comparisons
# ============================================================

RESULTS_BASE = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"

METHODS = ["LK", "JK_SS", "JK_MS"]
COLORS  = {"LK": "#cccccc", "JK_SS": "#888888", "JK_MS": "#333333"}   # light -> dark grey

GROUPS = ["flnc_31 vs 32\n(LR)", "flnc_51 vs 52\n(LR)", "ill_31 vs 32\n(SR)", "ill_51 vs 52\n(SR)"]

# Groups whose JK_MS values were produced together as one MS run pair get the same
# highlight color on their y-axis label:
#   Type A = Run1(flnc_31+ill_51) + Run2(flnc_32+ill_52)  -> pink
#   Type B = Run3(flnc_51+ill_31) + Run4(flnc_52+ill_32)  -> blue
PAIR_COLORS = {
    "flnc_31 vs 32\n(LR)": "#ff69b4",  # pink — Type A
    "ill_51 vs 52\n(SR)":  "#ff69b4",  # pink — Type A
    "flnc_51 vs 52\n(LR)": "#5b9bd5",  # blue — Type B
    "ill_31 vs 32\n(SR)":  "#5b9bd5",  # blue — Type B
}


def _highlight_ylabel(ax, group: str, fontsize: int = 8, labelpad: int = 6) -> None:
    """
    Set the row's y-axis label with a colored highlight box matching the
    JK_MS run pair that produced it (see PAIR_COLORS).

    Args:
        ax       : matplotlib Axes
        group    : str  -- group key, e.g. "flnc_31 vs 32\\n(LR)"
        fontsize : int
        labelpad : int

    Returns:
        None
    """
    color = PAIR_COLORS.get(group)
    kwargs = dict(fontsize=fontsize, fontweight="bold", labelpad=labelpad)
    if color:
        kwargs["bbox"] = dict(facecolor=color, alpha=0.35, edgecolor=color,
                               boxstyle="round,pad=0.35")
    ax.set_ylabel(group, **kwargs)

# report path for each (group, method) cell
REPORTS = {
    "flnc_31 vs 32\n(LR)": {
        "LK":    f"{RESULTS_BASE}/exprmnt_2026_07_19__01_26_46/replicate_comparison_2026_07_19__01_32_07.txt",
        "JK_SS": f"{RESULTS_BASE}/exprmnt_2026_07_19__01_34_25/replicate_comparison_2026_07_19__01_37_09.txt",
        "JK_MS": f"{RESULTS_BASE}/exprmnt_2026_07_19__17_34_14/replicate_comparison_2026_07_19__17_48_23.txt",
    },
    "flnc_51 vs 52\n(LR)": {
        "LK":    f"{RESULTS_BASE}/exprmnt_2026_07_19__01_26_46/replicate_comparison_2026_07_19__01_32_08.txt",
        "JK_SS": f"{RESULTS_BASE}/exprmnt_2026_07_19__01_34_25/replicate_comparison_2026_07_19__01_37_10.txt",
        "JK_MS": f"{RESULTS_BASE}/exprmnt_2026_07_19__17_34_16/replicate_comparison_2026_07_19__17_48_50.txt",
    },
    "ill_31 vs 32\n(SR)": {
        "LK":    f"{RESULTS_BASE}/exprmnt_2026_07_19__02_29_08/replicate_comparison_2026_07_19__03_15_05.txt",
        "JK_SS": f"{RESULTS_BASE}/exprmnt_2026_07_19__01_44_26/replicate_comparison_2026_07_19__02_28_51.txt",
        "JK_MS": f"{RESULTS_BASE}/exprmnt_2026_07_19__17_34_16/replicate_comparison_2026_07_19__17_49_16.txt",
    },
    "ill_51 vs 52\n(SR)": {
        "LK":    f"{RESULTS_BASE}/exprmnt_2026_07_19__02_29_08/replicate_comparison_2026_07_19__03_15_07.txt",
        "JK_SS": f"{RESULTS_BASE}/exprmnt_2026_07_19__01_44_26/replicate_comparison_2026_07_19__02_28_52.txt",
        "JK_MS": f"{RESULTS_BASE}/exprmnt_2026_07_19__17_34_14/replicate_comparison_2026_07_19__17_48_36.txt",
    },
}

# Every folder that contributed a report -- each gets a copy of all figures.
ALL_FOLDERS = [
    "exprmnt_2026_07_19__01_26_46",   # LK  real PacBio LR
    "exprmnt_2026_07_19__02_29_08",   # LK  real Illumina SR
    "exprmnt_2026_07_19__01_34_25",   # JK_SS real PacBio LR
    "exprmnt_2026_07_19__01_44_26",   # JK_SS real Illumina SR
    "exprmnt_2026_07_19__17_34_14",   # JK_MS Run1: flnc_31 + ill_51
    "exprmnt_2026_07_19__17_34_15",   # JK_MS Run2: flnc_32 + ill_52
    "exprmnt_2026_07_19__17_34_16",   # JK_MS Run3: flnc_51 + ill_31
    "exprmnt_2026_07_19__17_34_17",   # JK_MS Run4: flnc_52 + ill_32
]

# ============================================================

PERF_METRICS = ["Spearman", "Pearson", "MAE", "RMSE", "Bray-Curtis"]
CONCORDANCE_METRICS = ["Concordant Positive", "Rep1-only", "Rep2-only", "Concordant Negative"]

# Metrics/categories where lower = better (winner gets a red-edge highlight)
LOWER_BETTER = {"MAE", "RMSE", "Bray-Curtis", "Rep1-only", "Rep2-only", "Concordant Negative"}


# ============================================================
# Report parsing
# ============================================================

_METRIC_LINE_RE = re.compile(r"(Spearman|Pearson|MAE|RMSE|Bray-Curtis)\s*:\s*([-\d.]+)")
_UNIVERSE_HEADER_RE = re.compile(r"--- Metrics \(([^)]+)\) ---")
_CONCORDANCE_RES = {
    "Concordant Positive": re.compile(r"Concordant positive\s+\(both > 0\)\s*:\s*([\d,]+)"),
    "Rep1-only":           re.compile(r"Rep1-only[^\n]*:\s*([\d,]+)"),
    "Rep2-only":           re.compile(r"Rep2-only[^\n]*:\s*([\d,]+)"),
    "Concordant Negative": re.compile(r"Concordant negative\s+\(both = 0\)\s*:\s*([\d,]+)"),
}


def _classify_universe(header_text: str) -> str:
    """
    Description:
        Map a "--- Metrics (...) ---" header string to a short universe key.
    Args:
        header_text (str): the text captured inside the header parentheses.
    Returns:
        str: one of "all", "active", "concordant".
    """
    if header_text.startswith("all transcripts"):
        return "all"
    if header_text.startswith("active universe"):
        return "active"
    if header_text.startswith("concordant positive only"):
        return "concordant"
    raise ValueError(f"Unrecognized metrics universe header: {header_text!r}")


def parse_report(path: str) -> dict:
    """
    Description:
        Parse a replicate_comparison_*.txt report (written by
        run_replicate_comparison.py) into a metrics dict.
    Args:
        path (str): path to the report .txt file.
    Returns:
        dict: {
            "concordance": {category -> int count},
            "all":         {metric -> float},
            "active":      {metric -> float},
            "concordant":  {metric -> float},
        }
    """
    text = Path(path).read_text()

    concordance = {}
    for key, pattern in _CONCORDANCE_RES.items():
        m = pattern.search(text)
        if not m:
            raise ValueError(f"Could not find '{key}' in {path}")
        concordance[key] = int(m.group(1).replace(",", ""))

    parts = _UNIVERSE_HEADER_RE.split(text)
    universes = {}
    for i in range(1, len(parts), 2):
        header, block = parts[i], parts[i + 1]
        uni_key = _classify_universe(header)
        universes[uni_key] = {m: float(v) for m, v in _METRIC_LINE_RE.findall(block)}

    for uni_key in ("all", "active", "concordant"):
        if uni_key not in universes:
            raise ValueError(f"Missing '{uni_key}' universe metrics in {path}")

    return {"concordance": concordance, **universes}


def build_data() -> dict:
    """
    Description:
        Parse every (group, method) report and assemble the DATA structure
        used for plotting: DATA[group][method][universe_or_concordance][metric].
    Returns:
        dict: nested DATA structure.
    """
    data = {}
    for group in GROUPS:
        data[group] = {}
        for method in METHODS:
            data[group][method] = parse_report(REPORTS[group][method])
    return data


# ============================================================
# Plotting helpers (same structure as plot_method_comparison.py)
# ============================================================

def _bar_group(ax, values: dict, title: str, ylabel: str = "") -> None:
    """
    Draw a grouped bar chart with one group per method.

    Args:
        ax     : matplotlib Axes
        values : dict method -> float
        title  : str -- subplot title
        ylabel : str -- y-axis label
    """
    x     = np.arange(len(METHODS))
    width = 0.55
    bars  = [values[m] for m in METHODS]
    cols  = [COLORS[m] for m in METHODS]

    rects = ax.bar(x, bars, width, color=cols, edgecolor="white", linewidth=0.5)

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


def make_figure(data: dict, set_key: str, metric_keys: list, title: str) -> plt.Figure:
    """
    Build a figure: rows=groups, cols=metrics, bars=methods.

    Args:
        data        : dict  -- DATA structure from build_data()
        set_key     : str   -- "all" | "active" | "concordant" | "concordance"
        metric_keys : list  -- metric names to plot (columns)
        title       : str   -- figure suptitle

    Returns:
        plt.Figure
    """
    n_rows = len(GROUPS)
    n_cols = len(metric_keys)
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(2.8 * n_cols, 3.0 * n_rows),
        squeeze=False,
    )
    fig.suptitle(title, fontsize=12, fontweight="bold", y=1.01)

    for r, group in enumerate(GROUPS):
        for c, metric in enumerate(metric_keys):
            ax = axes[r][c]
            values = {m: data[group][m][set_key][metric] for m in METHODS}
            _bar_group(ax, values, metric)
            if c == 0:
                _highlight_ylabel(ax, group)

    method_handles = [plt.Rectangle((0, 0), 1, 1, color=COLORS[m], label=m) for m in METHODS]
    fig.legend(handles=method_handles, loc="upper right", fontsize=9,
               framealpha=0.9, title="Method", title_fontsize=9)

    pair_handles = [
        plt.Rectangle((0, 0), 1, 1, color="#ff69b4", alpha=0.35,
                       label="Type A: flnc_31+ill_51 / flnc_32+ill_52"),
        plt.Rectangle((0, 0), 1, 1, color="#5b9bd5", alpha=0.35,
                       label="Type B: flnc_51+ill_31 / flnc_52+ill_32"),
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
    Parse all 12 replicate comparison reports (4 groups x 3 methods) and
    generate 4 comparison figures, saved to figures/ in every contributing
    experiment folder.

    Run:
        cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
        conda activate NanoCount_5
        python analysis/plot_method_comparison_realdataReplicate.py
    """
    timestamp = datetime.now().strftime("%Y_%m_%d__%H_%M_%S")

    print("[plot_method_comparison_realdataReplicate] Parsing reports...")
    data = build_data()

    figures = [
        ("all", PERF_METRICS,
         "All Transcripts Metrics — LK vs JK_SS vs JK_MS (Real Data Replicates)",
         f"all_transcripts_{timestamp}"),
        ("active", PERF_METRICS,
         "Active Universe Metrics — LK vs JK_SS vs JK_MS (Real Data Replicates)",
         f"active_universe_{timestamp}"),
        ("concordant", PERF_METRICS,
         "Concordant Positive Metrics — LK vs JK_SS vs JK_MS (Real Data Replicates)",
         f"concordant_positive_{timestamp}"),
        ("concordance", CONCORDANCE_METRICS,
         "Concordance Counts — LK vs JK_SS vs JK_MS (Real Data Replicates)",
         f"concordance_counts_{timestamp}"),
    ]

    figs = []
    for set_key, metrics, title, fname in figures:
        print(f"  Building: {fname}")
        fig = make_figure(data, set_key, metrics, title)
        figs.append((fname, fig))

    for folder in ALL_FOLDERS:
        fig_dir = Path(RESULTS_BASE) / folder / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        for fname, fig in figs:
            out = fig_dir / f"{fname}.png"
            fig.savefig(out, dpi=150, bbox_inches="tight")
            print(f"  Saved: {out}")

    for _, fig in figs:
        plt.close(fig)

    print("[plot_method_comparison_realdataReplicate] Done.")


if __name__ == "__main__":
    main()
