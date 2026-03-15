#!/usr/bin/env python3
"""
Compare two transcript abundance TSV files and quantify how similar or dissimilar they are.

This script is designed to handle abundance tables that may have different column names,
for example:

File A:
    transcript_id    tpm

File B:
    target_id    length    eff_length    est_counts    tpm

Main features
-------------
1. Detects transcript ID and abundance columns automatically (or lets you specify them).
2. Aligns both files by transcript ID.
3. Reports overlap statistics.
4. Computes similarity metrics on the shared transcripts:
   - Pearson correlation
   - Spearman correlation
   - Mean absolute difference
   - RMSE
   - Bray-Curtis distance
   - Jensen-Shannon distance (after normalization)
5. Computes metrics on raw abundance and log1p(abundance).
6. Writes merged tables and summary outputs.

Example
-------
python compare_abundance_files.py \
    --file1 abundance_A.tsv \
    --file2 abundance_B.tsv \
    --outdir compare_abundance

If you want to force the columns:
python compare_abundance_files.py \
    --file1 abundance_A.tsv \
    --file2 abundance_B.tsv \
    --id-col1 transcript_id \
    --id-col2 target_id \
    --value-col1 tpm \
    --value-col2 tpm \
    --outdir compare_abundance
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd


DEFAULT_ID_CANDIDATES = [
    "transcript_id",
    "target_id",
    "Name",
    "name",
    "transcript",
    "target",
    "id",
    "ID",
]

DEFAULT_VALUE_CANDIDATES = [
    "tpm",
    "TPM",
    "est_counts",
    "NumReads",
    "count",
    "counts",
    "abundance",
    "expr",
    "expression",
]


# ---------- utility functions ----------

def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def detect_column(columns: Iterable[str], preferred_names: list[str], kind: str) -> str:
    cols = list(columns)
    colset = set(cols)

    for name in preferred_names:
        if name in colset:
            return name

    lowered = {c.lower(): c for c in cols}
    for name in preferred_names:
        if name.lower() in lowered:
            return lowered[name.lower()]

    raise ValueError(
        f"Could not automatically detect {kind} column. Available columns: {cols}"
    )


def read_abundance_table(
    filepath: Path,
    id_col: Optional[str],
    value_col: Optional[str],
) -> Tuple[pd.DataFrame, str, str]:
    df = pd.read_csv(filepath, sep="\t", comment="#")

    if df.empty:
        raise ValueError(f"File is empty: {filepath}")

    chosen_id = id_col or detect_column(df.columns, DEFAULT_ID_CANDIDATES, "ID")
    chosen_value = value_col or detect_column(df.columns, DEFAULT_VALUE_CANDIDATES, "abundance")

    if chosen_id not in df.columns:
        raise ValueError(f"ID column '{chosen_id}' not found in {filepath}")
    if chosen_value not in df.columns:
        raise ValueError(f"Value column '{chosen_value}' not found in {filepath}")

    out = df[[chosen_id, chosen_value]].copy()
    out.columns = ["transcript_id", "value"]
    out["transcript_id"] = out["transcript_id"].astype(str)
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    out = out.dropna(subset=["transcript_id", "value"])

    # If transcript IDs are duplicated, sum them.
    out = out.groupby("transcript_id", as_index=False)["value"].sum()
    return out, chosen_id, chosen_value


def pearson_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    if np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    xr = pd.Series(x).rank(method="average").to_numpy()
    yr = pd.Series(y).rank(method="average").to_numpy()
    return pearson_corr(xr, yr)


def rmse(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.sqrt(np.mean((x - y) ** 2)))


def mae(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean(np.abs(x - y)))


def bray_curtis_distance(x: np.ndarray, y: np.ndarray) -> float:
    denom = np.sum(np.abs(x + y))
    if denom == 0:
        return 0.0
    return float(np.sum(np.abs(x - y)) / denom)


def normalize_nonnegative(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = np.clip(x, a_min=0.0, a_max=None)
    s = x.sum()
    if s == 0:
        return np.full_like(x, 1.0 / len(x)) if len(x) > 0 else x
    return x / s


def kl_divergence(p: np.ndarray, q: np.ndarray, eps: float = 1e-12) -> float:
    p = np.clip(p, eps, None)
    q = np.clip(q, eps, None)
    return float(np.sum(p * np.log2(p / q)))


def jensen_shannon_distance(x: np.ndarray, y: np.ndarray) -> float:
    p = normalize_nonnegative(x)
    q = normalize_nonnegative(y)
    m = 0.5 * (p + q)
    js_div = 0.5 * kl_divergence(p, m) + 0.5 * kl_divergence(q, m)
    return float(np.sqrt(max(js_div, 0.0)))


def summarize_pair(x: np.ndarray, y: np.ndarray, prefix: str) -> dict:
    return {
        f"{prefix}_pearson": pearson_corr(x, y),
        f"{prefix}_spearman": spearman_corr(x, y),
        f"{prefix}_mae": mae(x, y),
        f"{prefix}_rmse": rmse(x, y),
        f"{prefix}_bray_curtis": bray_curtis_distance(x, y),
        f"{prefix}_js_distance": jensen_shannon_distance(x, y),
    }


def format_float(x: float) -> str:
    if pd.isna(x):
        return "nan"
    return f"{x:.6f}"


# ---------- main comparison ----------

def compare_tables(df1: pd.DataFrame, df2: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    merged = df1.merge(df2, on="transcript_id", how="outer", suffixes=("_1", "_2"))
    merged["present_in_file1"] = merged["value_1"].notna()
    merged["present_in_file2"] = merged["value_2"].notna()
    merged["shared"] = merged["present_in_file1"] & merged["present_in_file2"]
    merged["value_1_filled"] = merged["value_1"].fillna(0.0)
    merged["value_2_filled"] = merged["value_2"].fillna(0.0)
    merged["abs_diff"] = np.abs(merged["value_1_filled"] - merged["value_2_filled"])
    merged["log1p_value_1"] = np.log1p(merged["value_1_filled"])
    merged["log1p_value_2"] = np.log1p(merged["value_2_filled"])
    merged["log1p_abs_diff"] = np.abs(merged["log1p_value_1"] - merged["log1p_value_2"])

    shared = merged.loc[merged["shared"]].copy()
    x_shared = shared["value_1"].to_numpy(dtype=float)
    y_shared = shared["value_2"].to_numpy(dtype=float)

    # full vectors with missing treated as zero
    x_all = merged["value_1_filled"].to_numpy(dtype=float)
    y_all = merged["value_2_filled"].to_numpy(dtype=float)

    results = {
        "n_file1": int(df1.shape[0]),
        "n_file2": int(df2.shape[0]),
        "n_union": int(merged.shape[0]),
        "n_shared": int(shared.shape[0]),
        "n_only_file1": int((merged["present_in_file1"] & ~merged["present_in_file2"]).sum()),
        "n_only_file2": int((~merged["present_in_file1"] & merged["present_in_file2"]).sum()),
        "shared_fraction_of_file1": float(shared.shape[0] / df1.shape[0]) if df1.shape[0] else float("nan"),
        "shared_fraction_of_file2": float(shared.shape[0] / df2.shape[0]) if df2.shape[0] else float("nan"),
        "jaccard_ids": float(shared.shape[0] / merged.shape[0]) if merged.shape[0] else float("nan"),
    }

    # Metrics on shared transcripts only
    if len(x_shared) > 0:
        results.update(summarize_pair(x_shared, y_shared, prefix="shared_raw"))
        results.update(summarize_pair(np.log1p(x_shared), np.log1p(y_shared), prefix="shared_log1p"))
    else:
        for key in [
            "shared_raw_pearson",
            "shared_raw_spearman",
            "shared_raw_mae",
            "shared_raw_rmse",
            "shared_raw_bray_curtis",
            "shared_raw_js_distance",
            "shared_log1p_pearson",
            "shared_log1p_spearman",
            "shared_log1p_mae",
            "shared_log1p_rmse",
            "shared_log1p_bray_curtis",
            "shared_log1p_js_distance",
        ]:
            results[key] = float("nan")

    # Metrics on union, filling absent transcripts with 0
    results.update(summarize_pair(x_all, y_all, prefix="union_zero_filled_raw"))
    results.update(summarize_pair(np.log1p(x_all), np.log1p(y_all), prefix="union_zero_filled_log1p"))

    return merged, results


# ---------- output ----------

def write_summary_text(
    outpath: Path,
    results: dict,
    file1: Path,
    file2: Path,
    id_col1: str,
    id_col2: str,
    value_col1: str,
    value_col2: str,
) -> None:
    lines = []
    lines.append("Abundance comparison summary")
    lines.append("=" * 80)
    lines.append(f"File 1: {file1}")
    lines.append(f"  ID column: {id_col1}")
    lines.append(f"  Value column: {value_col1}")
    lines.append(f"File 2: {file2}")
    lines.append(f"  ID column: {id_col2}")
    lines.append(f"  Value column: {value_col2}")
    lines.append("")

    sections = [
        (
            "Overlap statistics",
            [
                "n_file1",
                "n_file2",
                "n_union",
                "n_shared",
                "n_only_file1",
                "n_only_file2",
                "shared_fraction_of_file1",
                "shared_fraction_of_file2",
                "jaccard_ids",
            ],
        ),
        (
            "Shared transcripts only (raw abundance)",
            [
                "shared_raw_pearson",
                "shared_raw_spearman",
                "shared_raw_mae",
                "shared_raw_rmse",
                "shared_raw_bray_curtis",
                "shared_raw_js_distance",
            ],
        ),
        (
            "Shared transcripts only (log1p abundance)",
            [
                "shared_log1p_pearson",
                "shared_log1p_spearman",
                "shared_log1p_mae",
                "shared_log1p_rmse",
                "shared_log1p_bray_curtis",
                "shared_log1p_js_distance",
            ],
        ),
        (
            "Union of transcripts, missing values filled with 0 (raw abundance)",
            [
                "union_zero_filled_raw_pearson",
                "union_zero_filled_raw_spearman",
                "union_zero_filled_raw_mae",
                "union_zero_filled_raw_rmse",
                "union_zero_filled_raw_bray_curtis",
                "union_zero_filled_raw_js_distance",
            ],
        ),
        (
            "Union of transcripts, missing values filled with 0 (log1p abundance)",
            [
                "union_zero_filled_log1p_pearson",
                "union_zero_filled_log1p_spearman",
                "union_zero_filled_log1p_mae",
                "union_zero_filled_log1p_rmse",
                "union_zero_filled_log1p_bray_curtis",
                "union_zero_filled_log1p_js_distance",
            ],
        ),
    ]

    for title, keys in sections:
        lines.append(title)
        lines.append("-" * len(title))
        for key in keys:
            value = results.get(key, float("nan"))
            if isinstance(value, (int, np.integer)):
                lines.append(f"{key}: {value}")
            else:
                lines.append(f"{key}: {format_float(value)}")
        lines.append("")

    lines.append("Interpretation guide")
    lines.append("-------------------")
    lines.append("- Pearson/Spearman closer to 1 means more similar ranking/scale.")
    lines.append("- MAE/RMSE closer to 0 means smaller differences.")
    lines.append("- Bray-Curtis closer to 0 means more similar composition.")
    lines.append("- Jensen-Shannon distance closer to 0 means more similar normalized distributions.")

    outpath.write_text("\n".join(lines))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare two abundance TSV files.")
    parser.add_argument("--file1", required=True, help="Path to first abundance TSV file")
    parser.add_argument("--file2", required=True, help="Path to second abundance TSV file")
    parser.add_argument("--id-col1", default=None, help="Transcript ID column in file1")
    parser.add_argument("--id-col2", default=None, help="Transcript ID column in file2")
    parser.add_argument("--value-col1", default=None, help="Abundance column in file1")
    parser.add_argument("--value-col2", default=None, help="Abundance column in file2")
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of most different transcripts to write out",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for merged table and summary",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    file1 = Path(args.file1)
    file2 = Path(args.file2)
    outdir = Path(args.outdir)
    safe_mkdir(outdir)

    df1, detected_id1, detected_val1 = read_abundance_table(
        file1, args.id_col1, args.value_col1
    )
    df2, detected_id2, detected_val2 = read_abundance_table(
        file2, args.id_col2, args.value_col2
    )

    merged, results = compare_tables(df1, df2)

    merged_sorted = merged.sort_values("abs_diff", ascending=False)
    top_diff = merged_sorted.head(args.top_n)
    top_log_diff = merged.sort_values("log1p_abs_diff", ascending=False).head(args.top_n)

    merged.to_csv(outdir / "merged_abundance_comparison.tsv", sep="\t", index=False)
    top_diff.to_csv(outdir / "top_differences_raw.tsv", sep="\t", index=False)
    top_log_diff.to_csv(outdir / "top_differences_log1p.tsv", sep="\t", index=False)
    pd.DataFrame([results]).to_csv(outdir / "summary_metrics.tsv", sep="\t", index=False)
    (outdir / "summary_metrics.json").write_text(json.dumps(results, indent=2))

    write_summary_text(
        outdir / "summary_report.txt",
        results,
        file1,
        file2,
        detected_id1,
        detected_id2,
        detected_val1,
        detected_val2,
    )

    print(f"Done. Results written to: {outdir}")
    print(f"Shared transcripts: {results['n_shared']}")
    print(f"Jaccard ID overlap: {format_float(results['jaccard_ids'])}")
    print(f"Shared raw Pearson: {format_float(results['shared_raw_pearson'])}")
    print(f"Shared raw Spearman: {format_float(results['shared_raw_spearman'])}")
    print(f"Union zero-filled log1p JS distance: {format_float(results['union_zero_filled_log1p_js_distance'])}")


if __name__ == "__main__":
    main()
