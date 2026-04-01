"""
run_gt_gt_comparison.py
=======================
Compute and print the inter-sample correlation between two ground truth files.

Answers the question: how similar are the two GT distributions to each other?
This is the upper-bound baseline — the best any predictor could theoretically achieve
if the two samples were identical.

Computes:
  - Spearman correlation (GT sim1 vs GT sim2)
  - Pearson correlation  (GT sim1 vs GT sim2)
  - On transcripts present in either GT (missing filled with 0)
  - On TP-only subset (transcripts with GT > 0 in both samples)

Output is printed to stdout.

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/JOLI_Kallisto
    conda activate NanoCount_5
    python analysis/run_gt_gt_comparison.py
"""

import os
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# ============================================================
# CONFIG — edit before running; do not edit below
# ============================================================
GT_PATHS = {
    "sim1": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample1_gt.tsv",
    "sim2": "/gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/PB_sample2_gt.tsv",
}
# ============================================================
# END CONFIG
# ============================================================


def load_gt(path: str) -> pd.DataFrame:
    """
    Load a ground truth TSV/CSV file.

    Args:
        path : str -- path to GT file.

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


def corr_str(a, b, label: str) -> None:
    """
    Compute and print Spearman + Pearson for two arrays.

    Args:
        a     : array-like -- first values.
        b     : array-like -- second values.
        label : str        -- description printed in the header.
    """
    sp_result = spearmanr(a, b)
    pe_result = pearsonr(a, b)
    sp = sp_result.statistic if hasattr(sp_result, "statistic") else sp_result[0]
    pe = pe_result.statistic if hasattr(pe_result, "statistic") else pe_result[0]
    print(f"  {label}")
    print(f"    Spearman : {sp:.6f}")
    print(f"    Pearson  : {pe:.6f}")


def main() -> None:
    """
    Load both GT files, merge, and report inter-sample correlation.
    """
    names = list(GT_PATHS.keys())
    assert len(names) == 2, "GT_PATHS must contain exactly 2 samples."
    s1, s2 = names

    print("=" * 55)
    print("GT vs GT inter-sample correlation")
    print("=" * 55)
    print(f"  Sample 1 : {s1}  ({GT_PATHS[s1]})")
    print(f"  Sample 2 : {s2}  ({GT_PATHS[s2]})")
    print()

    df1 = load_gt(GT_PATHS[s1])
    df2 = load_gt(GT_PATHS[s2])
    print(f"  {s1}: {len(df1)} transcripts  (nonzero: {(df1['tpm_gt']>0).sum()})")
    print(f"  {s2}: {len(df2)} transcripts  (nonzero: {(df2['tpm_gt']>0).sum()})")
    print()

    # Merge on transcript_id (outer join — missing → 0)
    merged = df1.merge(df2, on="transcript_id", how="outer",
                       suffixes=(f"_{s1}", f"_{s2}")).fillna(0.0)
    v1 = merged[f"tpm_gt_{s1}"].values
    v2 = merged[f"tpm_gt_{s2}"].values
    print(f"  Union transcript set : {len(merged)} transcripts")

    corr_str(v1, v2, "All transcripts (union, zeros included)")

    # TP-only: nonzero in both
    both_nonzero = (v1 > 0) & (v2 > 0)
    print(f"\n  Transcripts nonzero in both : {both_nonzero.sum()}")
    corr_str(v1[both_nonzero], v2[both_nonzero], "TP-only (GT > 0 in both samples)")

    # Nonzero in either
    either_nonzero = (v1 > 0) | (v2 > 0)
    print(f"\n  Transcripts nonzero in either : {either_nonzero.sum()}")
    corr_str(v1[either_nonzero], v2[either_nonzero], "Nonzero in either sample")

    print()
    print("=" * 55)


if __name__ == "__main__":
    main()
