"""
test_generate_abundances.py — Tests for src/generate_abundances.py

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations
    conda activate lrgsp_simulation
    python tests/test_generate_abundances.py
"""

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from generate_abundances import (
    read_transcript_ids,
    mode_lognormal,
    mode_custom,
    mode_uniform,
)

TX_IDS = [f"TX{i}" for i in range(100)]

FASTA_CONTENT = "\n".join(f">{tid}\nACGT" for tid in TX_IDS) + "\n"


def test_read_transcript_ids():
    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / "tx.fasta"
        fa.write_text(FASTA_CONTENT)
        ids = read_transcript_ids(str(fa))
        assert ids == TX_IDS, "IDs do not match"
        assert len(ids) == 100
    print("[PASS] read_transcript_ids")


def test_lognormal_sum_to_1m():
    df = mode_lognormal(TX_IDS, seed=42, dropout=0.0)
    assert abs(df["TPM"].sum() - 1e6) < 1.0, f"TPM sum off: {df['TPM'].sum()}"
    print("[PASS] mode_lognormal — sums to 1e6")


def test_lognormal_dropout():
    dropout = 0.3
    df = mode_lognormal(TX_IDS, seed=42, dropout=dropout)
    zero_frac = (df["TPM"] == 0).mean()
    # Allow ±5% tolerance
    assert abs(zero_frac - dropout) < 0.05, f"Dropout fraction off: {zero_frac:.2f}"
    print(f"[PASS] mode_lognormal — dropout ({zero_frac:.2f} ≈ {dropout})")


def test_lognormal_reproducible():
    df1 = mode_lognormal(TX_IDS, seed=7)
    df2 = mode_lognormal(TX_IDS, seed=7)
    assert df1["TPM"].equals(df2["TPM"]), "Same seed produced different results"
    print("[PASS] mode_lognormal — reproducible")


def test_uniform_equal():
    df = mode_uniform(TX_IDS)
    assert abs(df["TPM"].sum() - 1e6) < 1.0
    assert df["TPM"].nunique() == 1, "Uniform mode should give equal TPM"
    print("[PASS] mode_uniform")


def test_custom_basic():
    with tempfile.TemporaryDirectory() as tmp:
        custom_tsv = Path(tmp) / "custom.tsv"
        # Half the transcripts have TPM=1, half have TPM=0
        rows = [{"transcript_id": tid, "TPM": 1.0 if i < 50 else 0.0}
                for i, tid in enumerate(TX_IDS)]
        pd.DataFrame(rows).to_csv(custom_tsv, sep="\t", index=False)

        df = mode_custom(TX_IDS, str(custom_tsv))
        assert abs(df["TPM"].sum() - 1e6) < 1.0
        assert (df["TPM"] > 0).sum() == 50
    print("[PASS] mode_custom — basic")


def test_custom_cancer_upregulation():
    with tempfile.TemporaryDirectory() as tmp:
        custom_tsv   = Path(tmp) / "custom.tsv"
        cancer_file  = Path(tmp) / "cancer.txt"

        # All transcripts get TPM=1
        rows = [{"transcript_id": tid, "TPM": 1.0} for tid in TX_IDS]
        pd.DataFrame(rows).to_csv(custom_tsv, sep="\t", index=False)

        # Mark first 10 as cancer genes
        cancer_ids = TX_IDS[:10]
        cancer_file.write_text("\n".join(cancer_ids))

        df_base   = mode_custom(TX_IDS, str(custom_tsv))
        df_cancer = mode_custom(TX_IDS, str(custom_tsv),
                                cancer_genes_file=str(cancer_file), fold_change=5.0)

        # Cancer transcripts should have higher TPM than non-cancer
        cancer_mean     = df_cancer[df_cancer["transcript_id"].isin(cancer_ids)]["TPM"].mean()
        noncancer_mean  = df_cancer[~df_cancer["transcript_id"].isin(cancer_ids)]["TPM"].mean()
        assert cancer_mean > noncancer_mean, "Cancer genes not up-regulated"
    print("[PASS] mode_custom — cancer up-regulation")


def test_custom_missing_transcripts_get_zero():
    """Transcripts not in the custom table should receive TPM=0."""
    with tempfile.TemporaryDirectory() as tmp:
        custom_tsv = Path(tmp) / "custom.tsv"
        # Only provide TPM for first 50 transcripts
        rows = [{"transcript_id": tid, "TPM": 1.0} for tid in TX_IDS[:50]]
        pd.DataFrame(rows).to_csv(custom_tsv, sep="\t", index=False)

        df = mode_custom(TX_IDS, str(custom_tsv))
        missing_tpm = df[df["transcript_id"].isin(TX_IDS[50:])]["TPM"]
        assert (missing_tpm == 0).all(), "Missing transcripts should have TPM=0"
    print("[PASS] mode_custom — missing transcripts get TPM=0")


if __name__ == "__main__":
    print("=== test_generate_abundances.py ===")
    test_read_transcript_ids()
    test_lognormal_sum_to_1m()
    test_lognormal_dropout()
    test_lognormal_reproducible()
    test_uniform_equal()
    test_custom_basic()
    test_custom_cancer_upregulation()
    test_custom_missing_transcripts_get_zero()
    print("=== All tests passed ===")
