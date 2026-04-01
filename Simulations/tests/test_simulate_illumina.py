"""
test_simulate_illumina.py — Tests for src/simulate_illumina.py

Tests all helper functions. Subprocess calls (rsem-prepare-reference,
rsem-simulate-reads) are mocked — actual RSEM execution is tested via the
SLURM job (submit_illumina.sh).

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations
    conda activate lrgsp_simulation
    python tests/test_simulate_illumina.py
"""

import sys
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from simulate_illumina import (
    abundances_to_rsem_results,
    rename_outputs,
    write_read_to_isoform,
    prepare_rsem_reference,
)

# ── Synthetic data ─────────────────────────────────────────────────────────────

TX_IDS = ["TX1", "TX2", "TX3"]

FASTA_CONTENT = ">TX1\nACGTACGT\n>TX2\nACGT\n>TX3\nACGTACGTACGT\n"

ABUNDANCES = pd.DataFrame({
    "transcript_id": TX_IDS,
    "TPM": [500000.0, 300000.0, 200000.0],
})


# ── Tests ──────────────────────────────────────────────────────────────────────

def test_abundances_to_rsem_results():
    with tempfile.TemporaryDirectory() as tmp:
        fa   = Path(tmp) / "tx.fasta"
        abnd = Path(tmp) / "abundances.tsv"
        out  = Path(tmp) / "out.isoforms.results"

        fa.write_text(FASTA_CONTENT)
        ABUNDANCES.to_csv(abnd, sep="\t", index=False)

        abundances_to_rsem_results(str(abnd), str(fa), str(out))

        assert out.exists(), "isoforms.results not created"
        df = pd.read_csv(out, sep="\t")

        # Check required columns
        for col in ["transcript_id", "gene_id", "length", "effective_length", "TPM"]:
            assert col in df.columns, f"Missing column: {col}"

        # Check lengths extracted correctly
        assert df.loc[df["transcript_id"] == "TX1", "length"].iloc[0] == 8
        assert df.loc[df["transcript_id"] == "TX2", "length"].iloc[0] == 4
        assert df.loc[df["transcript_id"] == "TX3", "length"].iloc[0] == 12

        # Check TPM values preserved
        assert df.loc[df["transcript_id"] == "TX1", "TPM"].iloc[0] == 500000.0

    print("[PASS] abundances_to_rsem_results")


def test_rename_outputs():
    with tempfile.TemporaryDirectory() as tmp:
        prefix = str(Path(tmp) / "Illumina_raw")

        # Create fake RSEM output files
        Path(f"{prefix}_1.fq").write_text("@read1\nACGT\n+\nIIII\n")
        Path(f"{prefix}_2.fq").write_text("@read1\nACGT\n+\nIIII\n")
        Path(f"{prefix}.sim.isoforms.results").write_text(
            "transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n"
            "TX1\tTX1\t8\t8\t100\t500000\t0\t0\n"
        )

        rename_outputs(prefix, tmp)

        assert (Path(tmp) / "Illumina.simulated_1.fq").exists()
        assert (Path(tmp) / "Illumina.simulated_2.fq").exists()
        assert (Path(tmp) / "Illumina.isoform_counts.tsv").exists()

    print("[PASS] rename_outputs")


def test_write_read_to_isoform():
    with tempfile.TemporaryDirectory() as tmp:
        counts_tsv = Path(tmp) / "Illumina.isoform_counts.tsv"
        counts_tsv.write_text(
            "transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n"
            "TX1\tTX1\t8\t8\t100\t500000\t0\t0\n"
            "TX2\tTX2\t4\t4\t60\t300000\t0\t0\n"
            "TX3\tTX3\t12\t12\t0\t0\t0\t0\n"   # zero count — should be excluded
        )
        out = Path(tmp) / "Illumina.read_to_isoform.tsv"
        write_read_to_isoform(str(counts_tsv), str(out))

        df = pd.read_csv(out, sep="\t")
        assert list(df.columns) == ["transcript_id", "simulated_count"]
        assert len(df) == 2,       "TX3 (count=0) should be excluded"
        assert df["simulated_count"].sum() == 160

    print("[PASS] write_read_to_isoform")


def test_prepare_rsem_reference_skips_if_exists():
    """If RSEM reference already exists, rsem-prepare-reference should NOT be called."""
    with tempfile.TemporaryDirectory() as tmp:
        prefix = str(Path(tmp) / "ref")
        # Simulate existing reference
        Path(f"{prefix}.transcripts.fa").write_text(">TX1\nACGT\n")

        with patch("simulate_illumina._run") as mock_run:
            prepare_rsem_reference("tx.fasta", prefix, threads=4)
            mock_run.assert_not_called()

    print("[PASS] prepare_rsem_reference — skips if already exists")


def test_prepare_rsem_reference_calls_rsem():
    """If RSEM reference does not exist, rsem-prepare-reference should be called."""
    with tempfile.TemporaryDirectory() as tmp:
        prefix = str(Path(tmp) / "ref")

        with patch("simulate_illumina._run") as mock_run:
            prepare_rsem_reference("tx.fasta", prefix, threads=4)
            mock_run.assert_called_once()
            cmd = mock_run.call_args[0][0]
            assert "rsem-prepare-reference" in cmd[0]
            assert "tx.fasta" in cmd
            assert prefix in cmd

    print("[PASS] prepare_rsem_reference — calls rsem-prepare-reference")


if __name__ == "__main__":
    print("=== test_simulate_illumina.py ===")
    test_abundances_to_rsem_results()
    test_rename_outputs()
    test_write_read_to_isoform()
    test_prepare_rsem_reference_skips_if_exists()
    test_prepare_rsem_reference_calls_rsem()
    print("=== All tests passed ===")
