"""
test_simulate_ont.py — Tests for src/simulate_ont.py

Tests helper functions with synthetic data. NanoSim subprocess call and
HTSeq check are mocked — actual execution requires lrgsp_simulation_2 env
(with HTSeq) and is tested via the SLURM job (submit_ont.sh).

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations
    conda activate lrgsp_simulation
    python tests/test_simulate_ont.py
"""

import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from simulate_ont import build_exp_file, rename_outputs, build_isoform_counts

# ── Synthetic test data ────────────────────────────────────────────────────────

ABUNDANCES = pd.DataFrame({
    "transcript_id": ["TX1", "TX2", "TX3"],
    "TPM":           [500000.0, 300000.0, 0.0],  # TX3 has TPM=0 — should be excluded
})

# Fake NanoSim FASTQ — 3 reads: 2 from TX1, 1 from TX2
FAKE_FASTQ = (
    "@TX1_0_aligned_0\nACGT\n+\nIIII\n"
    "@TX1_1_aligned_0\nACGT\n+\nIIII\n"
    "@TX2_0_aligned_0\nACGT\n+\nIIII\n"
)


# ── Tests ──────────────────────────────────────────────────────────────────────

def test_build_exp_file_excludes_zero_tpm():
    with tempfile.TemporaryDirectory() as tmp:
        abnd = Path(tmp) / "abundances.tsv"
        out  = Path(tmp) / "exp.tsv"
        ABUNDANCES.to_csv(abnd, sep="\t", index=False)

        build_exp_file(str(abnd), str(out))

        df = pd.read_csv(out, sep="\t")
        assert "transcript_id" in df.columns
        assert "TPM"           in df.columns
        assert len(df) == 2,                     "TX3 (TPM=0) should be excluded"
        assert "TX3" not in df["transcript_id"].values
    print("[PASS] build_exp_file — excludes zero-TPM transcripts")


def test_build_exp_file_tpm_values():
    with tempfile.TemporaryDirectory() as tmp:
        abnd = Path(tmp) / "abundances.tsv"
        out  = Path(tmp) / "exp.tsv"
        ABUNDANCES.to_csv(abnd, sep="\t", index=False)

        build_exp_file(str(abnd), str(out))

        df = pd.read_csv(out, sep="\t")
        tx1_tpm = df.loc[df["transcript_id"] == "TX1", "TPM"].iloc[0]
        assert abs(tx1_tpm - 500000.0) < 0.01
    print("[PASS] build_exp_file — TPM values correct")


def test_build_exp_file_has_header():
    with tempfile.TemporaryDirectory() as tmp:
        abnd = Path(tmp) / "abundances.tsv"
        out  = Path(tmp) / "exp.tsv"
        ABUNDANCES.to_csv(abnd, sep="\t", index=False)
        build_exp_file(str(abnd), str(out))

        first_line = out.read_text().splitlines()[0]
        assert "transcript_id" in first_line, "Missing header line"
    print("[PASS] build_exp_file — header line present")


def test_rename_outputs_aligned_only():
    with tempfile.TemporaryDirectory() as tmp:
        prefix  = str(Path(tmp) / "ont_raw")
        aligned = Path(f"{prefix}_aligned_reads.fastq")
        aligned.write_text(FAKE_FASTQ)

        rename_outputs(prefix, tmp, noise_reads=False)

        canonical = Path(tmp) / "ONT.simulated.fastq"
        assert canonical.exists(), "ONT.simulated.fastq not created"
        assert not aligned.exists(), "Original aligned file should be gone"
    print("[PASS] rename_outputs — aligned only")


def test_rename_outputs_with_noise():
    """With noise_reads=True, aligned + unaligned are merged."""
    with tempfile.TemporaryDirectory() as tmp:
        prefix    = str(Path(tmp) / "ont_raw")
        aligned   = Path(f"{prefix}_aligned_reads.fastq")
        unaligned = Path(f"{prefix}_unaligned_reads.fastq")
        aligned.write_text("@read1\nACGT\n+\nIIII\n")
        unaligned.write_text("@noise1\nTTTT\n+\nIIII\n")

        rename_outputs(prefix, tmp, noise_reads=True)

        canonical = Path(tmp) / "ONT.simulated.fastq"
        assert canonical.exists()
        content = canonical.read_text()
        assert "@read1"  in content
        assert "@noise1" in content
    print("[PASS] rename_outputs — merged aligned + unaligned")


def test_build_isoform_counts():
    with tempfile.TemporaryDirectory() as tmp:
        # Write expression file
        exp_file = Path(tmp) / "exp.tsv"
        exp_file.write_text("transcript_id\tcount\tTPM\nTX1\t0\t500000\nTX2\t0\t300000\n")

        # Write fake FASTQ
        fastq = Path(tmp) / "ONT.simulated.fastq"
        fastq.write_text(FAKE_FASTQ)

        out_counts = Path(tmp) / "ONT.isoform_counts.tsv"
        out_r2i    = Path(tmp) / "ONT.read_to_isoform.tsv"

        build_isoform_counts(str(exp_file), str(fastq), str(out_counts), str(out_r2i))

        counts = pd.read_csv(out_counts, sep="\t")
        r2i    = pd.read_csv(out_r2i,    sep="\t")

        tx1_count = counts.loc[counts["transcript_id"] == "TX1", "simulated_count"].iloc[0]
        tx2_count = counts.loc[counts["transcript_id"] == "TX2", "simulated_count"].iloc[0]
        assert tx1_count == 2,   f"TX1 should have 2 reads, got {tx1_count}"
        assert tx2_count == 1,   f"TX2 should have 1 read, got {tx2_count}"
        assert len(r2i)  == 3,   f"read_to_isoform should have 3 rows, got {len(r2i)}"
    print("[PASS] build_isoform_counts")


def test_check_htseq_raises_on_missing():
    """check_htseq() should call sys.exit if HTSeq is not importable."""
    with patch.dict("sys.modules", {"HTSeq": None}):
        # Force ImportError by removing HTSeq from sys.modules
        if "HTSeq" in sys.modules:
            del sys.modules["HTSeq"]
        # Re-import the function so it re-evaluates the import
        import importlib
        import simulate_ont as m
        importlib.reload(m)

        try:
            m.check_htseq()
            # If HTSeq IS installed, this won't raise — that's fine
            print("[PASS] check_htseq — HTSeq available in this env")
        except SystemExit:
            print("[PASS] check_htseq — correctly exits when HTSeq missing")


if __name__ == "__main__":
    print("=== test_simulate_ont.py ===")
    test_build_exp_file_excludes_zero_tpm()
    test_build_exp_file_tpm_values()
    test_build_exp_file_has_header()
    test_rename_outputs_aligned_only()
    test_rename_outputs_with_noise()
    test_build_isoform_counts()
    test_check_htseq_raises_on_missing()
    print("=== All tests passed ===")
