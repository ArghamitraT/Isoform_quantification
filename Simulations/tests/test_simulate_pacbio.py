"""
test_simulate_pacbio.py — Tests for src/simulate_pacbio.py

Tests helper functions with synthetic data. IsoSeqSim subprocess call is
mocked — actual execution is tested via the SLURM job (submit_pacbio.sh).

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations
    conda activate lrgsp_simulation
    python tests/test_simulate_pacbio.py
"""

import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from simulate_pacbio import gtf_to_gpd, build_expr_file, parse_transcript_output

# ── Synthetic test data ────────────────────────────────────────────────────────
#
# chr1 = 20 bp
# TX1 (+): exons 1-5, 11-15  → 2 exons
# TX2 (-): exons 6-10, 16-20 → 2 exons

GTF_CONTENT = """\
chr1\t.\texon\t1\t5\t.\t+\t.\tgene_id "G1"; transcript_id "TX1";
chr1\t.\texon\t11\t15\t.\t+\t.\tgene_id "G1"; transcript_id "TX1";
chr1\t.\texon\t6\t10\t.\t-\t.\tgene_id "G2"; transcript_id "TX2";
chr1\t.\texon\t16\t20\t.\t-\t.\tgene_id "G2"; transcript_id "TX2";
"""

ABUNDANCES = pd.DataFrame({
    "transcript_id": ["TX1", "TX2"],
    "TPM": [600000.0, 400000.0],
})


def test_gtf_to_gpd_keys():
    with tempfile.TemporaryDirectory() as tmp:
        gtf = Path(tmp) / "test.gtf"
        gtf.write_text(GTF_CONTENT)
        gpd = gtf_to_gpd(str(gtf))
        assert set(gpd.keys()) == {"TX1", "TX2"}
    print("[PASS] gtf_to_gpd — correct transcript IDs")


def test_gtf_to_gpd_coordinates():
    """TX1 (+): exons (0,5),(10,15) → txStart=0, txEnd=15, exonCount=2."""
    with tempfile.TemporaryDirectory() as tmp:
        gtf = Path(tmp) / "test.gtf"
        gtf.write_text(GTF_CONTENT)
        gpd = gtf_to_gpd(str(gtf))
        parts = gpd["TX1"].split("\t")
        assert parts[0] == "TX1"          # name
        assert parts[1] == "chr1"         # chrom
        assert parts[2] == "+"            # strand
        assert parts[3] == "0"            # txStart (0-based: 1-1=0)
        assert parts[4] == "15"           # txEnd   (0-based end: 15)
        assert parts[7] == "2"            # exonCount
    print("[PASS] gtf_to_gpd — coordinates correct")


def test_gtf_to_gpd_minus_strand():
    """TX2 (-): sorted exons (5,10),(15,20) → txStart=5, txEnd=20."""
    with tempfile.TemporaryDirectory() as tmp:
        gtf = Path(tmp) / "test.gtf"
        gtf.write_text(GTF_CONTENT)
        gpd = gtf_to_gpd(str(gtf))
        parts = gpd["TX2"].split("\t")
        assert parts[2] == "-"
        assert parts[3] == "5"   # 6-1=5
        assert parts[4] == "20"
    print("[PASS] gtf_to_gpd — minus strand correct")


def test_build_expr_file_read_counts():
    """TPM 60/40 split → ~60%/40% of total reads."""
    with tempfile.TemporaryDirectory() as tmp:
        gtf  = Path(tmp) / "test.gtf"
        abnd = Path(tmp) / "abundances.tsv"
        out  = Path(tmp) / "expr.txt"

        gtf.write_text(GTF_CONTENT)
        ABUNDANCES.to_csv(abnd, sep="\t", index=False)

        gpd = gtf_to_gpd(str(gtf))
        build_expr_file(str(abnd), gpd, total_reads=1000, out_path=str(out))

        assert out.exists()
        lines = [l for l in out.read_text().splitlines() if l.strip()]
        assert len(lines) == 2

        counts = {l.split("\t")[0]: int(l.split("\t")[-1]) for l in lines}
        assert counts["TX1"] == 600   # 60% of 1000
        assert counts["TX2"] == 400   # 40% of 1000
    print("[PASS] build_expr_file — read counts correct")


def test_build_expr_file_skips_missing():
    """Transcripts not in GPD should be silently skipped."""
    with tempfile.TemporaryDirectory() as tmp:
        abnd = Path(tmp) / "abundances.tsv"
        out  = Path(tmp) / "expr.txt"

        extra = pd.DataFrame({
            "transcript_id": ["TX1", "TX2", "TX_MISSING"],
            "TPM": [500000.0, 300000.0, 200000.0],
        })
        extra.to_csv(abnd, sep="\t", index=False)

        gtf = Path(tmp) / "test.gtf"
        gtf.write_text(GTF_CONTENT)
        gpd = gtf_to_gpd(str(gtf))

        build_expr_file(str(abnd), gpd, total_reads=1000, out_path=str(out))

        lines = [l for l in out.read_text().splitlines() if l.strip()]
        tx_ids = [l.split("\t")[0] for l in lines]
        assert "TX_MISSING" not in tx_ids
    print("[PASS] build_expr_file — missing transcripts skipped")


def test_parse_transcript_output():
    """Parse IsoSeqSim -t output → isoform_counts + read_to_isoform."""
    with tempfile.TemporaryDirectory() as tmp:
        gtf = Path(tmp) / "test.gtf"
        gtf.write_text(GTF_CONTENT)
        gpd = gtf_to_gpd(str(gtf))

        # Simulate IsoSeqSim -t output: GPD line + read count
        tx_file = Path(tmp) / "transcript.txt"
        tx_file.write_text(
            gpd["TX1"] + "\t600\n" +
            gpd["TX2"] + "\t0\n"    # zero count
        )

        counts_out = Path(tmp) / "PacBio.isoform_counts.tsv"
        r2i_out    = Path(tmp) / "PacBio.read_to_isoform.tsv"
        parse_transcript_output(str(tx_file), str(counts_out), str(r2i_out))

        counts = pd.read_csv(counts_out, sep="\t")
        r2i    = pd.read_csv(r2i_out,    sep="\t")

        assert len(counts) == 2                        # all transcripts in counts
        assert len(r2i)    == 1                        # only TX1 has reads
        assert r2i.iloc[0]["transcript_id"] == "TX1"
        assert r2i.iloc[0]["simulated_count"] == 600
    print("[PASS] parse_transcript_output")


if __name__ == "__main__":
    print("=== test_simulate_pacbio.py ===")
    test_gtf_to_gpd_keys()
    test_gtf_to_gpd_coordinates()
    test_gtf_to_gpd_minus_strand()
    test_build_expr_file_read_counts()
    test_build_expr_file_skips_missing()
    test_parse_transcript_output()
    print("=== All tests passed ===")
