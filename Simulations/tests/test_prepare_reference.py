"""
test_prepare_reference.py — Tests for src/prepare_reference.py

Uses a tiny synthetic genome + GTF in a temp directory — no real data needed.

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations
    conda activate lrgsp_simulation
    python tests/test_prepare_reference.py
"""

import sys
import tempfile
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from prepare_reference import (
    parse_gtf,
    extract_transcripts,
    write_fasta,
    write_filtered_gtf,
    _revcomp,
)

# ── Synthetic test data ───────────────────────────────────────────────────────
#
# Genome: chr1 = "AAACCCGGGTTT" (12 bp)
#
# Transcripts:
#   TX1  (+) exon 1-3  + exon 7-9  → "AAA" + "GGG" = "AAAGGG" + polyA
#   TX2  (-) exon 4-6  + exon 10-12 → genomic "CCG" + "TTT" revcomp = revcomp("CCCGGG") wait
#
# Let me be explicit:
#   chr1 = AAACCCGGGTTT
#          123456789012
#
#   TX1 (+): exons (1,3) and (7,9)  → seq = "AAA" + "GGG" = "AAAGGG"
#   TX2 (-): exons (4,6) and (10,12) → genomic "CCC" + "TTT"
#            concat in genomic order → "CCCTTT"
#            revcomp → "AAAGGG"
#
# Both transcripts should produce "AAAGGG" + polyA

GENOME_SEQ = "AAACCCGGGTTT"  # 12 bp

GTF_LINES = """\
chr1\t.\texon\t1\t3\t.\t+\t.\tgene_id "G1"; transcript_id "TX1";
chr1\t.\texon\t7\t9\t.\t+\t.\tgene_id "G1"; transcript_id "TX1";
chr1\t.\texon\t4\t6\t.\t-\t.\tgene_id "G2"; transcript_id "TX2";
chr1\t.\texon\t10\t12\t.\t-\t.\tgene_id "G2"; transcript_id "TX2";
"""


def _write_genome_fasta(path: Path) -> None:
    """Write the synthetic genome FASTA and index it with pysam."""
    path.write_text(f">chr1\n{GENOME_SEQ}\n")
    pysam.faidx(str(path))


def _write_gtf(path: Path) -> None:
    path.write_text(GTF_LINES)


# ── Tests ─────────────────────────────────────────────────────────────────────

def test_revcomp():
    assert _revcomp("AAAGGG") == "CCCTTT"
    assert _revcomp("CCCTTT") == "AAAGGG"
    assert _revcomp("ACGT")   == "ACGT"
    print("[PASS] _revcomp")


def test_parse_gtf():
    with tempfile.TemporaryDirectory() as tmp:
        gtf = Path(tmp) / "test.gtf"
        _write_gtf(gtf)
        df = parse_gtf(str(gtf))
        assert len(df) == 4,                   f"Expected 4 exon rows, got {len(df)}"
        assert df["transcript_id"].nunique() == 2, "Expected 2 transcripts"
        assert set(df["transcript_id"]) == {"TX1", "TX2"}
    print("[PASS] parse_gtf")


def test_extract_transcripts_plus_strand():
    """TX1 on + strand: exons (1,3) + (7,9) → AAAGGG + polyA."""
    with tempfile.TemporaryDirectory() as tmp:
        fa  = Path(tmp) / "genome.fa"
        gtf = Path(tmp) / "test.gtf"
        _write_genome_fasta(fa)
        _write_gtf(gtf)
        df   = parse_gtf(str(gtf))
        seqs = extract_transcripts(df, str(fa), polya_len=5)

        expected = "AAAGGG" + "A" * 5
        assert seqs["TX1"] == expected, f"TX1 mismatch: {seqs['TX1']!r} != {expected!r}"
    print("[PASS] extract_transcripts (+ strand)")


def test_extract_transcripts_minus_strand():
    """TX2 on - strand: genomic exons (4,6)=CCC + (10,12)=TTT → revcomp(CCCTTT) = AAAGGG + polyA."""
    with tempfile.TemporaryDirectory() as tmp:
        fa  = Path(tmp) / "genome.fa"
        gtf = Path(tmp) / "test.gtf"
        _write_genome_fasta(fa)
        _write_gtf(gtf)
        df   = parse_gtf(str(gtf))
        seqs = extract_transcripts(df, str(fa), polya_len=5)

        expected = "AAAGGG" + "A" * 5
        assert seqs["TX2"] == expected, f"TX2 mismatch: {seqs['TX2']!r} != {expected!r}"
    print("[PASS] extract_transcripts (- strand)")


def test_write_fasta():
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "out.fasta"
        write_fasta({"TX1": "AAAGGG", "TX2": "CCCTTT"}, str(out))
        content = out.read_text()
        assert ">TX1\n" in content
        assert ">TX2\n" in content
        assert "AAAGGG" in content
    print("[PASS] write_fasta")


def test_write_filtered_gtf():
    with tempfile.TemporaryDirectory() as tmp:
        gtf     = Path(tmp) / "test.gtf"
        out_gtf = Path(tmp) / "out.gtf"
        _write_gtf(gtf)
        write_filtered_gtf(str(gtf), {"TX1"}, str(out_gtf))
        lines = [l for l in out_gtf.read_text().splitlines() if "TX2" in l]
        assert len(lines) == 0, "TX2 should have been filtered out"
        lines_tx1 = [l for l in out_gtf.read_text().splitlines() if "TX1" in l]
        assert len(lines_tx1) == 2, "TX1 should have 2 exon lines"
    print("[PASS] write_filtered_gtf")


if __name__ == "__main__":
    print("=== test_prepare_reference.py ===")
    test_revcomp()
    test_parse_gtf()
    test_extract_transcripts_plus_strand()
    test_extract_transcripts_minus_strand()
    test_write_fasta()
    test_write_filtered_gtf()
    print("=== All tests passed ===")
