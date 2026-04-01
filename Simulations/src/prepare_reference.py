"""
prepare_reference.py — Phase 1: GTF + genome FASTA → simulation-ready reference files.

What it does:
    1. Parses a GTF annotation with pandas (no gffutils) to build a transcript→exon table.
    2. Extracts transcript sequences from the genome FASTA using pysam.
    3. Appends a poly-A tail to every transcript sequence.
    4. Optionally integrates novel isoforms from a SQANTI classification TSV
       (skipping 'not_in_catalog' and 'novel_in_catalog' categories).
    5. Writes simulation-ready output files.

Inputs:
    --genome    : path to genome FASTA (e.g. GRCh38.primary_assembly.genome.fa)
    --gtf       : path to GTF annotation (e.g. gencode.v46.annotation.gtf)
    --output    : output path prefix (e.g. "files/refs/human_sim")
    --sqanti    : (optional) SQANTI classification TSV for novel isoforms
    --n_novel   : (optional) max novel isoforms to add from SQANTI (default: 0)
    --polya_len : (optional) poly-A tail length in bp (default: 100)

Outputs:
    <output>.transcripts.fasta   — transcript sequences with poly-A tails
    <output>.annotation.gtf      — GTF filtered to successfully extracted transcripts
    <output>.novel_isoforms.tsv  — novel isoforms added (empty header if none)

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing
    conda activate lrgsp_simulation
    python code/Simulations/src/prepare_reference.py \\
        --genome  /path/to/GRCh38.fa \\
        --gtf     /path/to/gencode.vXX.annotation.gtf \\
        --output  files/refs/human_sim
"""

import argparse
import re
from pathlib import Path

import pandas as pd
import pysam


# ── GTF parsing ──────────────────────────────────────────────────────────────

def parse_gtf(gtf_path: str) -> pd.DataFrame:
    """
    Parse a GTF file into a DataFrame of exon records.

    Only 'exon' feature rows are kept. Extracts gene_id and transcript_id
    from the attributes column via regex. GTF coordinates are 1-based inclusive.

    Args:
        gtf_path (str): Path to the GTF annotation file.

    Returns:
        pd.DataFrame: Columns — chrom, start (1-based), end (1-based),
                      strand, gene_id, transcript_id.
    """
    print(f"[prepare_reference] Parsing GTF: {gtf_path}")

    rows = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields
            rows.append({
                "chrom":         chrom,
                "start":         int(start),
                "end":           int(end),
                "strand":        strand,
                "gene_id":       _extract_attr(attrs, "gene_id"),
                "transcript_id": _extract_attr(attrs, "transcript_id"),
            })

    df = pd.DataFrame(rows)
    print(f"[prepare_reference] GTF parsed: {len(df)} exon rows, "
          f"{df['transcript_id'].nunique()} transcripts")
    return df


def _extract_attr(attrs: str, key: str) -> str:
    """
    Extract a quoted attribute value from a GTF attributes string.

    Args:
        attrs (str): Raw GTF attributes field (9th column).
        key   (str): Attribute name to extract (e.g. 'transcript_id').

    Returns:
        str: Attribute value, or empty string if not found.
    """
    m = re.search(rf'{key}\s+"([^"]+)"', attrs)
    return m.group(1) if m else ""


# ── Sequence extraction ──────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _revcomp(seq: str) -> str:
    """
    Return the reverse complement of a DNA string.

    Args:
        seq (str): DNA sequence (ACGT/N, any case).

    Returns:
        str: Reverse-complemented sequence.
    """
    return seq.translate(_COMPLEMENT)[::-1]


def extract_transcripts(
    exon_df: pd.DataFrame,
    genome_path: str,
    polya_len: int = 100,
) -> dict:
    """
    Fetch transcript sequences from the genome FASTA for every transcript
    in exon_df, then append a poly-A tail.

    GTF coordinates are 1-based inclusive; pysam.fetch uses 0-based half-open.
    Minus-strand transcripts: exons concatenated in ascending genomic order,
    then the whole sequence is reverse-complemented.

    Args:
        exon_df     (pd.DataFrame): Output of parse_gtf().
        genome_path (str):          Path to genome FASTA (pysam-indexable).
        polya_len   (int):          Number of 'A' bases to append.

    Returns:
        dict: {transcript_id: sequence_with_polya}
    """
    print(f"[prepare_reference] Opening genome FASTA: {genome_path}")
    print(f"[prepare_reference] Poly-A tail: {polya_len} bp")

    fasta   = pysam.FastaFile(genome_path)
    polya   = "A" * polya_len
    seqs    = {}
    skipped = 0

    for tx_id, group in exon_df.groupby("transcript_id"):
        chrom  = group["chrom"].iloc[0]
        strand = group["strand"].iloc[0]
        exons  = group.sort_values("start")

        try:
            # pysam uses 0-based half-open: (start-1, end)
            seq = "".join(
                fasta.fetch(chrom, row.start - 1, row.end)
                for row in exons.itertuples()
            )
        except (KeyError, ValueError) as e:
            skipped += 1
            print(f"  [WARN] Skipping {tx_id}: {e}")
            continue

        if strand == "-":
            seq = _revcomp(seq)

        seqs[tx_id] = seq + polya

    fasta.close()

    if skipped:
        print(f"[prepare_reference] Skipped {skipped} transcripts (chrom not in FASTA)")
    print(f"[prepare_reference] Extracted {len(seqs)} transcript sequences")
    return seqs


# ── SQANTI novel isoforms ────────────────────────────────────────────────────

_SQANTI_SKIP = {"not_in_catalog", "novel_in_catalog"}


def load_novel_isoforms(sqanti_tsv: str, n_novel: int) -> pd.DataFrame:
    """
    Load novel isoforms from a SQANTI classification TSV, excluding
    'not_in_catalog' and 'novel_in_catalog' structural categories.

    Args:
        sqanti_tsv (str): Path to SQANTI classification TSV.
        n_novel    (int): Max number of novel isoforms to return.

    Returns:
        pd.DataFrame: Filtered novel isoform rows (up to n_novel).
    """
    print(f"[prepare_reference] Loading SQANTI novel isoforms: {sqanti_tsv}")
    df = pd.read_csv(sqanti_tsv, sep="\t")

    required = {"isoform", "structural_category"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"SQANTI TSV missing expected columns: {missing}")

    df = df[~df["structural_category"].isin(_SQANTI_SKIP)].head(n_novel)
    print(f"[prepare_reference] Novel isoforms retained: {len(df)}")
    return df


# ── Writers ──────────────────────────────────────────────────────────────────

def write_fasta(sequences: dict, out_path: str) -> None:
    """
    Write {id: sequence} to a FASTA file with 80-character line wrapping.

    Args:
        sequences (dict): {transcript_id: sequence}
        out_path  (str):  Destination file path.

    Returns:
        None
    """
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for tx_id, seq in sequences.items():
            fh.write(f">{tx_id}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    print(f"[prepare_reference] Wrote FASTA: {out_path}  ({len(sequences)} entries)")


def write_filtered_gtf(gtf_path: str, kept_txs: set, out_path: str) -> None:
    """
    Copy the input GTF to out_path, keeping only exon lines whose
    transcript_id is in kept_txs.

    Args:
        gtf_path  (str):  Path to original GTF.
        kept_txs  (set):  Set of transcript IDs to retain.
        out_path  (str):  Destination GTF path.

    Returns:
        None
    """
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with open(gtf_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fields = line.split("\t")
            if len(fields) >= 9:
                tx_id = _extract_attr(fields[8], "transcript_id")
                if tx_id in kept_txs:
                    fout.write(line)
                    written += 1
    print(f"[prepare_reference] Wrote annotation GTF: {out_path}  ({written} lines)")


# ── Entry point ──────────────────────────────────────────────────────────────

def main():
    # ============================================================
    # CONFIG — edit these variables before running; do not edit below
    # ============================================================
    default_genome    = ""
    default_gtf       = ""
    default_output    = "simulation_ref"
    default_sqanti    = ""
    default_n_novel   = 0
    default_polya_len = 100
    # ============================================================

    parser = argparse.ArgumentParser(
        description="Phase 1: prepare simulation-ready references from GTF + genome FASTA."
    )
    parser.add_argument("--genome",    required=True,                help="Path to genome FASTA")
    parser.add_argument("--gtf",       required=True,                help="Path to GTF annotation")
    parser.add_argument("--output",    default=default_output,       help="Output path prefix")
    parser.add_argument("--sqanti",    default=default_sqanti,       help="(optional) SQANTI TSV")
    parser.add_argument("--n_novel",   default=default_n_novel,      type=int)
    parser.add_argument("--polya_len", default=default_polya_len,    type=int)
    args = parser.parse_args()

    print("=" * 60)
    print("Phase 1 — Reference Preparation")
    print(f"  genome    : {args.genome}")
    print(f"  gtf       : {args.gtf}")
    print(f"  output    : {args.output}")
    print(f"  sqanti    : {args.sqanti or '(none)'}")
    print(f"  n_novel   : {args.n_novel}")
    print(f"  polya_len : {args.polya_len}")
    print("=" * 60)

    exon_df   = parse_gtf(args.gtf)
    sequences = extract_transcripts(exon_df, args.genome, polya_len=args.polya_len)

    novel_df = pd.DataFrame()
    if args.sqanti and args.n_novel > 0:
        novel_df = load_novel_isoforms(args.sqanti, args.n_novel)

    write_fasta(sequences, f"{args.output}.transcripts.fasta")
    write_filtered_gtf(args.gtf, set(sequences.keys()), f"{args.output}.annotation.gtf")

    novel_tsv = f"{args.output}.novel_isoforms.tsv"
    if novel_df.empty:
        Path(novel_tsv).parent.mkdir(parents=True, exist_ok=True)
        Path(novel_tsv).write_text("isoform\tstructural_category\n")
        print(f"[prepare_reference] Wrote empty novel isoforms TSV: {novel_tsv}")
    else:
        Path(novel_tsv).parent.mkdir(parents=True, exist_ok=True)
        novel_df.to_csv(novel_tsv, sep="\t", index=False)
        print(f"[prepare_reference] Wrote novel isoforms TSV: {novel_tsv}  ({len(novel_df)} rows)")

    print("=" * 60)
    print("Phase 1 complete.")
    print(f"  {args.output}.transcripts.fasta  ({len(sequences)} transcripts)")
    print(f"  {args.output}.annotation.gtf")
    print(f"  {args.output}.novel_isoforms.tsv  ({len(novel_df)} novel isoforms)")
    print("=" * 60)


if __name__ == "__main__":
    main()
