"""
simulate_pacbio.py — Phase 3b: simulate PacBio CCS reads using IsoSeqSim.

What it does:
    1. Converts our abundances TSV (transcript_id, TPM) into IsoSeqSim's
       expression file format: GenePred (GPD) + read-count column.
    2. Runs IsoSeqSim in 'normal' mode with the Sequel completeness profiles
       and PacBio CCS error rates.
    3. Renames outputs to the canonical names used by the rest of the pipeline.

Inputs:
    --reference_prefix : path prefix from Phase 1 (provides .annotation.gtf)
    --genome           : path to genome FASTA (Phase 1 --genome input)
    --abundances       : abundances TSV from Phase 2 (transcript_id, TPM)
    --output_dir       : directory for all output files
    --threads          : CPU threads (default: 16)
    --read_count       : total reads to simulate (default: 10_000_000)
    --seed             : random seed (default: 42)
    --sub_rate         : substitution error rate (default: 0.004)
    --ins_rate         : insertion error rate    (default: 0.006)
    --del_rate         : deletion error rate     (default: 0.006)
    --keep_isoform_ids : flag, preserve original isoform IDs in read names

Outputs (inside --output_dir):
    PacBio.simulated.fasta        — simulated PacBio CCS reads
    PacBio.isoform_counts.tsv     — per-isoform simulated read counts
    PacBio.read_to_isoform.tsv    — read ID → isoform mapping

Run:
    conda activate lrgsp_simulation
    module load minimap2/2.29 SAMtools/1.21
    python code/Simulations/src/simulate_pacbio.py \\
        --reference_prefix files/refs/human_sim \\
        --genome           /path/to/GRCh38.fa \\
        --abundances       files/refs/human_sim.abundances.tsv \\
        --output_dir       files/results/exprmnt_XXXX \\
        --read_count       10000000
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# IsoSeqSim install path (read-only — do not modify)
ISOSEQSIM_PATH = (
    "/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation"
    "/lrgasp-simulation/src/isoseqsim"
)
ISOSEQSIM_BIN      = f"{ISOSEQSIM_PATH}/bin/isoseqsim.py"
ISOSEQSIM_C5_SEQUEL = f"{ISOSEQSIM_PATH}/utilities/5_end_completeness.PacBio-Sequel.tab"
ISOSEQSIM_C3_SEQUEL = f"{ISOSEQSIM_PATH}/utilities/3_end_completeness.PacBio-Sequel.tab"


# ── GTF → GenePred conversion ─────────────────────────────────────────────────

def gtf_to_gpd(gtf_path: str) -> dict:
    """
    Convert a GTF annotation to a dict of GenePred (GPD) strings, keyed by
    transcript_id. Used to build IsoSeqSim's expression profile input.

    GenePred columns (0-based half-open coordinates):
        name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd,
        exonCount, exonStarts, exonEnds

    Args:
        gtf_path (str): Path to GTF annotation file.

    Returns:
        dict: {transcript_id: gpd_line_string}
    """
    print(f"[simulate_pacbio] Converting GTF → GenePred: {gtf_path}")

    # Collect exons per transcript
    tx_exons: dict = {}
    tx_meta: dict  = {}

    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            chrom, _, _, start, end, _, strand, _, attrs = fields
            tx_id   = _extract_attr(attrs, "transcript_id")
            if not tx_id:
                continue
            tx_meta[tx_id] = (chrom, strand)
            tx_exons.setdefault(tx_id, []).append((int(start) - 1, int(end)))  # 0-based half-open

    gpd = {}
    for tx_id, exons in tx_exons.items():
        chrom, strand = tx_meta[tx_id]
        exons_sorted  = sorted(exons, key=lambda e: e[0])
        tx_start      = exons_sorted[0][0]
        tx_end        = exons_sorted[-1][1]
        exon_count    = len(exons_sorted)
        starts        = ",".join(str(e[0]) for e in exons_sorted) + ","
        ends          = ",".join(str(e[1]) for e in exons_sorted) + ","
        gpd[tx_id] = "\t".join([
            tx_id, chrom, strand,
            str(tx_start), str(tx_end),
            str(tx_start), str(tx_end),  # cdsStart = cdsEnd = txStart/txEnd (non-coding)
            str(exon_count),
            starts, ends,
        ])

    print(f"[simulate_pacbio] Built GenePred entries for {len(gpd)} transcripts")
    return gpd


def _extract_attr(attrs: str, key: str) -> str:
    """
    Extract a quoted attribute value from a GTF attributes string.

    Args:
        attrs (str): Raw GTF attributes field.
        key   (str): Attribute name (e.g. 'transcript_id').

    Returns:
        str: Attribute value, or empty string if not found.
    """
    m = re.search(rf'{key}\s+"([^"]+)"', attrs)
    return m.group(1) if m else ""


# ── Abundances → IsoSeqSim expression file ───────────────────────────────────

def build_expr_file(
    abundances_tsv: str,
    gpd: dict,
    total_reads: int,
    out_path: str,
) -> None:
    """
    Convert abundances TSV to IsoSeqSim expression file (GPD + read count).

    Read counts are computed as: count_i = round(TPM_i / sum(TPM) * total_reads).
    Transcripts with no GPD entry (not in GTF) are skipped with a warning.

    Args:
        abundances_tsv (str):  Path to Phase 2 abundances TSV.
        gpd            (dict): {transcript_id: gpd_line} from gtf_to_gpd().
        total_reads    (int):  Total read count to distribute across transcripts.
        out_path       (str):  Destination expression file path.

    Returns:
        None
    """
    print(f"[simulate_pacbio] Building IsoSeqSim expression file → {out_path}")

    abund = pd.read_csv(abundances_tsv, sep="\t")
    tpm_total = abund["TPM"].sum()
    if tpm_total == 0:
        raise ValueError("All TPM values are zero — cannot compute read counts.")

    # Convert TPM → read counts
    abund["read_count"] = (abund["TPM"] / tpm_total * total_reads).round().astype(int)

    skipped = 0
    written = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for _, row in abund.iterrows():
            tid = row["transcript_id"]
            cnt = row["read_count"]
            if tid not in gpd:
                skipped += 1
                continue
            if cnt == 0:
                continue
            fh.write(f"{gpd[tid]}\t{cnt}\n")
            written += 1

    if skipped:
        print(f"[simulate_pacbio] Skipped {skipped} transcripts not found in GTF")
    print(f"[simulate_pacbio] Expression file: {written} transcripts, "
          f"{abund['read_count'].sum():,} total reads")


# ── IsoSeqSim runner ─────────────────────────────────────────────────────────

def run_isoseqsim(
    genome: str,
    gtf: str,
    expr_file: str,
    output_fasta: str,
    output_transcript: str,
    tempdir: str,
    threads: int,
    read_count_millions: float,
    sub_rate: float,
    ins_rate: float,
    del_rate: float,
    keep_isoform_ids: bool,
    seed: int,
) -> None:
    """
    Run IsoSeqSim in 'normal' mode with Sequel completeness profiles.

    Args:
        genome              (str):   Path to genome FASTA.
        gtf                 (str):   Path to GTF annotation.
        expr_file           (str):   Path to expression file (GPD + count).
        output_fasta        (str):   Path for output FASTA reads.
        output_transcript   (str):   Path for output transcript annotation.
        tempdir             (str):   Temp directory for intermediate files.
        threads             (int):   CPU threads.
        read_count_millions (float): Total reads in millions.
        sub_rate            (float): Substitution error rate.
        ins_rate            (float): Insertion error rate.
        del_rate            (float): Deletion error rate.
        keep_isoform_ids    (bool):  Preserve original isoform IDs in read names.
        seed                (int):   Unused by IsoSeqSim directly; kept for API consistency.

    Returns:
        None
    """
    import subprocess
    cmd = [
        "python", ISOSEQSIM_BIN,
        "-g", genome,
        "-a", gtf,
        "--expr", expr_file,
        "-o", output_fasta,
        "-t", output_transcript,
        "--tempdir", tempdir,
        "--cpu", str(threads),
        "--read_number", str(read_count_millions),
        "--c5", ISOSEQSIM_C5_SEQUEL,
        "--c3", ISOSEQSIM_C3_SEQUEL,
        "--es", str(sub_rate),
        "--ei", str(ins_rate),
        "--ed", str(del_rate),
        "-m", "normal",
    ]
    if keep_isoform_ids:
        cmd.append("--keep_isoform_ids")

    print(f"[simulate_pacbio] Running IsoSeqSim ({read_count_millions}M reads)...")
    print(f"  Command: {' '.join(cmd)}")
    _run(cmd)
    print(f"[simulate_pacbio] IsoSeqSim complete")


# ── Output processing ─────────────────────────────────────────────────────────

def parse_transcript_output(transcript_file: str, out_counts: str, out_r2i: str) -> None:
    """
    Parse IsoSeqSim's transcript annotation output into isoform_counts and
    read_to_isoform TSVs.

    IsoSeqSim's -t output is a modified GenePred table where the last column
    is the simulated read count for that isoform.

    Args:
        transcript_file (str): Path to IsoSeqSim -t output file.
        out_counts      (str): Destination for PacBio.isoform_counts.tsv.
        out_r2i         (str): Destination for PacBio.read_to_isoform.tsv.

    Returns:
        None
    """
    rows = []
    with open(transcript_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            tx_id      = parts[0]
            read_count = int(parts[-1]) if parts[-1].isdigit() else 0
            rows.append({"transcript_id": tx_id, "simulated_count": read_count})

    df = pd.DataFrame(rows)

    # isoform_counts — all transcripts
    df.to_csv(out_counts, sep="\t", index=False)
    print(f"[simulate_pacbio] Wrote isoform counts: {out_counts} ({len(df)} rows)")

    # read_to_isoform — only expressed transcripts
    expressed = df[df["simulated_count"] > 0]
    expressed.to_csv(out_r2i, sep="\t", index=False)
    print(f"[simulate_pacbio] Wrote read-to-isoform: {out_r2i} ({len(expressed)} entries)")


# ── Subprocess helper ─────────────────────────────────────────────────────────

def _run(cmd: list) -> None:
    """
    Run a shell command, streaming output to terminal.

    Args:
        cmd (list): Command and arguments.

    Returns:
        None

    Raises:
        SystemExit: If the command returns a non-zero exit code.
    """
    import subprocess
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        print(f"[simulate_pacbio] ERROR: command failed (exit {result.returncode})")
        sys.exit(result.returncode)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    # ============================================================
    # CONFIG — edit these variables before running; do not edit below
    # ============================================================
    default_threads    = 16
    default_read_count = 10_000_000
    default_seed       = 42
    default_sub_rate   = 0.004
    default_ins_rate   = 0.006
    default_del_rate   = 0.006
    # ============================================================

    parser = argparse.ArgumentParser(
        description="Phase 3b: simulate PacBio CCS reads with IsoSeqSim."
    )
    parser.add_argument("--reference_prefix", required=True,
                        help="Path prefix from Phase 1 (provides .annotation.gtf)")
    parser.add_argument("--genome",           required=True,
                        help="Path to genome FASTA (same as Phase 1 --genome)")
    parser.add_argument("--abundances",       required=True,
                        help="Abundances TSV from Phase 2 (transcript_id, TPM)")
    parser.add_argument("--output_dir",       required=True)
    parser.add_argument("--threads",          default=default_threads,    type=int)
    parser.add_argument("--read_count",       default=default_read_count, type=int)
    parser.add_argument("--seed",             default=default_seed,       type=int)
    parser.add_argument("--sub_rate",         default=default_sub_rate,   type=float)
    parser.add_argument("--ins_rate",         default=default_ins_rate,   type=float)
    parser.add_argument("--del_rate",         default=default_del_rate,   type=float)
    parser.add_argument("--keep_isoform_ids", action="store_true")
    args = parser.parse_args()

    out_dir  = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    gtf_path = f"{args.reference_prefix}.annotation.gtf"
    read_count_millions = args.read_count / 1_000_000

    print("=" * 60)
    print("Phase 3b — PacBio Read Simulation (IsoSeqSim)")
    print(f"  reference_prefix : {args.reference_prefix}")
    print(f"  genome           : {args.genome}")
    print(f"  abundances       : {args.abundances}")
    print(f"  output_dir       : {args.output_dir}")
    print(f"  read_count       : {args.read_count:,}")
    print(f"  threads          : {args.threads}")
    print(f"  error rates      : sub={args.sub_rate}, ins={args.ins_rate}, del={args.del_rate}")
    print("=" * 60)

    # Step 1 — GTF → GenePred
    gpd = gtf_to_gpd(gtf_path)

    # Step 2 — Build expression file
    expr_file = str(out_dir / "pacbio_expr.txt")
    build_expr_file(args.abundances, gpd, args.read_count, expr_file)

    # Step 3 — Run IsoSeqSim
    output_fasta      = str(out_dir / "PacBio.simulated.fasta")
    output_transcript = str(out_dir / "pacbio_transcript.txt")
    tempdir           = str(out_dir / "isoseqsim_tmp")

    run_isoseqsim(
        genome=args.genome,
        gtf=gtf_path,
        expr_file=expr_file,
        output_fasta=output_fasta,
        output_transcript=output_transcript,
        tempdir=tempdir,
        threads=args.threads,
        read_count_millions=read_count_millions,
        sub_rate=args.sub_rate,
        ins_rate=args.ins_rate,
        del_rate=args.del_rate,
        keep_isoform_ids=args.keep_isoform_ids,
        seed=args.seed,
    )

    # Step 4 — Parse transcript output → isoform counts + read-to-isoform
    parse_transcript_output(
        output_transcript,
        str(out_dir / "PacBio.isoform_counts.tsv"),
        str(out_dir / "PacBio.read_to_isoform.tsv"),
    )

    print("=" * 60)
    print("Phase 3b complete. Outputs:")
    for name in ["PacBio.simulated.fasta", "PacBio.isoform_counts.tsv",
                 "PacBio.read_to_isoform.tsv"]:
        p = out_dir / name
        print(f"  [{'OK' if p.exists() else 'MISSING'}] {p}")
    print("=" * 60)


if __name__ == "__main__":
    main()
