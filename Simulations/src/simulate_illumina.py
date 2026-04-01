"""
simulate_illumina.py — Phase 3a: simulate paired-end Illumina reads using RSEM.

What it does:
    1. Builds a RSEM reference from the transcript FASTA (Phase 1 output).
    2. Converts the abundances TSV (Phase 2 output) into RSEM isoforms.results format.
    3. Runs rsem-simulate-reads to produce paired-end FASTQ files.
    4. Renames outputs to the canonical names used by the rest of the pipeline.

Inputs:
    --reference_prefix : path prefix from Phase 1 (e.g. files/refs/human_sim)
    --abundances       : path to abundances TSV from Phase 2 (transcript_id, TPM)
    --output_dir       : directory for all output files
    --model_file       : path to RSEM *.model file (learned from real data via
                         rsem-calculate-expression; see note below)
    --threads          : CPU threads (default: 16)
    --read_count       : read pairs to simulate (default: 100_000_000)
    --seed             : random seed (default: 42)
    --theta0           : background noise fraction 0–1 (default: 0.05)

Outputs (inside --output_dir):
    Illumina.simulated_1.fq          — forward reads
    Illumina.simulated_2.fq          — reverse reads
    Illumina.isoform_counts.tsv      — per-isoform simulated read counts
    Illumina.read_to_isoform.tsv     — read ID → isoform ID mapping (from RSEM header)

NOTE — RSEM model file:
    rsem-simulate-reads requires a model file that encodes the sequencing error profile,
    fragment length distribution, and read start position distribution. This file is
    produced by running rsem-calculate-expression on a real Illumina RNA-seq dataset:

        rsem-calculate-expression --paired-end \\
            read_1.fastq read_2.fastq \\
            <rsem_reference_prefix> sample_name

    The model file will be at: sample_name.stat/sample_name.model
    Pass that path as --model_file.

Run:
    conda activate lrgsp_simulation
    module load RSEM/1.3.3-foss-2022a SAMtools/1.21
    python code/Simulations/src/simulate_illumina.py \\
        --reference_prefix files/refs/human_sim \\
        --abundances       files/refs/human_sim.abundances.tsv \\
        --output_dir       files/results/exprmnt_XXXX \\
        --model_file       /path/to/sample.stat/sample.model \\
        --read_count       100000000
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


# ── RSEM reference builder ────────────────────────────────────────────────────

def prepare_rsem_reference(transcript_fasta: str, rsem_ref_prefix: str, threads: int) -> None:
    """
    Build a RSEM reference from a transcript FASTA using rsem-prepare-reference.

    Args:
        transcript_fasta (str): Path to .transcripts.fasta from Phase 1.
        rsem_ref_prefix  (str): Output prefix for RSEM reference files.
        threads          (int): Number of CPU threads.

    Returns:
        None
    """
    # Skip if reference already built (idempotent)
    if Path(f"{rsem_ref_prefix}.transcripts.fa").exists():
        print(f"[simulate_illumina] RSEM reference already exists, skipping: {rsem_ref_prefix}")
        return

    cmd = [
        "rsem-prepare-reference",
        "--num-threads", str(threads),
        transcript_fasta,
        rsem_ref_prefix,
    ]
    print(f"[simulate_illumina] Building RSEM reference: {' '.join(cmd)}")
    _run(cmd)
    print(f"[simulate_illumina] RSEM reference built: {rsem_ref_prefix}")


# ── Abundances → RSEM isoforms.results ───────────────────────────────────────

def abundances_to_rsem_results(
    abundances_tsv: str,
    transcript_fasta: str,
    out_path: str,
) -> None:
    """
    Convert our abundances TSV (transcript_id, TPM) into RSEM isoforms.results format.

    RSEM's simulator only reads the TPM column, but the full file format is required.
    Lengths are computed from the transcript FASTA. Other numeric columns are set to 0.

    Args:
        abundances_tsv   (str): Path to Phase 2 abundances TSV.
        transcript_fasta (str): Path to Phase 1 .transcripts.fasta (for lengths).
        out_path         (str): Destination isoforms.results file path.

    Returns:
        None
    """
    print(f"[simulate_illumina] Building RSEM isoforms.results from {abundances_tsv}")

    # Read abundances
    abund = pd.read_csv(abundances_tsv, sep="\t")

    # Compute transcript lengths from FASTA
    lengths = {}
    current_id = None
    current_len = 0
    with open(transcript_fasta) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    lengths[current_id] = current_len
                current_id = line[1:].split()[0]
                current_len = 0
            else:
                current_len += len(line)
        if current_id is not None:
            lengths[current_id] = current_len

    print(f"[simulate_illumina] Computed lengths for {len(lengths)} transcripts")

    # Build isoforms.results table
    rows = []
    for _, row in abund.iterrows():
        tid = row["transcript_id"]
        tpm = row["TPM"]
        length = lengths.get(tid, 0)
        rows.append({
            "transcript_id":    tid,
            "gene_id":          tid,          # use transcript_id as gene_id placeholder
            "length":           length,
            "effective_length": max(1, length),
            "expected_count":   0.0,
            "TPM":              tpm,
            "FPKM":             0.0,
            "IsoPct":           0.0,
        })

    df = pd.DataFrame(rows, columns=[
        "transcript_id", "gene_id", "length", "effective_length",
        "expected_count", "TPM", "FPKM", "IsoPct",
    ])
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"[simulate_illumina] Wrote RSEM isoforms.results: {out_path} ({len(df)} rows)")


# ── RSEM simulation ───────────────────────────────────────────────────────────

def run_rsem_simulate(
    rsem_ref_prefix: str,
    model_file: str,
    isoforms_results: str,
    output_prefix: str,
    read_count: int,
    theta0: float,
    seed: int,
) -> None:
    """
    Run rsem-simulate-reads to generate paired-end FASTQ files.

    Args:
        rsem_ref_prefix  (str):   RSEM reference prefix (from prepare_rsem_reference).
        model_file       (str):   Path to RSEM *.model file (from real data run).
        isoforms_results (str):   Path to isoforms.results file.
        output_prefix    (str):   Prefix for rsem-simulate-reads output files.
        read_count       (int):   Total read pairs to simulate.
        theta0           (float): Fraction of background noise reads (0–1).
        seed             (int):   Random seed.

    Returns:
        None
    """
    cmd = [
        "rsem-simulate-reads",
        rsem_ref_prefix,
        model_file,
        isoforms_results,
        str(theta0),
        str(read_count),
        output_prefix,
        "--seed", str(seed),
    ]
    print(f"[simulate_illumina] Running rsem-simulate-reads ({read_count:,} reads)...")
    print(f"  Command: {' '.join(cmd)}")
    _run(cmd)
    print(f"[simulate_illumina] rsem-simulate-reads complete")


# ── Output renaming and summary ───────────────────────────────────────────────

def rename_outputs(rsem_prefix: str, output_dir: str) -> None:
    """
    Rename RSEM output files to the canonical names used by this pipeline.

    RSEM produces: <prefix>_1.fq, <prefix>_2.fq, <prefix>.sim.isoforms.results
    Pipeline expects: Illumina.simulated_1.fq, Illumina.simulated_2.fq,
                      Illumina.isoform_counts.tsv

    Args:
        rsem_prefix (str): Prefix passed to rsem-simulate-reads.
        output_dir  (str): Directory containing the output files.

    Returns:
        None
    """
    out = Path(output_dir)

    renames = {
        f"{rsem_prefix}_1.fq":              out / "Illumina.simulated_1.fq",
        f"{rsem_prefix}_2.fq":              out / "Illumina.simulated_2.fq",
        f"{rsem_prefix}.sim.isoforms.results": out / "Illumina.isoform_counts.tsv",
    }
    for src_name, dst in renames.items():
        src = Path(src_name)
        if src.exists():
            src.rename(dst)
            print(f"[simulate_illumina] Renamed: {src.name} → {dst.name}")
        else:
            print(f"[simulate_illumina] [WARN] Expected output not found: {src}")


def write_read_to_isoform(isoform_counts_tsv: str, out_path: str) -> None:
    """
    Build a read_to_isoform TSV from RSEM isoform counts (approximation).

    Since RSEM does not produce a per-read assignment file directly, this writes
    a summary table of (transcript_id, simulated_count) from the sim.isoforms.results.

    Args:
        isoform_counts_tsv (str): Path to Illumina.isoform_counts.tsv.
        out_path           (str): Destination TSV path.

    Returns:
        None
    """
    df = pd.read_csv(isoform_counts_tsv, sep="\t")
    # expected_count in .sim.isoforms.results is the simulated read count
    summary = df[["transcript_id", "expected_count"]].rename(
        columns={"expected_count": "simulated_count"}
    )
    summary = summary[summary["simulated_count"] > 0]
    summary.to_csv(out_path, sep="\t", index=False)
    print(f"[simulate_illumina] Wrote read-to-isoform summary: {out_path} ({len(summary)} entries)")


# ── Subprocess helper ─────────────────────────────────────────────────────────

def _run(cmd: list) -> None:
    """
    Run a shell command, streaming stdout/stderr to the terminal.

    Args:
        cmd (list): Command and arguments as a list of strings.

    Returns:
        None

    Raises:
        SystemExit: If the command returns a non-zero exit code.
    """
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        print(f"[simulate_illumina] ERROR: command failed (exit {result.returncode})")
        print(f"  {' '.join(cmd)}")
        sys.exit(result.returncode)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    # ============================================================
    # CONFIG — edit these variables before running; do not edit below
    # ============================================================
    default_threads    = 16
    default_read_count = 100_000_000
    default_seed       = 42
    default_theta0     = 0.05       # 5% background noise reads
    # ============================================================

    parser = argparse.ArgumentParser(
        description="Phase 3a: simulate paired-end Illumina reads with RSEM."
    )
    parser.add_argument("--reference_prefix", required=True,
                        help="Path prefix from Phase 1 (e.g. files/refs/human_sim)")
    parser.add_argument("--abundances",       required=True,
                        help="Abundances TSV from Phase 2 (transcript_id, TPM)")
    parser.add_argument("--output_dir",       required=True,
                        help="Directory for output files")
    parser.add_argument("--model_file",       required=True,
                        help="RSEM *.model file from rsem-calculate-expression on real data")
    parser.add_argument("--threads",     default=default_threads,    type=int)
    parser.add_argument("--read_count",  default=default_read_count, type=int)
    parser.add_argument("--seed",        default=default_seed,       type=int)
    parser.add_argument("--theta0",      default=default_theta0,     type=float)
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    transcript_fasta = f"{args.reference_prefix}.transcripts.fasta"
    rsem_ref_prefix  = str(out_dir / "rsem_ref" / "ref")
    rsem_out_prefix  = str(out_dir / "Illumina_raw")
    isoforms_results = str(out_dir / "rsem_input.isoforms.results")

    print("=" * 60)
    print("Phase 3a — Illumina Read Simulation (RSEM)")
    print(f"  reference_prefix : {args.reference_prefix}")
    print(f"  abundances       : {args.abundances}")
    print(f"  output_dir       : {args.output_dir}")
    print(f"  model_file       : {args.model_file}")
    print(f"  read_count       : {args.read_count:,}")
    print(f"  threads          : {args.threads}")
    print(f"  seed             : {args.seed}")
    print(f"  theta0           : {args.theta0}")
    print("=" * 60)

    # Step 1 — Build RSEM reference
    Path(rsem_ref_prefix).parent.mkdir(parents=True, exist_ok=True)
    prepare_rsem_reference(transcript_fasta, rsem_ref_prefix, args.threads)

    # Step 2 — Convert abundances to RSEM isoforms.results
    abundances_to_rsem_results(args.abundances, transcript_fasta, isoforms_results)

    # Step 3 — Simulate reads
    run_rsem_simulate(
        rsem_ref_prefix, args.model_file, isoforms_results,
        rsem_out_prefix, args.read_count, args.theta0, args.seed,
    )

    # Step 4 — Rename outputs to canonical names
    rename_outputs(rsem_out_prefix, str(out_dir))

    # Step 5 — Write read-to-isoform summary
    isoform_counts = str(out_dir / "Illumina.isoform_counts.tsv")
    if Path(isoform_counts).exists():
        write_read_to_isoform(isoform_counts, str(out_dir / "Illumina.read_to_isoform.tsv"))

    print("=" * 60)
    print("Phase 3a complete. Outputs:")
    for name in ["Illumina.simulated_1.fq", "Illumina.simulated_2.fq",
                 "Illumina.isoform_counts.tsv", "Illumina.read_to_isoform.tsv"]:
        p = out_dir / name
        status = "OK" if p.exists() else "MISSING"
        print(f"  [{status}] {p}")
    print("=" * 60)


if __name__ == "__main__":
    main()
