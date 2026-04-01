"""
simulate_ont.py — Phase 3c: simulate ONT reads using Trans-NanoSim.

What it does:
    1. Builds an expression profile file from our abundances TSV in the format
       NanoSim expects (transcript_id, col2, TPM with header line).
    2. Runs NanoSim in transcriptome mode using a pre-trained human cDNA or
       dRNA error model.
    3. Renames outputs to the canonical names used by the rest of the pipeline.

Inputs:
    --reference_prefix : path prefix from Phase 1 (provides .transcripts.fasta)
    --abundances       : abundances TSV from Phase 2 (transcript_id, TPM)
    --output_dir       : directory for all output files
    --ont_type         : "cDNA" or "dRNA" (default: "cDNA")
    --threads          : CPU threads (default: 16)
    --read_count       : total reads to simulate (default: 30_000_000)
    --seed             : random seed (default: 42)
    --noise_reads      : flag — include unaligned/noise reads in output

Outputs (inside --output_dir):
    ONT.simulated.fastq          — simulated ONT reads (FASTQ)
    ONT.isoform_counts.tsv       — per-isoform simulated read counts
    ONT.read_to_isoform.tsv      — read ID → isoform mapping

Pre-trained models used (read-only):
    cDNA: .../NanoSim/pre-trained_models/human_NA12878_cDNA_Bham1_guppy/training
    dRNA: .../NanoSim/pre-trained_models/human_NA12878_dRNA_Bham1_guppy/training

⚠ DEPENDENCY NOTE:
    NanoSim requires the 'HTSeq' package, which is NOT in the lrgsp_simulation env.
    To install it, create a new env:
        conda create -n lrgsp_simulation_2 --clone lrgsp_simulation
        conda activate lrgsp_simulation_2
        pip install HTSeq
        conda env export > env/lrgsp_simulation_2.yml

Run:
    conda activate lrgsp_simulation_2   # env with HTSeq
    python code/Simulations/src/simulate_ont.py \\
        --reference_prefix files/refs/human_sim \\
        --abundances       files/refs/human_sim.abundances.tsv \\
        --output_dir       files/results/exprmnt_XXXX \\
        --ont_type         cDNA \\
        --read_count       30000000
"""

import sys
from pathlib import Path

import pandas as pd


# NanoSim paths (read-only — do not modify)
_NANOSIM_BASE = (
    "/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation"
    "/lrgasp-simulation/src/NanoSim"
)
NANOSIM_BIN = f"{_NANOSIM_BASE}/src/simulator.py"
NANOSIM_MODELS = {
    "cDNA": f"{_NANOSIM_BASE}/pre-trained_models/human_NA12878_cDNA_Bham1_guppy/training",
    "dRNA": f"{_NANOSIM_BASE}/pre-trained_models/human_NA12878_dRNA_Bham1_guppy/training",
}


# ── Pre-flight checks ─────────────────────────────────────────────────────────

def check_htseq() -> None:
    """
    Verify HTSeq is importable. Raises SystemExit with installation instructions
    if it is not available (NanoSim requires it).

    Returns:
        None
    """
    try:
        import HTSeq  # noqa: F401
        print("[simulate_ont] HTSeq found — OK")
    except ImportError:
        print(
            "[simulate_ont] ERROR: HTSeq is not installed in the current environment.\n"
            "  NanoSim requires HTSeq. To fix, create a new env:\n"
            "    conda create -n lrgsp_simulation_2 --clone lrgsp_simulation\n"
            "    conda activate lrgsp_simulation_2\n"
            "    pip install HTSeq\n"
            "    conda env export > env/lrgsp_simulation_2.yml\n"
            "  Then re-run with: conda activate lrgsp_simulation_2"
        )
        sys.exit(1)


def check_model(ont_type: str) -> str:
    """
    Verify the pre-trained NanoSim model exists for the given ONT type.

    Args:
        ont_type (str): "cDNA" or "dRNA".

    Returns:
        str: Path to the model prefix.
    """
    if ont_type not in NANOSIM_MODELS:
        print(f"[simulate_ont] ERROR: unknown ont_type '{ont_type}'. Use 'cDNA' or 'dRNA'.")
        sys.exit(1)
    model_prefix = NANOSIM_MODELS[ont_type]
    profile = model_prefix + "_model_profile"
    if not Path(profile).exists():
        print(f"[simulate_ont] ERROR: NanoSim model not found: {profile}")
        print(f"  Expected at: {_NANOSIM_BASE}/pre-trained_models/")
        sys.exit(1)
    print(f"[simulate_ont] Using pre-trained model: {model_prefix}")
    return model_prefix


# ── Expression profile builder ────────────────────────────────────────────────

def build_exp_file(abundances_tsv: str, out_path: str) -> None:
    """
    Convert abundances TSV to NanoSim expression profile format.

    NanoSim format (tab-separated, with header):
        transcript_id  <col2>  TPM
    Only transcripts with TPM > 0 are included.
    NanoSim internally strips version suffixes (splits on '.'), so IDs must
    match the transcript FASTA headers after that strip.

    Args:
        abundances_tsv (str): Path to Phase 2 abundances TSV.
        out_path       (str): Destination expression profile file path.

    Returns:
        None
    """
    print(f"[simulate_ont] Building NanoSim expression profile from {abundances_tsv}")
    df = pd.read_csv(abundances_tsv, sep="\t")
    expressed = df[df["TPM"] > 0].copy()

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("transcript_id\tcount\tTPM\n")
        for _, row in expressed.iterrows():
            fh.write(f"{row['transcript_id']}\t0\t{row['TPM']:.6f}\n")

    print(f"[simulate_ont] Expression profile: {len(expressed)} transcripts with TPM>0 "
          f"→ {out_path}")


# ── NanoSim runner ────────────────────────────────────────────────────────────

def run_nanosim(
    transcript_fasta: str,
    exp_file: str,
    model_prefix: str,
    output_prefix: str,
    read_count: int,
    ont_type: str,
    threads: int,
    seed: int,
    noise_reads: bool,
) -> None:
    """
    Run NanoSim in transcriptome mode.

    Args:
        transcript_fasta (str):  Path to .transcripts.fasta from Phase 1.
        exp_file         (str):  Path to expression profile file.
        model_prefix     (str):  Pre-trained model prefix.
        output_prefix    (str):  Prefix for NanoSim output files.
        read_count       (int):  Total reads to simulate.
        ont_type         (str):  "cDNA" or "dRNA".
        threads          (int):  CPU threads.
        seed             (int):  Random seed.
        noise_reads      (bool): If False, pass --aligned_only to suppress noise reads.

    Returns:
        None
    """
    import subprocess

    basecaller = "guppy"
    read_type  = "dRNA" if ont_type == "dRNA" else "cDNA"

    cmd = [
        "python", NANOSIM_BIN, "transcriptome",
        "-rt",  transcript_fasta,
        "-e",   exp_file,
        "-c",   model_prefix,
        "-o",   output_prefix,
        "-n",   str(read_count),
        "-b",   basecaller,
        "-r",   read_type,
        "-t",   str(threads),
        "--seed", str(seed),
        "--fastq",
        "--no_model_ir",   # skip intron retention (requires genome ref; not provided here)
    ]
    if not noise_reads:
        cmd.append("--aligned_only")

    print(f"[simulate_ont] Running NanoSim ({read_count:,} reads, {ont_type})...")
    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        print(f"[simulate_ont] ERROR: NanoSim failed (exit {result.returncode})")
        sys.exit(result.returncode)
    print(f"[simulate_ont] NanoSim complete")


# ── Output renaming and post-processing ──────────────────────────────────────

def rename_outputs(nanosim_prefix: str, output_dir: str, noise_reads: bool) -> None:
    """
    Rename NanoSim output files to the canonical pipeline names.

    NanoSim produces:
      <prefix>_aligned_reads.fastq
      <prefix>_unaligned_reads.fastq  (if noise_reads=True)
      <prefix>_aligned_error_profile

    Pipeline expects: ONT.simulated.fastq (aligned reads only, or merged if noise_reads)

    Args:
        nanosim_prefix (str):  Prefix passed to NanoSim -o.
        output_dir     (str):  Output directory path.
        noise_reads    (bool): Whether unaligned reads were also simulated.

    Returns:
        None
    """
    out = Path(output_dir)
    aligned   = Path(f"{nanosim_prefix}_aligned_reads.fastq")
    unaligned = Path(f"{nanosim_prefix}_unaligned_reads.fastq")
    canonical = out / "ONT.simulated.fastq"

    if noise_reads and unaligned.exists():
        # Merge aligned + unaligned into single output
        print(f"[simulate_ont] Merging aligned + unaligned reads → {canonical}")
        with open(canonical, "w") as fout:
            for src in [aligned, unaligned]:
                if src.exists():
                    fout.write(src.read_text())
        aligned.unlink(missing_ok=True)
        unaligned.unlink(missing_ok=True)
    elif aligned.exists():
        aligned.rename(canonical)
        print(f"[simulate_ont] Renamed: {aligned.name} → {canonical.name}")
    else:
        print(f"[simulate_ont] [WARN] Expected aligned reads not found: {aligned}")


def build_isoform_counts(exp_file: str, canonical_fasta: str,
                         out_counts: str, out_r2i: str) -> None:
    """
    Parse NanoSim's FASTQ headers to count simulated reads per transcript,
    then write isoform_counts and read_to_isoform TSVs.

    NanoSim FASTQ header format: @{transcript_id}_{read_idx}_{aligned|unaligned}_{N}

    Transcript IDs are resolved by matching FASTQ headers against the known
    transcript list from the expression file (longest-prefix match), which
    correctly handles transcript IDs that contain underscores.

    Args:
        exp_file        (str): Expression profile used (to get full transcript list).
        canonical_fasta (str): Path to ONT.simulated.fastq.
        out_counts      (str): Destination for ONT.isoform_counts.tsv.
        out_r2i         (str): Destination for ONT.read_to_isoform.tsv.

    Returns:
        None
    """
    from collections import Counter

    # Load known transcript IDs from expression file, sorted longest-first
    # so we always match the most specific (longest) ID first.
    exp_df = pd.read_csv(exp_file, sep="\t")
    known_ids = sorted(exp_df["transcript_id"].tolist(), key=len, reverse=True)

    def _parse_tx_id(header: str) -> str:
        """Match header to the known transcript ID it starts with."""
        for tid in known_ids:
            if header == tid or header.startswith(tid + "_"):
                return tid
        # Fallback: strip last three underscore-separated tokens (read_idx, status, N)
        parts = header.rsplit("_", 3)
        return parts[0] if len(parts) >= 3 else header

    # Count reads per transcript from FASTQ headers
    counts: Counter = Counter()
    read_records = []
    with open(canonical_fasta) as fh:
        for line in fh:
            if line.startswith("@"):
                header = line[1:].split()[0]
                tx_id  = _parse_tx_id(header)
                counts[tx_id] += 1
                read_records.append({"read_id": header, "transcript_id": tx_id})

    # Build counts TSV from expression profile (all expressed transcripts)
    exp_df = pd.read_csv(exp_file, sep="\t")
    rows = []
    for _, row in exp_df.iterrows():
        tid = row["transcript_id"]
        rows.append({"transcript_id": tid, "simulated_count": counts.get(tid, 0)})
    counts_df = pd.DataFrame(rows)
    counts_df.to_csv(out_counts, sep="\t", index=False)
    print(f"[simulate_ont] Wrote isoform counts: {out_counts} ({len(counts_df)} rows)")

    # Build read-to-isoform TSV
    r2i_df = pd.DataFrame(read_records)
    r2i_df.to_csv(out_r2i, sep="\t", index=False)
    print(f"[simulate_ont] Wrote read-to-isoform: {out_r2i} ({len(r2i_df)} entries)")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    import argparse
    # ============================================================
    # CONFIG — edit these variables before running; do not edit below
    # ============================================================
    default_ont_type   = "cDNA"
    default_threads    = 16
    default_read_count = 30_000_000
    default_seed       = 42
    # ============================================================

    parser = argparse.ArgumentParser(
        description="Phase 3c: simulate ONT reads with Trans-NanoSim."
    )
    parser.add_argument("--reference_prefix", required=True,
                        help="Path prefix from Phase 1 (provides .transcripts.fasta)")
    parser.add_argument("--abundances",       required=True,
                        help="Abundances TSV from Phase 2 (transcript_id, TPM)")
    parser.add_argument("--output_dir",       required=True)
    parser.add_argument("--ont_type",         default=default_ont_type,
                        choices=["cDNA", "dRNA"])
    parser.add_argument("--threads",          default=default_threads,    type=int)
    parser.add_argument("--read_count",       default=default_read_count, type=int)
    parser.add_argument("--seed",             default=default_seed,       type=int)
    parser.add_argument("--noise_reads",      action="store_true",
                        help="Include unaligned/noise reads in output")
    args = parser.parse_args()

    out_dir           = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    transcript_fasta  = f"{args.reference_prefix}.transcripts.fasta"
    nanosim_out_prefix = str(out_dir / "ont_raw")

    print("=" * 60)
    print("Phase 3c — ONT Read Simulation (NanoSim)")
    print(f"  reference_prefix : {args.reference_prefix}")
    print(f"  abundances       : {args.abundances}")
    print(f"  output_dir       : {args.output_dir}")
    print(f"  ont_type         : {args.ont_type}")
    print(f"  read_count       : {args.read_count:,}")
    print(f"  threads          : {args.threads}")
    print(f"  seed             : {args.seed}")
    print(f"  noise_reads      : {args.noise_reads}")
    print("=" * 60)

    # Step 0 — pre-flight checks
    check_htseq()
    model_prefix = check_model(args.ont_type)

    # Step 1 — build expression profile
    exp_file = str(out_dir / "ont_expression.tsv")
    build_exp_file(args.abundances, exp_file)

    # Step 2 — run NanoSim
    run_nanosim(
        transcript_fasta=transcript_fasta,
        exp_file=exp_file,
        model_prefix=model_prefix,
        output_prefix=nanosim_out_prefix,
        read_count=args.read_count,
        ont_type=args.ont_type,
        threads=args.threads,
        seed=args.seed,
        noise_reads=args.noise_reads,
    )

    # Step 3 — rename outputs
    rename_outputs(nanosim_out_prefix, str(out_dir), noise_reads=args.noise_reads)

    # Step 4 — isoform counts + read-to-isoform
    canonical_fastq = str(out_dir / "ONT.simulated.fastq")
    if Path(canonical_fastq).exists():
        build_isoform_counts(
            exp_file, canonical_fastq,
            str(out_dir / "ONT.isoform_counts.tsv"),
            str(out_dir / "ONT.read_to_isoform.tsv"),
        )

    print("=" * 60)
    print("Phase 3c complete. Outputs:")
    for name in ["ONT.simulated.fastq", "ONT.isoform_counts.tsv",
                 "ONT.read_to_isoform.tsv"]:
        p = out_dir / name
        print(f"  [{'OK' if p.exists() else 'MISSING'}] {p}")
    print("=" * 60)


if __name__ == "__main__":
    main()
