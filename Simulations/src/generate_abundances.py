"""
generate_abundances.py — Phase 2: generate a synthetic TPM abundance profile.

What it does:
    Reads transcript IDs from a FASTA file and assigns TPM values using one
    of three modes:
      - lognormal : realistic long-tail distribution with optional dropout
      - custom    : user-supplied TPM table, optionally up-regulating cancer genes
      - uniform   : equal TPM for every transcript (benchmarking / null condition)

Inputs:
    --mode         : "lognormal" | "custom" | "uniform"
    --transcripts  : path to .transcripts.fasta (produced by prepare_reference.py)
    --output       : path for output TSV
    --seed         : random seed (default: 42)
    --dropout      : fraction of transcripts set to TPM=0 (lognormal only, default: 0.3)
    --mu           : log-normal mean (default: 2.0)
    --sigma        : log-normal std  (default: 2.0)
    --custom_tpm   : TSV with columns transcript_id, TPM (required for mode=custom)
    --cancer_genes : file with one gene/transcript ID per line to up-regulate
    --fold_change  : multiplier applied to cancer gene TPMs (default: 5.0)

Outputs:
    <output>  — TSV with columns: transcript_id, TPM

Run:
    cd /gpfs/commons/home/atalukder/RNA_Splicing
    conda activate lrgsp_simulation
    python code/Simulations/src/generate_abundances.py \\
        --mode lognormal \\
        --transcripts files/refs/human_sim.transcripts.fasta \\
        --output      files/refs/human_sim.abundances.tsv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ── Transcript ID reader ─────────────────────────────────────────────────────

def read_transcript_ids(fasta_path: str) -> list:
    """
    Parse transcript IDs from a FASTA file (lines starting with '>').

    Args:
        fasta_path (str): Path to the transcripts FASTA.

    Returns:
        list: Ordered list of transcript ID strings.
    """
    ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                # Take everything after '>' up to first whitespace
                ids.append(line[1:].split()[0])
    print(f"[generate_abundances] Read {len(ids)} transcript IDs from {fasta_path}")
    return ids


# ── Abundance modes ──────────────────────────────────────────────────────────

def mode_lognormal(
    tx_ids: list,
    mu: float = 2.0,
    sigma: float = 2.0,
    dropout: float = 0.3,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Assign TPM values drawn from a log-normal distribution with optional dropout.

    Steps:
      1. Draw raw values from LogNormal(mu, sigma).
      2. Zero out a random `dropout` fraction of transcripts.
      3. Renormalise so values sum to 1e6 (TPM convention).

    Args:
        tx_ids  (list):  Transcript ID strings.
        mu      (float): Log-normal mean parameter.
        sigma   (float): Log-normal std parameter.
        dropout (float): Fraction of transcripts to set TPM=0 (0–1).
        seed    (int):   Random seed for reproducibility.

    Returns:
        pd.DataFrame: Columns — transcript_id, TPM.
    """
    rng = np.random.default_rng(seed)
    n   = len(tx_ids)

    raw = rng.lognormal(mean=mu, sigma=sigma, size=n)

    # Apply dropout
    n_drop = int(n * dropout)
    drop_idx = rng.choice(n, size=n_drop, replace=False)
    raw[drop_idx] = 0.0

    # Normalise to TPM (sum = 1e6)
    total = raw.sum()
    tpm = (raw / total * 1e6) if total > 0 else raw

    df = pd.DataFrame({"transcript_id": tx_ids, "TPM": tpm})
    expressed = (df["TPM"] > 0).sum()
    print(f"[generate_abundances] lognormal: {expressed}/{n} transcripts expressed "
          f"(dropout={dropout}, mu={mu}, sigma={sigma})")
    return df


def mode_custom(
    tx_ids: list,
    custom_tpm_file: str,
    cancer_genes_file: str = "",
    fold_change: float = 5.0,
    seed: int = 42,
) -> pd.DataFrame:
    """
    Load a user-supplied TPM table and optionally up-regulate cancer genes.

    Any transcript in tx_ids but not in the custom table is assigned TPM=0.
    After cancer-gene fold-change, values are renormalised to sum=1e6.

    Args:
        tx_ids            (list):  Transcript IDs from the FASTA.
        custom_tpm_file   (str):   Path to TSV with columns transcript_id, TPM.
        cancer_genes_file (str):   Path to file with one gene/transcript ID per line.
        fold_change       (float): Multiplier for cancer gene entries.
        seed              (int):   Unused (for API consistency).

    Returns:
        pd.DataFrame: Columns — transcript_id, TPM.
    """
    print(f"[generate_abundances] Loading custom TPM table: {custom_tpm_file}")
    custom = pd.read_csv(custom_tpm_file, sep="\t")

    required = {"transcript_id", "TPM"}
    missing = required - set(custom.columns)
    if missing:
        raise ValueError(f"Custom TPM file missing columns: {missing}")

    tpm_map = dict(zip(custom["transcript_id"], custom["TPM"]))

    # Up-regulate cancer genes if provided
    if cancer_genes_file:
        print(f"[generate_abundances] Applying fold-change {fold_change}x to cancer genes: {cancer_genes_file}")
        cancer_ids = set(Path(cancer_genes_file).read_text().splitlines())
        for tid in list(tpm_map.keys()):
            if tid in cancer_ids:
                tpm_map[tid] = tpm_map[tid] * fold_change
        print(f"[generate_abundances] Cancer IDs matched in table: "
              f"{sum(1 for t in tpm_map if t in cancer_ids)}")

    tpm = np.array([tpm_map.get(tid, 0.0) for tid in tx_ids], dtype=float)

    # Renormalise
    total = tpm.sum()
    tpm = (tpm / total * 1e6) if total > 0 else tpm

    df = pd.DataFrame({"transcript_id": tx_ids, "TPM": tpm})
    expressed = (df["TPM"] > 0).sum()
    print(f"[generate_abundances] custom: {expressed}/{len(tx_ids)} transcripts expressed")
    return df


def mode_uniform(tx_ids: list) -> pd.DataFrame:
    """
    Assign equal TPM to all transcripts (sum = 1e6).

    Args:
        tx_ids (list): Transcript ID strings.

    Returns:
        pd.DataFrame: Columns — transcript_id, TPM.
    """
    n   = len(tx_ids)
    tpm = np.full(n, 1e6 / n)
    df  = pd.DataFrame({"transcript_id": tx_ids, "TPM": tpm})
    print(f"[generate_abundances] uniform: {n} transcripts, TPM={1e6/n:.4f} each")
    return df


# ── Summary printer ───────────────────────────────────────────────────────────

def print_summary(df: pd.DataFrame) -> None:
    """
    Print a brief summary: total transcripts, expressed count, top 10 by TPM.

    Args:
        df (pd.DataFrame): Columns — transcript_id, TPM.

    Returns:
        None
    """
    expressed = (df["TPM"] > 0).sum()
    print(f"\n[generate_abundances] Summary:")
    print(f"  Total transcripts : {len(df)}")
    print(f"  Expressed (TPM>0) : {expressed}")
    print(f"  Top 10 by TPM:")
    top10 = df.nlargest(10, "TPM")
    for _, row in top10.iterrows():
        print(f"    {row['transcript_id']:40s}  {row['TPM']:.2f}")
    print()


# ── Entry point ──────────────────────────────────────────────────────────────

def main():
    # ============================================================
    # CONFIG — edit these variables before running; do not edit below
    # ============================================================
    default_mode         = "lognormal"
    default_seed         = 42
    default_dropout      = 0.3
    default_mu           = 2.0
    default_sigma        = 2.0
    default_fold_change  = 5.0
    # ============================================================

    parser = argparse.ArgumentParser(
        description="Phase 2: generate synthetic TPM abundance profile."
    )
    parser.add_argument("--mode",         required=True,
                        choices=["lognormal", "custom", "uniform"])
    parser.add_argument("--transcripts",  required=True,
                        help="Path to .transcripts.fasta from Phase 1")
    parser.add_argument("--output",       required=True,
                        help="Output TSV path (transcript_id, TPM)")
    parser.add_argument("--seed",         default=default_seed,        type=int)
    parser.add_argument("--dropout",      default=default_dropout,     type=float)
    parser.add_argument("--mu",           default=default_mu,          type=float)
    parser.add_argument("--sigma",        default=default_sigma,       type=float)
    parser.add_argument("--custom_tpm",   default="",
                        help="TSV with transcript_id, TPM (required for mode=custom)")
    parser.add_argument("--cancer_genes", default="",
                        help="File with one ID per line to up-regulate (mode=custom)")
    parser.add_argument("--fold_change",  default=default_fold_change, type=float)
    args = parser.parse_args()

    print("=" * 60)
    print("Phase 2 — Abundance Generation")
    print(f"  mode         : {args.mode}")
    print(f"  transcripts  : {args.transcripts}")
    print(f"  output       : {args.output}")
    print(f"  seed         : {args.seed}")
    print("=" * 60)

    tx_ids = read_transcript_ids(args.transcripts)

    if args.mode == "lognormal":
        df = mode_lognormal(tx_ids, mu=args.mu, sigma=args.sigma,
                            dropout=args.dropout, seed=args.seed)
    elif args.mode == "custom":
        if not args.custom_tpm:
            raise ValueError("--custom_tpm is required for mode=custom")
        df = mode_custom(tx_ids, args.custom_tpm,
                         cancer_genes_file=args.cancer_genes,
                         fold_change=args.fold_change, seed=args.seed)
    else:  # uniform
        df = mode_uniform(tx_ids)

    print_summary(df)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"[generate_abundances] Wrote: {args.output}")


if __name__ == "__main__":
    main()
