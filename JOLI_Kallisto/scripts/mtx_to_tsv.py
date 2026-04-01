"""
mtx_to_tsv.py

Description:
    Reformats the kallisto quant-tcc output into a canonical abundance.tsv
    with columns:
        transcript_id   tpm
    Only transcripts with TPM > 0 are written.

    Strategy (in order of preference):
      1. Primary: read quant-tcc's own abundance.tsv (target_id, est_counts, tpm)
         and rename target_id → transcript_id.  This avoids any dimension
         mismatch between matrix.abundance.tpm.mtx and transcripts.txt, because
         quant-tcc's abundance.tsv already has the correct transcript IDs matched
         to TPM values. Works identically for long-read and short-read data.
      2. Fallback: parse matrix.abundance.tpm.mtx + transcripts.txt.
         Only used when quant-tcc did not produce abundance.tsv (unexpected).

    Called automatically by run_lr_kallisto.sh as Step 5 after quant-tcc.

Inputs:
    sample_out_dir (str): Path to a per-sample kallisto quant-tcc output
                          directory. Must contain abundance.tsv (preferred)
                          or matrix.abundance.tpm.mtx + transcripts.txt.

Outputs:
    <sample_out_dir>/abundance.tsv  — tab-separated transcript_id + tpm table

Usage:
    python mtx_to_tsv.py <sample_out_dir>
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy.io import mmread


def mtx_to_tsv(sample_out_dir):
    """
    Reformat kallisto quant-tcc output to canonical abundance.tsv.

    Prefers reading quant-tcc's own abundance.tsv (which already has correct
    transcript IDs matched to TPM values) over parsing the sparse MTX file.
    Falls back to MTX parsing only if abundance.tsv is absent.

    Args:
        sample_out_dir (str): Directory containing kallisto quant-tcc outputs.

    Returns:
        str: Path to the written abundance.tsv file.
    """
    output_tsv_path  = os.path.join(sample_out_dir, "abundance.tsv")
    tpm_mtx_path     = os.path.join(sample_out_dir, "matrix.abundance.tpm.mtx")
    transcripts_path = os.path.join(sample_out_dir, "transcripts.txt")

    # --- Primary path: use quant-tcc's own abundance.tsv ---
    # quant-tcc writes abundance.tsv with columns: target_id, est_counts, tpm.
    # Reading it directly avoids the MTX ↔ transcripts.txt dimension mismatch
    # that occurs when quant-tcc quantifies all index transcripts but bustools
    # transcripts.txt only lists transcripts seen in the BUS file.
    if os.path.isfile(output_tsv_path):
        print(f"  Reading quant-tcc abundance.tsv: {output_tsv_path}")
        df = pd.read_csv(output_tsv_path, sep="\t")

        # Normalise column names — quant-tcc uses 'target_id'
        if "target_id" in df.columns:
            df = df.rename(columns={"target_id": "transcript_id"})
        if "transcript_id" not in df.columns:
            raise ValueError(
                f"abundance.tsv has no 'target_id' or 'transcript_id' column. "
                f"Columns found: {list(df.columns)}"
            )
        if "tpm" not in df.columns:
            raise ValueError(
                f"abundance.tsv has no 'tpm' column. Columns found: {list(df.columns)}"
            )

        total_before = len(df)
        df = df[df["tpm"] > 0][["transcript_id", "tpm"]]
        df.to_csv(output_tsv_path, sep="\t", index=False)
        print(f"  Written {len(df)} / {total_before} transcripts (TPM > 0) → {output_tsv_path}")
        return output_tsv_path

    # --- Fallback: parse matrix.abundance.tpm.mtx + transcripts.txt ---
    print(f"  [WARNING] abundance.tsv not found; falling back to MTX parsing.")
    if not os.path.isfile(tpm_mtx_path):
        raise FileNotFoundError(f"TPM matrix not found: {tpm_mtx_path}")
    if not os.path.isfile(transcripts_path):
        raise FileNotFoundError(f"Transcripts file not found: {transcripts_path}")

    print(f"  Reading TPM matrix : {tpm_mtx_path}")
    print(f"  Reading transcripts: {transcripts_path}")

    tpm_matrix  = mmread(tpm_mtx_path)
    transcripts = pd.read_csv(transcripts_path, header=None, sep="\t")
    n_tx        = len(transcripts)

    # quant-tcc matrix shape depends on read type:
    #   long-read  (--long):  (n_transcripts, 1) → rows = transcripts, 1 sample col
    #   short-read (paired):  (2, n_transcripts)  → 2 rows (R1 + R2), cols = transcripts
    # Detect orientation from which dimension matches n_tx, then sum across the
    # read-file axis so we get one TPM value per transcript.
    dense = np.asarray(tpm_matrix.todense())   # (rows, cols) numpy matrix
    print(f"  MTX shape: {dense.shape},  transcripts: {n_tx}")

    if dense.shape[0] == n_tx:
        # (n_transcripts, n_read_files) — average across columns.
        # Each column is a separate TPM estimate (e.g. one per read file);
        # averaging gives one value per transcript that sums to ~1M.
        tpm_values = np.asarray(dense.mean(axis=1)).flatten()
    elif dense.shape[1] == n_tx:
        # (n_read_files, n_transcripts) — average across rows.
        # Short-read paired-end: row 0 = R1 TPM, row 1 = R2 TPM.
        # Average so the result sums to ~1M, not 2M.
        tpm_values = np.asarray(dense.mean(axis=0)).flatten()
    else:
        raise ValueError(
            f"Cannot reconcile MTX shape {dense.shape} with {n_tx} transcripts. "
            f"Neither dimension matches n_tx."
        )

    df = pd.DataFrame({
        "transcript_id": transcripts.iloc[:, 0].values,
        "tpm":           tpm_values,
    })
    df = df[df["tpm"] > 0]
    df.to_csv(output_tsv_path, sep="\t", columns=["transcript_id", "tpm"],
              header=True, index=False)
    print(f"  Written {len(df)} transcripts (TPM > 0) → {output_tsv_path}")
    return output_tsv_path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python mtx_to_tsv.py <sample_out_dir>")
        sys.exit(1)

    sample_out_dir = sys.argv[1]
    mtx_to_tsv(sample_out_dir)
