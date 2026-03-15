"""
mtx_to_tsv.py

Description:
    Converts the kallisto quant-tcc output matrix (matrix.abundance.tpm.mtx)
    into a human-readable TSV file (abundance.tsv) with columns:
        transcript_id   tpm
    Only transcripts with TPM > 0 are written (zero-expression transcripts
    are dropped to keep the file compact).

    Called automatically by run_lr_kallisto.sh as Step 5 after quant-tcc.

Inputs:
    sample_out_dir (str): Path to a per-sample kallisto output directory
                          containing:
                          - matrix.abundance.tpm.mtx  (sparse TPM matrix)
                          - transcripts.txt            (transcript ID list)

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
    Convert matrix.abundance.tpm.mtx + transcripts.txt → abundance.tsv.

    Args:
        sample_out_dir (str): Directory containing the kallisto quant-tcc outputs.

    Returns:
        str: Path to the written abundance.tsv file.
    """
    tpm_mtx_path   = os.path.join(sample_out_dir, "matrix.abundance.tpm.mtx")
    transcripts_path = os.path.join(sample_out_dir, "transcripts.txt")
    output_tsv_path  = os.path.join(sample_out_dir, "abundance.tsv")

    # --- Validate inputs ---
    if not os.path.isfile(tpm_mtx_path):
        raise FileNotFoundError(f"TPM matrix not found: {tpm_mtx_path}")
    if not os.path.isfile(transcripts_path):
        raise FileNotFoundError(f"Transcripts file not found: {transcripts_path}")

    print(f"  Reading TPM matrix : {tpm_mtx_path}")
    print(f"  Reading transcripts: {transcripts_path}")

    # --- Load sparse TPM matrix and transcript labels ---
    tpm_matrix   = mmread(tpm_mtx_path)
    transcripts  = pd.read_csv(transcripts_path, header=None, sep="\t")

    # Convert sparse matrix to dense column vector, then attach transcript IDs
    tpm_values = np.asarray(tpm_matrix.todense()).flatten()
    df = pd.DataFrame({
        "transcript_id": transcripts.iloc[:, 0].values,
        "tpm":           tpm_values,
    })

    # Drop zero-expression transcripts
    df = df[df["tpm"] > 0]

    # Write output TSV
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
