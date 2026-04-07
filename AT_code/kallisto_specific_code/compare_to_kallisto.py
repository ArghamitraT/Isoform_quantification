import os
import numpy as np
import pandas as pd
from kallisto_ec_loader import read_count_mtx

KALLISTO_DIR = "/gpfs/commons/home/tmehta/proj/kallisto_lr/ds_52_sub_2000"
OUT_DIR = "/gpfs/commons/home/tmehta/proj/JOLI/Results/ec_debug"

# Try to load raw counts (est_counts) instead of TPM
count_mtx_path = os.path.join(KALLISTO_DIR, "matrix.abundance.mtx")  # raw counts
tpm_mtx_path = os.path.join(KALLISTO_DIR, "matrix.abundance.tpm.mtx")
tx_path = os.path.join(KALLISTO_DIR, "transcripts.txt")
len_path = os.path.join(KALLISTO_DIR, "transcript_lengths.txt")

# Load transcripts and lengths
with open(tx_path) as f:
    k_tx_names = [ln.strip() for ln in f if ln.strip()]
with open(len_path) as f:
    tx_lengths = np.array([float(ln.strip()) for ln in f if ln.strip()])

len_df = pd.DataFrame({"transcript": k_tx_names, "length": tx_lengths})

# Load Kallisto raw counts
if os.path.exists(count_mtx_path):
    nrows, ncols, nnz, rows, cols, vals = read_count_mtx(count_mtx_path)
    k_counts = np.zeros(ncols, dtype=float)
    for r, c, v in zip(rows, cols, vals):
        k_counts[c - 1] += v
    k_df = pd.DataFrame({"transcript": k_tx_names, "est_count": k_counts})
    k_df["theta_kallisto_count"] = k_df["est_count"] / k_df["est_count"].sum()
    print(f"Loaded raw counts from {count_mtx_path}")
else:
    print(f"Raw counts file not found: {count_mtx_path}")
    k_df = None

# Load Kallisto TPM
nrows, ncols, nnz, rows, cols, vals = read_count_mtx(tpm_mtx_path)
k_tpm = np.zeros(ncols, dtype=float)
for r, c, v in zip(rows, cols, vals):
    k_tpm[c - 1] += v

tpm_df = pd.DataFrame({"transcript": k_tx_names, "tpm": k_tpm})
tpm_df["theta_kallisto_tpm"] = tpm_df["tpm"] / tpm_df["tpm"].sum()

# Load JOLI output
theta_path = os.path.join(OUT_DIR, "theta_joli.tsv")
joli_df = pd.read_csv(theta_path, sep="\t")
joli_df = joli_df.merge(len_df, on="transcript", how="left")

# Merge all data
merged = joli_df.merge(tpm_df[["transcript", "theta_kallisto_tpm"]], on="transcript", how="inner")
if k_df is not None:
    merged = merged.merge(k_df[["transcript", "theta_kallisto_count"]], on="transcript", how="left")

# Compute JOLI TPM (theta / length, normalized)
merged["joli_rho"] = merged["theta_joli"] / merged["length"]
merged["joli_tpm"] = merged["joli_rho"] / merged["joli_rho"].sum()

print(f"\n=== Comparison Results ===")
print(f"Common transcripts: {len(merged)}")

print(f"\n--- JOLI theta vs Kallisto TPM ---")
print(f"Pearson:  {merged['theta_joli'].corr(merged['theta_kallisto_tpm'], method='pearson'):.4f}")
print(f"Spearman: {merged['theta_joli'].corr(merged['theta_kallisto_tpm'], method='spearman'):.4f}")

print(f"\n--- JOLI TPM vs Kallisto TPM ---")
print(f"Pearson:  {merged['joli_tpm'].corr(merged['theta_kallisto_tpm'], method='pearson'):.4f}")
print(f"Spearman: {merged['joli_tpm'].corr(merged['theta_kallisto_tpm'], method='spearman'):.4f}")

if k_df is not None:
    print(f"\n--- JOLI theta vs Kallisto raw counts (normalized) ---")
    print(f"Pearson:  {merged['theta_joli'].corr(merged['theta_kallisto_count'], method='pearson'):.4f}")
    print(f"Spearman: {merged['theta_joli'].corr(merged['theta_kallisto_count'], method='spearman'):.4f}")

# Show top discrepancies
print(f"\n=== Top Discrepancies ===")
merged["diff"] = abs(merged["theta_joli"] - merged["theta_kallisto_tpm"])
top_diff = merged.nlargest(10, "diff")
print(top_diff[["transcript", "length", "theta_joli", "theta_kallisto_tpm", "n_joli"]].to_string())