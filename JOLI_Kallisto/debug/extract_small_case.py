"""
extract_small_case.py
=====================
Step 1 of debug plan: extract a small, self-contained test case from a full
bustools count output, for side-by-side EM debugging between JK and LK.

What it does:
  - Reads count.mtx + matrix.ec from a real run output folder
  - Picks the top N ECs by read count (default: 20)
  - Collects all transcripts that appear in those ECs
  - Re-indexes everything (new compact EC ids, new compact transcript ids)
  - Writes a minimal count.mtx, matrix.ec, transcripts.txt, run_info.json
    to the small_test/ output folder

Inputs (set in CONFIG):
  SOURCE_DIR   : folder with count.mtx, matrix.ec, transcripts.txt, run_info.json
  OUT_DIR      : where to write the small test case files
  TOP_N_ECS    : how many ECs to keep (by descending read count)

Outputs (in OUT_DIR):
  count.mtx        — MatrixMarket sparse matrix (1 row x TOP_N_ECS cols)
  matrix.ec        — EC-to-transcript mapping (re-indexed)
  transcripts.txt  — transcript names for the compact set
  run_info.json    — copy of original run_info (unchanged)
  extraction_log.txt — detailed log of what was extracted
"""

import json
import os
import sys

import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix

# ============================================================
# CONFIG
# ============================================================
SOURCE_DIR = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_2026_03_14__10_42_22/ds_52_furtherDownsampled"
OUT_DIR    = os.path.join(os.path.dirname(os.path.abspath(__file__)), "small_test")
TOP_N_ECS  = 20
# ============================================================

def log(msg, fh=None):
    print(msg)
    if fh:
        fh.write(msg + "\n")


def load_ec_map(ec_path):
    """Read matrix.ec -> list of (ec_idx, [tx_idx, ...])."""
    ecs = []
    with open(ec_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            idx, txs = line.split("\t")
            tx_list = [int(t) for t in txs.split(",")]
            ecs.append((int(idx), tx_list))
    return ecs  # list[(ec_id, [tx_ids])]


def load_counts(mtx_path, n_ecs):
    """Read count.mtx -> np.ndarray shape (n_ecs,) of counts."""
    sparse = mmread(mtx_path).tocsr()
    counts = np.asarray(sparse.todense()).flatten().astype(np.float64)
    # Pad/trim to n_ecs
    if len(counts) < n_ecs:
        counts = np.pad(counts, (0, n_ecs - len(counts)))
    elif len(counts) > n_ecs:
        counts = counts[:n_ecs]
    return counts


def load_transcripts(tx_path):
    """Read transcripts.txt -> list of transcript names."""
    with open(tx_path) as f:
        return [line.strip() for line in f if line.strip()]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    log_path = os.path.join(OUT_DIR, "extraction_log.txt")
    log_fh   = open(log_path, "w")

    log("=" * 60, log_fh)
    log("extract_small_case.py", log_fh)
    log(f"Source : {SOURCE_DIR}", log_fh)
    log(f"Out    : {OUT_DIR}", log_fh)
    log(f"Top N  : {TOP_N_ECS}", log_fh)
    log("=" * 60, log_fh)

    # ---- Load source data ----
    ec_list      = load_ec_map(os.path.join(SOURCE_DIR, "matrix.ec"))
    transcripts  = load_transcripts(os.path.join(SOURCE_DIR, "transcripts.txt"))
    counts       = load_counts(os.path.join(SOURCE_DIR, "count.mtx"), len(ec_list))

    log(f"\nLoaded {len(ec_list)} ECs, {len(transcripts)} transcripts", log_fh)
    log(f"Total reads in count.mtx: {counts.sum():.0f}", log_fh)

    # ---- Pick top N ECs by count ----
    # Sort EC indices by descending count
    sorted_by_count = sorted(range(len(counts)), key=lambda i: counts[i], reverse=True)
    top_ec_orig_ids = sorted_by_count[:TOP_N_ECS]

    log(f"\nTop {TOP_N_ECS} ECs selected:", log_fh)
    log(f"  {'orig_ec_id':>12}  {'count':>8}  {'n_transcripts':>14}  transcript_indices", log_fh)
    for orig_id in top_ec_orig_ids:
        tx_list = ec_list[orig_id][1]
        log(f"  {orig_id:>12}  {counts[orig_id]:>8.0f}  {len(tx_list):>14}  {tx_list}", log_fh)

    # ---- Collect all unique transcript indices from selected ECs ----
    used_tx_orig = sorted(set(
        tx for orig_id in top_ec_orig_ids
        for tx in ec_list[orig_id][1]
    ))
    log(f"\nUnique transcripts involved: {len(used_tx_orig)}", log_fh)
    log(f"  Original tx indices: {used_tx_orig}", log_fh)
    log(f"  Transcript names:", log_fh)
    for oi in used_tx_orig:
        log(f"    [{oi}] {transcripts[oi]}", log_fh)

    # ---- Build compact re-indexing ----
    # new tx index 0..K-1  →  original tx index
    orig_to_new_tx = {orig: new for new, orig in enumerate(used_tx_orig)}
    new_to_name    = [transcripts[oi] for oi in used_tx_orig]

    # new EC index 0..N-1  →  original EC index
    new_ec_list    = []
    new_counts     = []
    for new_ec_id, orig_ec_id in enumerate(top_ec_orig_ids):
        orig_txs  = ec_list[orig_ec_id][1]
        new_txs   = [orig_to_new_tx[t] for t in orig_txs]
        new_ec_list.append(new_txs)
        new_counts.append(counts[orig_ec_id])

    log(f"\nRe-indexed: {len(new_ec_list)} ECs, {len(new_to_name)} transcripts", log_fh)
    log("\nNew EC definitions:", log_fh)
    for i, (txs, cnt) in enumerate(zip(new_ec_list, new_counts)):
        names = [new_to_name[t] for t in txs]
        log(f"  EC {i:>3}  count={cnt:>6.0f}  txs={txs}  names={names}", log_fh)

    # ---- Write count.mtx ----
    # MatrixMarket: 1 row (sample) x n_ecs cols; each EC gets its count
    n_ecs = len(new_ec_list)
    row   = np.zeros(n_ecs, dtype=np.int32)   # all in row 0
    col   = np.arange(n_ecs, dtype=np.int32)
    data  = np.array(new_counts, dtype=np.float64)
    sparse_mat = csr_matrix((data, (row, col)), shape=(1, n_ecs))

    mtx_path = os.path.join(OUT_DIR, "count.mtx")
    mmwrite(mtx_path, sparse_mat, comment="", field="real")
    log(f"\nWrote count.mtx: {mtx_path}", log_fh)

    # ---- Write matrix.ec ----
    ec_path = os.path.join(OUT_DIR, "matrix.ec")
    with open(ec_path, "w") as f:
        for i, txs in enumerate(new_ec_list):
            f.write(f"{i}\t{','.join(str(t) for t in txs)}\n")
    log(f"Wrote matrix.ec: {ec_path}", log_fh)

    # Write count.ec.txt with ORIGINAL global transcript indices
    # (kallisto quant-tcc reads this and maps against its full 218K index,
    #  so indices must be the original global ones, not the re-indexed compact ones)
    count_ec_path = os.path.join(OUT_DIR, "count.ec.txt")
    with open(count_ec_path, "w") as f:
        for new_ec_id, orig_ec_id in enumerate(top_ec_orig_ids):
            orig_txs = ec_list[orig_ec_id][1]   # original global indices
            f.write(f"{new_ec_id}\t{','.join(str(t) for t in orig_txs)}\n")
    log(f"Wrote count.ec.txt (original global tx indices, for kallisto): {count_ec_path}", log_fh)

    # ---- Write transcripts.txt ----
    tx_path = os.path.join(OUT_DIR, "transcripts.txt")
    with open(tx_path, "w") as f:
        for name in new_to_name:
            f.write(name + "\n")
    log(f"Wrote transcripts.txt: {tx_path}", log_fh)

    # ---- Copy run_info.json (unchanged — total_reads stays as original) ----
    run_info_src = os.path.join(SOURCE_DIR, "run_info.json")
    run_info_dst = os.path.join(OUT_DIR, "run_info.json")
    if os.path.exists(run_info_src):
        with open(run_info_src) as f:
            run_info = json.load(f)
        with open(run_info_dst, "w") as f:
            json.dump(run_info, f, indent=2)
        log(f"Copied run_info.json: {run_info_dst}", log_fh)
        log(f"  n_pseudoaligned = {run_info.get('n_pseudoaligned', 'N/A')}", log_fh)

    # ---- Summary ----
    total_small = sum(new_counts)
    log("\n" + "=" * 60, log_fh)
    log("SUMMARY", log_fh)
    log(f"  ECs in small test:         {n_ecs}", log_fh)
    log(f"  Transcripts in small test: {len(new_to_name)}", log_fh)
    log(f"  Total reads in small test: {total_small:.0f}", log_fh)
    log(f"  (out of {counts.sum():.0f} total in full dataset)", log_fh)
    log(f"  Coverage: {100*total_small/counts.sum():.1f}% of all reads", log_fh)
    log("=" * 60, log_fh)
    log(f"\nOutput files in: {OUT_DIR}", log_fh)
    log("  count.mtx, matrix.ec, count.ec.txt, transcripts.txt, run_info.json", log_fh)

    log_fh.close()
    print(f"\nLog written to: {log_path}")


if __name__ == "__main__":
    main()
