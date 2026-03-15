"""
create_micro_test.py
====================
Creates a micro test case with exactly 2 transcripts and 3 ECs — the same
structure as the toy_em_example.py, but using REAL transcript names and
global indices so both JK and kallisto quant-tcc can actually run on it.

Structure (mirrors toy_em_example.py exactly):
  EC 0: only transcript A  (single-tx)  — 30 reads
  EC 1: only transcript B  (single-tx)  — 10 reads
  EC 2: both A and B       (multi-tx)   — 60 reads

Transcripts chosen from the small_test (real entries in the 218K kallisto index):
  A = NM_001402.6  (compact index 1 in small_test → global index 68660)
  B = NM_001101.5  (compact index 2 in small_test → global index 73705)

Two output files differ intentionally:
  matrix.ec      — uses compact indices 0,1 (for JK, which only knows these 2 tx)
  count.ec.txt   — uses global indices 68660,73705 (for kallisto, which has all 218K tx)

Outputs (in debug/micro_test/):
  count.mtx         — 1×3 sparse matrix (3 ECs, counts 30/10/60)
  matrix.ec         — EC definitions with compact tx indices (for JK)
  count.ec.txt      — EC definitions with global tx indices (for kallisto)
  transcripts.txt   — 2 transcript names
  run_info.json     — minimal run_info with n_pseudoaligned=100
"""

import json
import os

import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

# ============================================================
# CONFIG
# ============================================================
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "micro_test")

# Transcripts
TX_NAMES   = ["NM_001402.6", "NM_001101.5"]   # compact indices 0, 1
TX_GLOBALS = [68660, 73705]                    # global kallisto index positions

# ECs: (compact_tx_indices, count)
# Mirrors toy_em_example.py exactly
ECS = [
    ([0],    30),   # EC 0: only tx A (single-tx)
    ([1],    10),   # EC 1: only tx B (single-tx)
    ([0, 1], 60),   # EC 2: tx A and B (multi-tx)
]
# ============================================================


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f"Writing micro test to: {OUT_DIR}")

    n_ecs = len(ECS)

    # ---- count.mtx ----
    row  = np.zeros(n_ecs, dtype=np.int32)
    col  = np.arange(n_ecs, dtype=np.int32)
    data = np.array([cnt for _, cnt in ECS], dtype=np.float64)
    mat  = csr_matrix((data, (row, col)), shape=(1, n_ecs))
    mmwrite(os.path.join(OUT_DIR, "count.mtx"), mat, comment="", field="real")
    print(f"  count.mtx  : {n_ecs} ECs, counts={data.tolist()}")

    # ---- matrix.ec (compact indices — for JK) ----
    with open(os.path.join(OUT_DIR, "matrix.ec"), "w") as f:
        for ec_id, (txs, _) in enumerate(ECS):
            f.write(f"{ec_id}\t{','.join(str(t) for t in txs)}\n")
    print(f"  matrix.ec  : compact indices (for JK)")

    # ---- count.ec.txt (global indices — for kallisto) ----
    with open(os.path.join(OUT_DIR, "count.ec.txt"), "w") as f:
        for ec_id, (txs, _) in enumerate(ECS):
            global_txs = [TX_GLOBALS[t] for t in txs]
            f.write(f"{ec_id}\t{','.join(str(t) for t in global_txs)}\n")
    print(f"  count.ec.txt: global indices (for kallisto): {TX_GLOBALS}")

    # ---- transcripts.txt ----
    with open(os.path.join(OUT_DIR, "transcripts.txt"), "w") as f:
        for name in TX_NAMES:
            f.write(name + "\n")
    print(f"  transcripts.txt: {TX_NAMES}")

    # ---- run_info.json ----
    run_info = {"n_pseudoaligned": sum(cnt for _, cnt in ECS)}
    with open(os.path.join(OUT_DIR, "run_info.json"), "w") as f:
        json.dump(run_info, f, indent=2)
    print(f"  run_info.json: n_pseudoaligned={run_info['n_pseudoaligned']}")

    # ---- Summary ----
    print()
    print("=" * 50)
    print("MICRO TEST STRUCTURE")
    print("=" * 50)
    print(f"{'EC':>4}  {'type':<12}  {'count':>6}  compact_txs  global_txs")
    print("-" * 60)
    for ec_id, (txs, cnt) in enumerate(ECS):
        kind       = "single-tx" if len(txs) == 1 else "multi-tx"
        globals_   = [TX_GLOBALS[t] for t in txs]
        names      = [TX_NAMES[t] for t in txs]
        print(f"  {ec_id:>2}  {kind:<12}  {cnt:>6}  {txs}  {globals_}  {names}")
    print()
    print("Run JK on this with debug_joli_em.py (change SMALL_TEST_DIR to micro_test/)")
    print()
    print("Run kallisto:")
    print("  KALLISTO=/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/kallisto/build/src/kallisto")
    print(f"  $KALLISTO quant-tcc --long -P PacBio -t 4 \\")
    print(f"      {OUT_DIR}/count.mtx \\")
    print(f"      -i /gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/new_index.idx \\")
    print(f"      -e {OUT_DIR}/count.ec.txt \\")
    print(f"      -o {OUT_DIR}/kallisto_out/")


if __name__ == "__main__":
    main()
