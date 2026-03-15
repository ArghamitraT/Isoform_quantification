"""
compare_small_em.py
===================
Step 4 of debug plan: compare JK (JOLI EM) vs LK (kallisto quant-tcc) outputs
on the small test case, side by side with full logging.

Also performs Step 5: back-calculates the effective length that kallisto used
by dividing est_counts/total_reads (= theta) by tpm/1e6, then comparing to 1.0.

Inputs (set in CONFIG):
  JOLI_ABUNDANCE   : joli_abundance.tsv from debug_joli_em.py
  KALLISTO_ABUNDANCE: abundance.tsv from kallisto quant-tcc on the small test
                      (produced by running run_lr_kallisto.sh or manually)
  KALLISTO_MTX_DIR  : directory containing matrix.abundance.mtx and
                      matrix.abundance.tpm.mtx from kallisto quant-tcc
  TOTAL_READS       : n_pseudoaligned from run_info.json (for theta recovery)

Outputs:
  compare_log.txt   — full side-by-side comparison and eff_len analysis
"""

import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.io import mmread

# ============================================================
# CONFIG
# ============================================================
SMALL_TEST_DIR     = os.path.join(os.path.dirname(os.path.abspath(__file__)), "small_test")
JOLI_ABUNDANCE     = os.path.join(SMALL_TEST_DIR, "joli_abundance.tsv")
KALLISTO_MTX_DIR   = os.path.join(SMALL_TEST_DIR, "kallisto_out")
KALLISTO_ABUNDANCE = os.path.join(KALLISTO_MTX_DIR, "abundance.tsv")  # from mtx_to_tsv.py
TOTAL_READS        = None   # set to None to auto-load from run_info.json
# ============================================================

def log(msg, fh=None):
    print(msg)
    if fh:
        fh.write(msg + "\n")
        fh.flush()


def load_total_reads():
    run_info_path = os.path.join(SMALL_TEST_DIR, "run_info.json")
    with open(run_info_path) as f:
        info = json.load(f)
    return int(info.get("n_pseudoaligned", 0))


def recover_kallisto_theta_from_mtx(mtx_dir, transcripts_path, total_reads):
    """
    Recover per-transcript theta from kallisto's matrix.abundance.mtx
    (estimated counts matrix) by dividing by total_reads.

    Returns: dict {transcript_name: theta}
    """
    est_mtx_path = os.path.join(mtx_dir, "matrix.abundance.mtx")
    if not os.path.exists(est_mtx_path):
        return None

    tx_names = []
    with open(transcripts_path) as f:
        tx_names = [l.strip() for l in f if l.strip()]

    sparse = mmread(est_mtx_path).tocsr()
    est_counts = np.asarray(sparse.todense()).flatten().astype(np.float64)
    theta = est_counts / total_reads if total_reads > 0 else est_counts

    return {name: t for name, t in zip(tx_names, theta)}


def main():
    os.makedirs(KALLISTO_MTX_DIR, exist_ok=True)
    log_path = os.path.join(SMALL_TEST_DIR, "compare_log.txt")
    log_fh   = open(log_path, "w")

    log("=" * 70, log_fh)
    log("compare_small_em.py — JK vs LK EM comparison on small test case", log_fh)
    log("=" * 70, log_fh)

    total_reads = TOTAL_READS or load_total_reads()
    log(f"\ntotal_reads (n_pseudoaligned) = {total_reads}", log_fh)

    # ---- Load JK abundance ----
    if not os.path.exists(JOLI_ABUNDANCE):
        log(f"\nERROR: JOLI abundance not found: {JOLI_ABUNDANCE}", log_fh)
        log("  Run debug_joli_em.py first.", log_fh)
        sys.exit(1)

    jk_df = pd.read_csv(JOLI_ABUNDANCE, sep="\t")
    log(f"\nJK abundance loaded: {JOLI_ABUNDANCE}", log_fh)
    log(f"  Columns: {list(jk_df.columns)}", log_fh)
    log(f"  Rows: {len(jk_df)}", log_fh)

    # Normalise column names
    jk_id_col  = "target_id"  if "target_id"  in jk_df.columns else jk_df.columns[0]
    jk_tpm_col = "tpm"        if "tpm"        in jk_df.columns else "tpm"
    jk_ec_col  = "est_counts" if "est_counts" in jk_df.columns else None

    jk_map_tpm   = dict(zip(jk_df[jk_id_col], jk_df[jk_tpm_col]))
    jk_map_theta = {}
    if jk_ec_col:
        jk_map_theta = {row[jk_id_col]: row[jk_ec_col] / total_reads
                        for _, row in jk_df.iterrows()}

    # ---- Load LK abundance ----
    if not os.path.exists(KALLISTO_ABUNDANCE):
        log(f"\nWARNING: Kallisto abundance not found: {KALLISTO_ABUNDANCE}", log_fh)
        log("  Run kallisto quant-tcc on the small test case first, then", log_fh)
        log("  run mtx_to_tsv.py on the output directory.", log_fh)
        log("  Command:", log_fh)
        kallisto = "/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/kallisto/build/src/kallisto"
        idx      = "/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/SOTA/lr-kallisto/new_index.idx"
        log(f"    {kallisto} quant-tcc \\", log_fh)
        log(f"        --long -P PacBio -t 4 \\", log_fh)
        log(f"        {SMALL_TEST_DIR}/count.mtx \\", log_fh)
        log(f"        -i {idx} \\", log_fh)
        log(f"        -e {SMALL_TEST_DIR}/count.ec.txt \\", log_fh)
        log(f"        -o {KALLISTO_MTX_DIR}/", log_fh)
        log(f"    python mtx_to_tsv.py {KALLISTO_MTX_DIR}", log_fh)
        lk_available = False
    else:
        lk_df = pd.read_csv(KALLISTO_ABUNDANCE, sep="\t")
        log(f"\nLK abundance loaded: {KALLISTO_ABUNDANCE}", log_fh)
        log(f"  Columns: {list(lk_df.columns)}", log_fh)
        log(f"  Rows: {len(lk_df)}", log_fh)
        lk_id_col  = "transcript_id" if "transcript_id" in lk_df.columns else lk_df.columns[0]
        lk_tpm_col = "tpm"           if "tpm"           in lk_df.columns else lk_df.columns[1]
        lk_map_tpm = dict(zip(lk_df[lk_id_col], lk_df[lk_tpm_col]))
        lk_available = True

    # ---- Load kallisto theta from matrix.abundance.mtx (if available) ----
    tx_path = os.path.join(SMALL_TEST_DIR, "transcripts.txt")
    lk_theta = recover_kallisto_theta_from_mtx(KALLISTO_MTX_DIR, tx_path, total_reads)
    if lk_theta:
        log(f"\nKallisto theta recovered from matrix.abundance.mtx:", log_fh)
        for name, t in sorted(lk_theta.items(), key=lambda x: -x[1]):
            log(f"  {name:<35} theta={t:.8f}", log_fh)

    # ---- Side-by-side TPM comparison ----
    if lk_available:
        all_tx = sorted(set(jk_map_tpm) | set(lk_map_tpm))
        log(f"\n{'='*70}", log_fh)
        log(f"SIDE-BY-SIDE COMPARISON (all {len(all_tx)} transcripts in small test)", log_fh)
        log(f"{'='*70}", log_fh)
        header = f"{'transcript':<35}  {'JK_tpm':>12}  {'LK_tpm':>12}  {'abs_diff':>12}  {'rel_diff_%':>12}"
        log(header, log_fh)
        log("-" * len(header), log_fh)

        diffs = []
        for tx in all_tx:
            jk_tpm = jk_map_tpm.get(tx, 0.0)
            lk_tpm = lk_map_tpm.get(tx, 0.0)
            abs_d  = abs(jk_tpm - lk_tpm)
            ref    = max(jk_tpm, lk_tpm)
            rel_d  = 100.0 * abs_d / ref if ref > 0 else 0.0
            diffs.append(abs_d)
            log(f"  {tx:<35}  {jk_tpm:>12.4f}  {lk_tpm:>12.4f}  {abs_d:>12.4f}  {rel_d:>11.1f}%", log_fh)

        log(f"\nSummary:", log_fh)
        log(f"  Max absolute TPM diff : {max(diffs):.4f}", log_fh)
        log(f"  Mean absolute TPM diff: {np.mean(diffs):.4f}", log_fh)
        log(f"  Transcripts diff > 1% : {sum(1 for d, tx in zip(diffs, all_tx) if d > 0.01 * max(jk_map_tpm.get(tx,0), lk_map_tpm.get(tx,0)))}", log_fh)

    # ---- Step 5: back-calculate effective lengths from kallisto ----
    if lk_theta and lk_available:
        log(f"\n{'='*70}", log_fh)
        log("STEP 5: Back-calculate effective lengths used by kallisto", log_fh)
        log("  Formula: eff_len = (theta / tpm) * 1e6 * sum(theta/eff_len)", log_fh)
        log("  Simplified: if tpm = theta/eff_len / sum(theta/eff_len) * 1e6", log_fh)
        log("  Then: theta/tpm is proportional to eff_len", log_fh)
        log("  So: eff_len_ratio[t] = (theta[t]/tpm[t]) / mean(theta/tpm)", log_fh)
        log(f"{'='*70}", log_fh)

        ratios = {}
        for tx in lk_theta:
            tpm = lk_map_tpm.get(tx, 0.0)
            t   = lk_theta[tx]
            if tpm > 0 and t > 0:
                ratios[tx] = t / (tpm / 1e6)

        if ratios:
            mean_ratio = np.mean(list(ratios.values()))
            log(f"\n  theta / (tpm/1e6)  [proportional to eff_len]:", log_fh)
            for tx, r in sorted(ratios.items()):
                norm = r / mean_ratio
                log(f"  {tx:<35}  ratio={r:.4f}  normalised={norm:.4f}  "
                    f"{'~= 1.0 (uniform)' if abs(norm - 1.0) < 0.01 else 'DIFFERS from 1.0'}", log_fh)

            all_same = all(abs(r/mean_ratio - 1.0) < 0.01 for r in ratios.values())
            log(f"\n  Conclusion: kallisto eff_lens are "
                f"{'UNIFORM (all ~1.0)' if all_same else 'NON-UNIFORM (real lengths used)'}",
                log_fh)

    log(f"\nLog written to: {log_path}", log_fh)
    log_fh.close()
    print(f"\nLog: {log_path}")


if __name__ == "__main__":
    main()
