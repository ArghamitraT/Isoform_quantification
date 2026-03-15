"""
debug_joli_em.py
================
Step 2 of debug plan: run JK EM on the small test case with full verbose
logging — prints theta at EVERY round and per-EC assignments in round 1.

This is a debug-only script; it does NOT modify any production code.
It imports the production modules (load_tcc, weights, em_algorithm,
output_writer) directly and wraps them with extra instrumentation.

Inputs (set in CONFIG):
  SMALL_TEST_DIR : folder produced by extract_small_case.py
  OUT_DIR        : where to write joli_abundance.tsv and the debug log

Outputs:
  joli_abundance.tsv  — final abundance in kallisto 5-column format
  debug_em_log.txt    — full per-round theta trace + round-1 EC breakdown
"""

import os
import sys
import time

import numpy as np

# Add JOLI_Kallisto to path so we can import production modules
JOLI_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, JOLI_DIR)

from load_tcc   import load_tcc_data
from weights    import compute_weights
from output_writer import write_abundance

# ============================================================
# CONFIG
# ============================================================
SMALL_TEST_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "small_test")
OUT_DIR        = SMALL_TEST_DIR   # write output alongside small test files
MAX_EM_ROUNDS  = 10000
MIN_ROUNDS     = 50
EFF_LEN_MODE   = "uniform"        # must match JK production setting
# ============================================================

TOLERANCE          = np.finfo(np.float64).tiny
ALPHA_LIMIT        = 1e-7
ALPHA_CHANGE_LIMIT = 1e-2
ALPHA_CHANGE       = 1e-2


def log(msg, fh=None):
    print(msg)
    if fh:
        fh.write(msg + "\n")
        fh.flush()


def run_verbose_em(tcc_data, weight_data, log_fh):
    """
    Run EM on tcc_data/weight_data with per-round logging.

    Mirrors JoliEM.run() exactly but prints theta after every round.
    Also prints full per-EC breakdown for round 1 to show fractional
    read assignments.

    Returns:
        theta (np.ndarray) — final normalized abundances
        n_rounds (int)
        converged (bool)
    """
    n_tx         = len(tcc_data.transcript_names)
    ec_counts    = tcc_data.ec_counts
    ec_txs       = tcc_data.ec_transcripts
    ec_weights   = weight_data.ec_weights
    n_ecs        = len(ec_counts)
    tx_names     = tcc_data.transcript_names

    log(f"\n{'='*60}", log_fh)
    log(f"JoliEM verbose run", log_fh)
    log(f"  n_transcripts : {n_tx}", log_fh)
    log(f"  n_ecs         : {n_ecs}", log_fh)
    log(f"  total_reads   : {tcc_data.total_reads}", log_fh)
    log(f"  MAX_EM_ROUNDS : {MAX_EM_ROUNDS}", log_fh)
    log(f"  MIN_ROUNDS    : {MIN_ROUNDS}", log_fh)
    log(f"  eff_len_mode  : {EFF_LEN_MODE}", log_fh)
    log(f"  eff_lens      : {weight_data.eff_lens.tolist()}", log_fh)
    log(f"{'='*60}", log_fh)

    # ---- Initialise theta (uniform, same as kallisto) ----
    theta = np.full(n_tx, 1.0 / n_tx, dtype=np.float64)
    log(f"\nInitial theta (uniform 1/{n_tx} = {1.0/n_tx:.6e}):", log_fh)
    for i, (name, t) in enumerate(zip(tx_names, theta)):
        log(f"  tx[{i:>3}] {name:<30}  theta={t:.6e}", log_fh)

    converged = False

    for round_num in range(MAX_EM_ROUNDS):
        n = np.zeros(n_tx, dtype=np.float64)

        is_round1 = (round_num == 0)
        if is_round1:
            log(f"\n{'='*60}", log_fh)
            log(f"ROUND 1 — detailed EC breakdown:", log_fh)
            log(f"{'='*60}", log_fh)

        # ---- E + M step ----
        for ec_id, (txs, w_arr, cnt) in enumerate(
                zip(ec_txs, ec_weights, ec_counts)):

            if cnt == 0:
                continue

            if len(txs) == 1:
                # Single-tx EC: deterministic assignment
                n[txs[0]] += cnt
                if is_round1:
                    log(f"\n  EC {ec_id:>4}  count={cnt:>6.0f}  "
                        f"single-tx → all to tx[{txs[0]}] "
                        f"({tx_names[txs[0]]})", log_fh)
            else:
                # Multi-tx EC: fractional assignment
                theta_w  = np.array([theta[t] * w_arr[i]
                                     for i, t in enumerate(txs)])
                denom    = theta_w.sum()

                if is_round1:
                    log(f"\n  EC {ec_id:>4}  count={cnt:>6.0f}  "
                        f"multi-tx ({len(txs)} transcripts):", log_fh)
                    log(f"    denom = {denom:.6e}", log_fh)

                if denom >= TOLERANCE:
                    fracs = theta_w / denom
                    for i, t in enumerate(txs):
                        contrib = cnt * fracs[i]
                        n[t]   += contrib
                        if is_round1:
                            log(f"    tx[{t:>3}] {tx_names[t]:<30}  "
                                f"theta={theta[t]:.4e}  w={w_arr[i]:.4e}  "
                                f"frac={fracs[i]:.4f}  contrib={contrib:.4f}",
                                log_fh)
                else:
                    if is_round1:
                        log(f"    denom below TOLERANCE — skipping", log_fh)

        # ---- Normalise ----
        total = n.sum()
        theta_new = n / total if total > 0 else theta.copy()

        # ---- Convergence check ----
        changed = int(np.sum(
            (theta_new > ALPHA_CHANGE_LIMIT) &
            (np.abs(theta_new - theta) /
             np.maximum(theta_new, TOLERANCE) > ALPHA_CHANGE)
        ))

        # ---- Per-round log ----
        nonzero = int((theta_new > 0).sum())
        log(f"\nRound {round_num+1:>5}  |  changed={changed:>5}  "
            f"nonzero_tx={nonzero}  total_n={total:.1f}", log_fh)

        # Print full theta every round (small test case has few transcripts)
        for i, (name, t_old, t_new) in enumerate(
                zip(tx_names, theta, theta_new)):
            delta = t_new - t_old
            log(f"  tx[{i:>3}] {name:<30}  "
                f"theta_old={t_old:.6e}  theta_new={t_new:.6e}  "
                f"delta={delta:+.2e}", log_fh)

        theta = theta_new

        if changed == 0 and round_num >= MIN_ROUNDS:
            converged = True
            log(f"\n>>> Converged at round {round_num+1}", log_fh)
            break

    if not converged:
        log(f"\n>>> Reached max_em_rounds={MAX_EM_ROUNDS} without convergence.",
            log_fh)

    # ---- Zero out very small abundances ----
    n_zeroed = int(((theta > 0) & (theta < ALPHA_LIMIT / 10)).sum())
    theta[theta < ALPHA_LIMIT / 10] = 0.0
    if n_zeroed:
        log(f"\nZeroed {n_zeroed} transcripts below {ALPHA_LIMIT/10:.1e}", log_fh)

    log(f"\n{'='*60}", log_fh)
    log(f"FINAL THETA after EM:", log_fh)
    for i, (name, t) in enumerate(zip(tx_names, theta)):
        log(f"  tx[{i:>3}] {name:<30}  theta={t:.8f}", log_fh)

    return theta, round_num + 1, converged


def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    log_path = os.path.join(OUT_DIR, "debug_em_log.txt")
    log_fh   = open(log_path, "w")

    log("=" * 60, log_fh)
    log("debug_joli_em.py", log_fh)
    log(f"small_test_dir : {SMALL_TEST_DIR}", log_fh)
    log(f"out_dir        : {OUT_DIR}", log_fh)
    log("=" * 60, log_fh)

    # Step 1: load TCC data
    log("\n--- Loading TCC data ---", log_fh)
    tcc_data = load_tcc_data(SMALL_TEST_DIR)
    log(f"  ECs loaded          : {len(tcc_data.ec_counts)}", log_fh)
    log(f"  Transcripts loaded  : {len(tcc_data.transcript_names)}", log_fh)
    log(f"  Total reads         : {tcc_data.total_reads}", log_fh)
    log(f"  ECs with count > 0  : {int((tcc_data.ec_counts > 0).sum())}", log_fh)

    log("\n  EC counts (all ECs in small test):", log_fh)
    for i, cnt in enumerate(tcc_data.ec_counts):
        txs   = tcc_data.ec_transcripts[i]
        names = [tcc_data.transcript_names[t] for t in txs]
        log(f"    EC {i:>3}  count={cnt:>6.0f}  txs={txs}  names={names}", log_fh)

    # Step 2: compute weights
    log("\n--- Computing weights ---", log_fh)
    weight_data = compute_weights(tcc_data, mode=EFF_LEN_MODE)
    log(f"  eff_lens : {weight_data.eff_lens.tolist()}", log_fh)
    log(f"  ec_weights (per EC):", log_fh)
    for i, w in enumerate(weight_data.ec_weights):
        log(f"    EC {i:>3}  weights={w.tolist()}", log_fh)

    # Step 3: run verbose EM
    theta, n_rounds, converged = run_verbose_em(tcc_data, weight_data, log_fh)

    # Step 4: write abundance.tsv
    out_path = os.path.join(OUT_DIR, "joli_abundance.tsv")
    log(f"\n--- Writing abundance.tsv ---", log_fh)
    summary = write_abundance(
        theta            = theta,
        eff_lens         = weight_data.eff_lens,
        total_reads      = tcc_data.total_reads,
        transcript_names = tcc_data.transcript_names,
        output_path      = out_path,
    )
    log(f"  n_transcripts  : {summary['n_transcripts']}", log_fh)
    log(f"  n_nonzero      : {summary['n_nonzero']}", log_fh)
    log(f"  total_tpm      : {summary['total_tpm']:.2f}", log_fh)
    log(f"  total_counts   : {summary['total_est_counts']:.2f}", log_fh)

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.2f}s", log_fh)
    log(f"Log written to: {log_path}", log_fh)
    log(f"Abundance written to: {out_path}", log_fh)
    log_fh.close()

    print(f"\nLog    : {log_path}")
    print(f"Output : {out_path}")


if __name__ == "__main__":
    main()
