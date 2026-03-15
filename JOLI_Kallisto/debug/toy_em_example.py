"""
toy_em_example.py
=================
A minimal, fully-annotated EM example to build concrete intuition before
reading the real code.

Setup: 2 transcripts, 3 equivalence classes
  EC 0: only transcript A   (single-tx) — 30 reads
  EC 1: only transcript B   (single-tx) — 10 reads
  EC 2: both A and B        (multi-tx)  — 60 reads

True answer (if eff_lens are equal): theta_A should be higher because
EC 2's reads will be split ~proportional to theta, and EC 0 pushes A up.

We run TWO implementations side by side:
  - JK style: single-tx ECs inside loop, normalize every round, eff_len=1
  - LK style: single-tx ECs skipped in loop, added post-EM, eff_len varies

Inputs: none (hardcoded)
Outputs: printed trace of every EM round + final TPM comparison
"""

import numpy as np

# ============================================================
# CONFIG
# ============================================================
N_ROUNDS   = 200
MIN_ROUNDS = 50

# Equivalence classes
EC_COUNTS = [30, 10, 60]          # reads per EC
EC_TX     = [[0], [1], [0, 1]]    # transcripts per EC (0=A, 1=B)

# Effective lengths
# Try: EFF_LEN_UNIFORM = [1.0, 1.0]  → should match when single-tx handling also matches
# Try: EFF_LEN_REAL    = [2.0, 1.0]  → A is twice as long as B
EFF_LEN_JK = [1.0, 1.0]           # JK Phase 1: uniform
EFF_LEN_LK = [2.0, 1.0]           # LK: real lengths (A is longer)

TX_NAMES = ["transcript_A", "transcript_B"]

# Convergence constants (same in both implementations)
TOLERANCE          = np.finfo(np.float64).tiny
ALPHA_CHANGE_LIMIT = 1e-2
ALPHA_CHANGE       = 1e-2
# ============================================================


def convergence_check(theta_new, theta_old):
    """
    Count transcripts that changed significantly.
    Matches kallisto EMAlgorithm.h and JK em_algorithm.py exactly.

    Args:
        theta_new (np.ndarray): updated abundances.
        theta_old (np.ndarray): previous abundances.

    Returns:
        int: number of transcripts that changed significantly.
    """
    return int(np.sum(
        (theta_new > ALPHA_CHANGE_LIMIT) &
        (np.abs(theta_new - theta_old) /
         np.maximum(theta_new, TOLERANCE) > ALPHA_CHANGE)
    ))


# ===========================================================
# JK-style EM
# Single-tx ECs inside loop, normalize every round, eff_len=1
# ===========================================================

def run_jk_em(verbose=True):
    """
    JOLI-Kallisto style EM.

    Single-tx ECs contribute to every round's n[], which is then normalized.
    Effective lengths are uniform (1.0), so TPM = theta * 1e6 / sum(theta) = theta * 1e6.

    Args:
        verbose (bool): print per-round trace.

    Returns:
        np.ndarray: final theta (normalized probabilities).
    """
    n_tx    = 2
    eff_len = np.array(EFF_LEN_JK, dtype=np.float64)
    counts  = np.array(EC_COUNTS, dtype=np.float64)
    theta   = np.full(n_tx, 1.0 / n_tx)   # uniform init

    if verbose:
        print("=" * 60)
        print("JK-style EM (single-tx INSIDE loop, eff_len=1)")
        print("=" * 60)
        print(f"Init theta = {theta}")
        print()

    for rnd in range(N_ROUNDS):
        n = np.zeros(n_tx)

        for ec_id, (txs, cnt) in enumerate(zip(EC_TX, counts)):
            if cnt == 0:
                continue

            if len(txs) == 1:
                # Single-tx EC: deterministic assignment — included every round
                n[txs[0]] += cnt
                if verbose and rnd == 0:
                    print(f"  Round 1, EC {ec_id} (single-tx → tx{txs[0]}): "
                          f"cnt={cnt}, n[{txs[0]}] += {cnt}")
            else:
                # Multi-tx EC: fractional assignment
                weights   = np.array([1.0 / eff_len[t] for t in txs])
                theta_w   = np.array([theta[t] * weights[i]
                                      for i, t in enumerate(txs)])
                denom     = theta_w.sum()
                if verbose and rnd == 0:
                    print(f"  Round 1, EC {ec_id} (multi-tx {txs}): "
                          f"cnt={cnt}, theta={[theta[t] for t in txs]}, "
                          f"weights={weights.tolist()}, denom={denom:.6f}")
                if denom >= TOLERANCE:
                    fracs = theta_w / denom
                    for i, t in enumerate(txs):
                        contrib = cnt * fracs[i]
                        n[t]   += contrib
                        if verbose and rnd == 0:
                            print(f"    tx{t}: frac={fracs[i]:.4f}, contrib={contrib:.4f}")

        # Normalize
        total     = n.sum()
        theta_new = n / total if total > 0 else theta.copy()

        chg = convergence_check(theta_new, theta)
        if verbose and rnd == 0:
            print(f"\n  Round 1 result: n={n}, total={total:.1f}")
            print(f"  theta_new={theta_new}")
            print()

        theta = theta_new
        if chg == 0 and rnd >= MIN_ROUNDS:
            if verbose:
                print(f"Converged at round {rnd+1}")
            break

    # TPM: rho[t] = theta[t]/eff_len[t], then normalize
    rho = theta / eff_len
    tpm = rho / rho.sum() * 1e6

    if verbose:
        print("\nFinal theta:", theta)
        print("Final TPM  :", tpm)

    return theta, tpm


# ===========================================================
# LK-style EM
# Single-tx ECs SKIPPED in loop, added post-EM as raw counts
# alpha = raw expected counts (not normalized probabilities)
# ===========================================================

def run_lk_em(verbose=True):
    """
    lr-Kallisto style EM (long-read / PacBio branch).

    Single-tx ECs are SKIPPED during the EM loop entirely.
    After convergence, their raw counts are added to alpha.
    alpha is NOT normalized during the loop; it accumulates raw expected counts.
    Effective lengths can be non-uniform; TPM comes from alpha/eff_len normalization.

    Args:
        verbose (bool): print per-round trace.

    Returns:
        np.ndarray: final TPM values.
    """
    n_tx    = 2
    eff_len = np.array(EFF_LEN_LK, dtype=np.float64)
    counts  = np.array(EC_COUNTS, dtype=np.float64)
    alpha   = np.full(n_tx, 1.0 / n_tx)   # uniform init (same as kallisto)

    # Identify single-tx and multi-tx ECs
    single_ecs = [(ec_id, txs, cnt)
                  for ec_id, (txs, cnt) in enumerate(zip(EC_TX, counts))
                  if len(txs) == 1]
    multi_ecs  = [(ec_id, txs, cnt)
                  for ec_id, (txs, cnt) in enumerate(zip(EC_TX, counts))
                  if len(txs) > 1]

    if verbose:
        print("=" * 60)
        print("LK-style EM (single-tx SKIPPED in loop, added post-EM)")
        print(f"eff_lens = {eff_len.tolist()}")
        print("=" * 60)
        print(f"Init alpha = {alpha}")
        print(f"Single-tx ECs (skipped in loop): {[(e[0], e[1]) for e in single_ecs]}")
        print(f"Multi-tx ECs (in loop):          {[(e[0], e[1]) for e in multi_ecs]}")
        print()

    final_round = False
    for rnd in range(N_ROUNDS):
        next_alpha = np.zeros(n_tx)

        for ec_id, txs, cnt in multi_ecs:
            if cnt == 0:
                continue

            weights = np.array([1.0 / eff_len[t] for t in txs])
            denom   = sum(alpha[t] * weights[i] for i, t in enumerate(txs))

            if verbose and rnd == 0:
                print(f"  Round 1, EC {ec_id} (multi-tx {txs}): "
                      f"cnt={cnt}, alpha={[alpha[t] for t in txs]}, "
                      f"weights={weights.tolist()}, denom={denom:.6f}")

            if denom >= TOLERANCE:
                count_norm = cnt / denom
                for i, t in enumerate(txs):
                    contrib = weights[i] * alpha[t] * count_norm
                    next_alpha[t] += contrib
                    if verbose and rnd == 0:
                        print(f"    tx{t}: contrib={contrib:.4f}")

        # Convergence check (only multi-tx contributions visible here)
        chg = convergence_check(next_alpha, alpha)

        if verbose and rnd == 0:
            print(f"\n  Round 1: next_alpha (multi-tx only) = {next_alpha}")
            print(f"  Note: single-tx ECs not yet added (LK adds them after convergence)")
            print()

        # Update alpha (no normalization — raw counts)
        alpha = next_alpha.copy()

        if chg == 0 and rnd >= MIN_ROUNDS:
            if verbose:
                print(f"Converged at round {rnd+1}")
            final_round = True
            break

    # Post-EM: add single-tx counts to alpha (LK long-read behavior)
    if verbose:
        print("\nPost-EM: adding single-tx EC counts to alpha")
    for ec_id, txs, cnt in single_ecs:
        alpha[txs[0]] += cnt
        if verbose:
            print(f"  EC {ec_id} (tx{txs[0]}): alpha[{txs[0]}] += {cnt} → {alpha[txs[0]]:.4f}")

    # compute_rho: alpha / eff_len, then normalize
    rho = alpha / eff_len
    tpm = rho / rho.sum() * 1e6

    if verbose:
        print(f"\nFinal alpha (raw counts) : {alpha}")
        print(f"rho = alpha/eff_len      : {rho}")
        print(f"Final TPM                : {tpm}")

    return alpha, tpm


# ===========================================================
# Side-by-side comparison
# ===========================================================

def main():
    print("\n" + "=" * 60)
    print("SETUP")
    print("=" * 60)
    print("Transcripts: A (index 0), B (index 1)")
    print("ECs:")
    for ec_id, (txs, cnt) in enumerate(zip(EC_TX, EC_COUNTS)):
        names = [TX_NAMES[t] for t in txs]
        kind  = "single-tx" if len(txs) == 1 else "multi-tx"
        print(f"  EC {ec_id}: count={cnt}, transcripts={names} ({kind})")
    print()
    print("JK eff_lens:", EFF_LEN_JK)
    print("LK eff_lens:", EFF_LEN_LK)
    print()

    # Run both
    jk_theta, jk_tpm = run_jk_em(verbose=True)
    print()
    lk_alpha, lk_tpm = run_lk_em(verbose=True)

    # Side-by-side
    print("\n" + "=" * 60)
    print("SIDE-BY-SIDE COMPARISON")
    print("=" * 60)
    print(f"{'Transcript':<20} {'JK_TPM':>12} {'LK_TPM':>12} {'diff%':>10}")
    print("-" * 56)
    for i, name in enumerate(TX_NAMES):
        jk = jk_tpm[i]
        lk = lk_tpm[i]
        ref = max(jk, lk)
        diff_pct = 100.0 * abs(jk - lk) / ref if ref > 0 else 0.0
        print(f"  {name:<18} {jk:>12.2f} {lk:>12.2f} {diff_pct:>9.1f}%")

    print()
    print("KEY INSIGHT:")
    print("  If JK and LK give different TPM for the same input data,")
    print("  the cause is one (or both) of:")
    print("    1. Different eff_len (LK uses real lengths, JK uses 1.0)")
    print("    2. Different single-tx EC handling (LK: post-EM, JK: in-loop)")
    print()
    print("  Try setting EFF_LEN_LK = [1.0, 1.0] in CONFIG to isolate effect #2.")


if __name__ == "__main__":
    main()
