"""
em_algorithm.py
===============
Step 1.3 of JOLI-Kallisto Phase 1: EC-based EM algorithm.

Implements JoliEM — a vectorized EM that operates on equivalence classes (ECs)
rather than individual reads. Mathematically equivalent to kallisto's C++ EM
(EMAlgorithm.h) but written in Python/NumPy, and designed to extend to
Dirichlet MAP (Phase 2), multi-sample GD (Phase 3), and VI (Phase 4).

Single-transcript EC handling (Fix A — matches kallisto EMAlgorithm.h exactly):
  Single-tx ECs are excluded from the EM loop entirely.  After convergence,
  their raw read counts are added directly to alpha (expected counts).  This
  prevents single-tx counts from biasing the fractional assignment of
  multi-tx ECs in cases where a transcript appears in both EC types.

Variable mapping to kallisto EMAlgorithm.h:
  alpha       <->  alpha_          (raw expected counts; NOT normalized)
  theta       <->  (internal only) normalized abundances used for convergence
  n           <->  next_alpha      (accumulated expected counts per M-step)
  ec_counts   <->  counts_         (read count per EC)
  ec_weights  <->  weight_map_     (1/eff_len per transcript per EC)
  TOLERANCE   <->  TOLERANCE       (numeric floor for denominator)

Variable mapping to AT_code EM_VIorMAP_GD_vector.py:
  alpha       <->  all_theta[sample] (after scaling back to raw counts)
  n           <->  all_n[sample]
  phi         <->  all_Phi_ri      (fractional assignment; implicit here via n)

Inputs:
  - TCCData   (from load_tcc.py)
  - WeightData (from weights.py)
  - EM hyperparameters (max_em_rounds, min_rounds, convergence thresholds)

Outputs:
  - EMResult dataclass:
      alpha      : np.ndarray (float64, shape [n_transcripts]) -- raw expected counts
                   (multi-tx EM counts + single-tx raw counts; does NOT sum to 1)
      n_rounds   : int  -- number of EM iterations completed
      converged  : bool -- whether the convergence criterion was met
"""

from dataclasses import dataclass

import numpy as np

from load_tcc import TCCData
from weights import WeightData


# ============================================================
# Constants — match kallisto EMAlgorithm.h exactly
# ============================================================

TOLERANCE          = np.finfo(np.float64).tiny   # ~2.2e-308; floor for denominators
ALPHA_LIMIT        = 1e-7    # transcripts below this × 0.1 are zeroed out
ALPHA_CHANGE_LIMIT = 1e-2    # ignore transcripts below this abundance in convergence check
ALPHA_CHANGE       = 1e-2    # 1% relative change threshold for convergence
DEFAULT_MIN_ROUNDS = 50      # minimum EM iterations before convergence is checked


# ============================================================
# Result container
# ============================================================

@dataclass
class EMResult:
    """
    Output from one run of JoliEM.

    Attributes:
        alpha     (np.ndarray, float64, shape [n_transcripts]):
                   Raw expected read counts per transcript. Multi-tx EM expected
                   counts + single-tx raw counts added post-convergence.
                   Does NOT sum to 1. Matches kallisto's alpha_ after run().
                   Use alpha / eff_lens → rho → TPM for output (see output_writer.py).

        n_rounds  (int): Number of EM iterations completed.

        converged (bool): True if the convergence criterion was satisfied before
                          reaching max_em_rounds.
    """
    alpha: np.ndarray
    n_rounds: int
    converged: bool


# ============================================================
# JoliEM class
# ============================================================

class JoliEM:
    """
    EC-based EM algorithm for isoform quantification.

    The E-step and M-step are combined in each iteration (same structure as
    kallisto's EMAlgorithm.h::run()):

      For each EC with count > 0:
        denom = sum_t ( theta[t] * ec_weight[t] )
        For each t in EC:
          n[t] += ec_count * theta[t] * ec_weight[t] / denom

    Single-transcript ECs are handled deterministically (no division needed),
    matching kallisto's optimization for singleton ECs.

    Multi-transcript ECs are vectorized using pre-built flat index arrays to
    avoid a Python loop over EC positions. See _preprocess() for details.
    """

    def __init__(self, tcc_data: TCCData, weight_data: WeightData):
        """
        Pre-process TCC and weight data into vectorization-friendly flat arrays.

        Args:
            tcc_data    : TCCData    -- Parsed bustools output.
            weight_data : WeightData -- Effective lengths and EC weights.
        """
        self.n_transcripts  = len(tcc_data.transcript_names)
        self.ec_counts      = tcc_data.ec_counts          # shape (n_ecs,)
        self.ec_transcripts = tcc_data.ec_transcripts     # list[list[int]]
        self.eff_lens       = weight_data.eff_lens        # shape (n_transcripts,)
        self.ec_weights     = weight_data.ec_weights      # list[np.ndarray]

        self._preprocess()

    def _preprocess(self):
        """
        Split ECs into single-transcript and multi-transcript groups.
        Build flat index arrays for vectorized multi-tx EM computation.

        Single-tx ECs: deterministic assignment, handled with np.add.at once.
        Multi-tx ECs:  fractional assignment, vectorized via flat arrays.

        Flat arrays built here (all parallel, indexed by EC position):
          _multi_flat_tx      : transcript index for each (ec, position) pair
          _multi_flat_weights : weight for each (ec, position) pair
          _multi_flat_ec_idx  : which multi-EC each position belongs to
          _multi_ec_counts    : count for each multi-EC (one per EC, not per position)
        """
        # --- Separate single-tx and multi-tx ECs ---
        single_ec_ids = []
        multi_ec_ids  = []
        for ec_id, txs in enumerate(self.ec_transcripts):
            if len(txs) == 1:
                single_ec_ids.append(ec_id)
            else:
                multi_ec_ids.append(ec_id)

        # Single-tx: static arrays (independent of theta, computed once)
        if single_ec_ids:
            single_ec_ids = np.array(single_ec_ids, dtype=np.int64)
            self._single_tx_ids    = np.array(
                [self.ec_transcripts[i][0] for i in single_ec_ids], dtype=np.int64)
            self._single_ec_counts = self.ec_counts[single_ec_ids]  # shape (n_single,)
        else:
            self._single_tx_ids    = np.array([], dtype=np.int64)
            self._single_ec_counts = np.array([], dtype=np.float64)

        # Multi-tx: build flat arrays
        flat_tx      = []
        flat_weights = []
        flat_ec_idx  = []
        multi_counts = []

        for local_idx, ec_id in enumerate(multi_ec_ids):
            txs     = self.ec_transcripts[ec_id]
            weights = self.ec_weights[ec_id]           # np.ndarray, shape (len(txs),)
            for pos, (tx, w) in enumerate(zip(txs, weights)):
                flat_tx.append(tx)
                flat_weights.append(w)
                flat_ec_idx.append(local_idx)
            multi_counts.append(self.ec_counts[ec_id])

        if flat_tx:
            self._multi_flat_tx      = np.array(flat_tx,      dtype=np.int64)
            self._multi_flat_weights = np.array(flat_weights,  dtype=np.float64)
            self._multi_flat_ec_idx  = np.array(flat_ec_idx,   dtype=np.int64)
            self._multi_ec_counts    = np.array(multi_counts,   dtype=np.float64)
            self._n_multi_ecs        = len(multi_ec_ids)
        else:
            self._multi_flat_tx      = np.array([], dtype=np.int64)
            self._multi_flat_weights = np.array([], dtype=np.float64)
            self._multi_flat_ec_idx  = np.array([], dtype=np.int64)
            self._multi_ec_counts    = np.array([], dtype=np.float64)
            self._n_multi_ecs        = 0

        print(f"[JoliEM] Pre-processed: {len(single_ec_ids)} single-tx ECs, "
              f"{self._n_multi_ecs} multi-tx ECs, "
              f"{len(self._multi_flat_tx)} total EC-transcript positions")

    def _em_step(self, theta: np.ndarray) -> np.ndarray:
        """
        Run one combined E+M step and return updated n (expected counts per tx).

        Only multi-transcript ECs are processed here.  Single-tx ECs are excluded
        from the EM loop and added as raw counts post-convergence (Fix A — matches
        kallisto EMAlgorithm.h long-read branch exactly).

        The E-step computes fractional assignment weights; the M-step accumulates
        them into n[t]. Both happen in one pass, matching kallisto's structure.

        Args:
            theta : np.ndarray (float64, shape [n_transcripts]) -- current normalized
                    abundances (derived from multi-tx EM only; excludes single-tx).

        Returns:
            np.ndarray (float64, shape [n_transcripts]) -- n[t] = expected counts
            from multi-tx ECs only.
        """
        n = np.zeros(self.n_transcripts, dtype=np.float64)

        # --- Multi-transcript ECs: vectorized fractional assignment ---
        # Single-tx ECs are intentionally excluded here; added post-convergence.
        if self._n_multi_ecs > 0:
            # theta values and weights for every (ec, position) in flat arrays
            theta_per_pos   = theta[self._multi_flat_tx]              # shape (n_positions,)
            weighted_theta  = theta_per_pos * self._multi_flat_weights # theta[t] * w[t,ec]

            # Denominator per multi-EC: sum of weighted_theta over positions in each EC
            denominators = np.zeros(self._n_multi_ecs, dtype=np.float64)
            np.add.at(denominators, self._multi_flat_ec_idx, weighted_theta)

            # Broadcast denominator and EC count back to each position
            denom_per_pos = denominators[self._multi_flat_ec_idx]    # shape (n_positions,)
            count_per_pos = self._multi_ec_counts[self._multi_flat_ec_idx]

            # Mask positions where denominator is too small (avoids div-by-zero)
            valid = denom_per_pos >= TOLERANCE
            contributions = np.where(
                valid,
                count_per_pos * weighted_theta / np.where(valid, denom_per_pos, 1.0),
                0.0
            )

            # Accumulate contributions into n
            np.add.at(n, self._multi_flat_tx, contributions)

        return n

    def run(
        self,
        max_em_rounds: int = 10000,
        min_rounds: int    = DEFAULT_MIN_ROUNDS,
    ) -> EMResult:
        """
        Run EM until convergence or max_em_rounds, then zero small abundances.

        Convergence criterion (matches kallisto EMAlgorithm.h exactly):
          Count transcripts where:
            theta_new[t] > ALPHA_CHANGE_LIMIT  AND
            |theta_new[t] - theta[t]| / theta_new[t] > ALPHA_CHANGE
          Stop when this count == 0 AND round > min_rounds.

        After convergence, zero out small multi-tx abundances, then build alpha:
          alpha = theta * total_multi_reads  +  single_tx_raw_counts
        This matches kallisto's post-EM step where single-tx counts are added
        directly to alpha_ without going through the EM fractional assignment.

        Args:
            max_em_rounds : int -- Maximum EM iterations (default: 10000).
            min_rounds    : int -- Minimum iterations before convergence check (default: 50).

        Returns:
            EMResult -- alpha (raw expected counts), n_rounds, converged.
        """
        print(f"\n[JoliEM] Starting EM: max_rounds={max_em_rounds}, "
              f"min_rounds={min_rounds}, n_transcripts={self.n_transcripts}")

        # --- Initialization: uniform, matching kallisto's alpha_[i] = 1/num_targets ---
        theta     = np.full(self.n_transcripts, 1.0 / self.n_transcripts,
                            dtype=np.float64)
        converged = False

        for round_num in range(max_em_rounds):
            # One E+M step
            n = self._em_step(theta)

            # M-step: normalize to get new theta
            total = n.sum()
            if total > 0:
                theta_new = n / total
            else:
                # No reads assigned — keep current theta (degenerate case)
                theta_new = theta.copy()

            # Convergence check (matches kallisto exactly)
            changed = int(np.sum(
                (theta_new > ALPHA_CHANGE_LIMIT) &
                (np.abs(theta_new - theta) / np.maximum(theta_new, TOLERANCE)
                 > ALPHA_CHANGE)
            ))

            theta = theta_new

            if changed == 0 and round_num >= min_rounds:
                converged = True
                print(f"[JoliEM] Converged at round {round_num + 1}")
                break

            # Progress checkpoint every 100 rounds
            if (round_num + 1) % 100 == 0:
                print(f"[JoliEM] Round {round_num + 1}: "
                      f"changed={changed}, "
                      f"nonzero_tx={int((theta > 0).sum())}")

        if not converged:
            print(f"[JoliEM] Reached max_em_rounds={max_em_rounds} without convergence.")

        # Zero out very small multi-tx abundances on theta before scaling.
        # Threshold ALPHA_LIMIT/10 = 1e-8 on normalized theta is equivalent to
        # zeroing raw counts below ~1e-8 * total_multi_reads (matches kallisto).
        n_zeroed = int(((theta > 0) & (theta < ALPHA_LIMIT / 10)).sum())
        theta[theta < ALPHA_LIMIT / 10] = 0.0
        if n_zeroed:
            print(f"[JoliEM] Zeroed {n_zeroed} multi-tx transcripts below {ALPHA_LIMIT/10:.1e}")

        # --- Build alpha: raw expected counts (matches kallisto alpha_ after run()) ---
        # Scale multi-tx theta back to expected counts using total multi-tx reads.
        total_multi_reads = float(self._multi_ec_counts.sum()) if self._n_multi_ecs > 0 else 0.0
        alpha = theta * total_multi_reads

        # Add single-tx raw counts post-convergence (Fix A).
        # Each single-tx EC assigns its count deterministically to its one transcript.
        if len(self._single_tx_ids) > 0:
            np.add.at(alpha, self._single_tx_ids, self._single_ec_counts)

        print(f"[JoliEM] Done. Rounds={round_num + 1}, "
              f"nonzero_transcripts={int((alpha > 0).sum())}")

        return EMResult(alpha=alpha, n_rounds=round_num + 1, converged=converged)
