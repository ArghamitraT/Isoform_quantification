"""
multi_sample_em.py
==================
Step 2 of JOLI-Kallisto Phase 2: Multi-sample MAP EM with shared Dirichlet prior.

Implements MultiSampleJoliEM — an outer GD loop that learns a shared Dirichlet
hyperparameter alpha across all samples, with a per-sample inner MAP EM loop
(posterior-mean update).

Architecture:
  Outer loop (GD rounds):
    For each sample s (sequential):
      Run JoliEM with alpha_prior=alpha, init_theta=theta_s (warm-start)
      → theta_s updated in place
    Stack theta_matrix = np.stack(all_theta_s)    # shape (S, T)
    GD step: DirichletOptimizer.update(theta_matrix)
    → alpha updated; check GD convergence

The posterior-mean M-step per round:
  theta_new[t] = (n[t] + alpha[t]) / sum(n + alpha)

where n[t] is the expected count from the E-step (EC-level) and alpha is the
shared Dirichlet concentration vector. Single-tx ECs are handled identically
to single-sample mode (Fix A: raw counts added post-convergence).

Inputs:
  - List of sample directories, each containing bustools TCC output:
      count.mtx, matrix.ec, transcripts.txt, run_info.json
  - All samples must use the same reference transcriptome (same transcripts.txt)
  - eff_len_mode  : "uniform" (Phase 2) | "kallisto" (future)
  - convergence_mode : "joli" (recommended for MAP)
  - GD hyperparameters: gd_lr, alpha_initial, max_gd_rounds, gd_convergence_tol

Outputs (via write_results):
  - Per-sample abundance.tsv written to output_dir/sample_name/
  - alpha_final.npy : learned Dirichlet hyperparameters
  - gd_loss_history.pkl : GD loss per outer round (list of lists)
"""

import os
import pickle

import numpy as np

from load_tcc import load_tcc_data, load_flens, TCCData
from weights import compute_weights, WeightData
from em_algorithm import JoliEM, EMResult
from output_writer import write_abundance
from dirichlet_optimizer import DirichletOptimizer
from training_tracker import TrainingTracker


class MultiSampleJoliEM:
    """
    Multi-sample MAP EM with shared Dirichlet prior learned via gradient descent.

    All samples must share the same reference transcriptome (same transcript
    ordering in transcripts.txt). The transcript space is therefore already
    shared — no union step required.
    """

    def __init__(
        self,
        sample_dirs:      list,
        eff_len_mode:     str   = "uniform",
        convergence_mode: str   = "joli",
        max_em_rounds:    int   = 10000,
        min_em_rounds:    int   = 50,
        max_gd_rounds:    int   = 500,
        gd_lr:            float = 0.01,
        alpha_initial:    float = 1.0,
        gd_convergence_tol: float = 1e-6,
        gd_steps_per_round: int = 10,
    ):
        """
        Load TCC data for all samples and initialize the Dirichlet optimizer.

        Args:
            sample_dirs        : list[str] -- Paths to bustools output directories,
                                             one per sample.
            eff_len_mode       : str   -- "uniform" | "kallisto". Effective length
                                         mode passed to compute_weights().
            convergence_mode   : str   -- "joli" (recommended for MAP) or "kallisto".
                                         Passed to JoliEM.run().
            max_em_rounds      : int   -- Max inner EM rounds per sample per GD round.
            min_em_rounds      : int   -- Min inner EM rounds before convergence check.
            max_gd_rounds      : int   -- Max outer GD iterations.
            gd_lr              : float -- Adam learning rate for alpha.
            alpha_initial      : float -- Starting value for all alpha[t].
            gd_convergence_tol : float -- Stop outer loop when |gd_loss_change| < this.
            gd_steps_per_round : int   -- Adam steps per outer GD round.
        """
        if len(sample_dirs) < 2:
            raise ValueError(
                f"MultiSampleJoliEM requires at least 2 samples, got {len(sample_dirs)}"
            )

        self.sample_dirs        = sample_dirs
        self.eff_len_mode       = eff_len_mode
        self.convergence_mode   = convergence_mode
        self.max_em_rounds      = max_em_rounds
        self.min_em_rounds      = min_em_rounds
        self.max_gd_rounds      = max_gd_rounds
        self.gd_lr              = gd_lr
        self.alpha_initial      = alpha_initial
        self.gd_convergence_tol = gd_convergence_tol
        self.gd_steps_per_round = gd_steps_per_round

        # --- Load TCC data and build JoliEM instance for each sample ---
        print(f"\n[MultiSampleJoliEM] Loading {len(sample_dirs)} samples...")
        self.tcc_data_list    = []
        self.weight_data_list = []
        self.em_list          = []
        self.sample_names     = []

        for idx, sdir in enumerate(sample_dirs):
            sname = os.path.basename(os.path.normpath(sdir))
            self.sample_names.append(sname)
            print(f"\n  [{idx + 1}/{len(sample_dirs)}] {sname}")

            tcc_data = load_tcc_data(sdir)

            if eff_len_mode == "kallisto":
                flens       = load_flens(sdir, n_transcripts=len(tcc_data.transcript_names))
                weight_data = compute_weights(tcc_data, flens=flens, mode="kallisto")
            else:
                weight_data = compute_weights(tcc_data, mode=eff_len_mode)

            self.tcc_data_list.append(tcc_data)
            self.weight_data_list.append(weight_data)
            self.em_list.append(JoliEM(tcc_data, weight_data))

        # Verify all samples share the same number of transcripts
        n_tx_list = [len(td.transcript_names) for td in self.tcc_data_list]
        if len(set(n_tx_list)) != 1:
            raise ValueError(
                f"Samples have different transcript counts: {list(zip(self.sample_names, n_tx_list))}. "
                "All samples must use the same reference transcriptome."
            )
        self.n_transcripts = n_tx_list[0]

        print(f"\n[MultiSampleJoliEM] All samples loaded. "
              f"n_samples={len(sample_dirs)}, n_transcripts={self.n_transcripts}")

        # --- Initialize shared Dirichlet optimizer ---
        self.dirichlet_opt = DirichletOptimizer(
            n_transcripts = self.n_transcripts,
            gd_lr         = gd_lr,
            alpha_initial  = alpha_initial,
        )

        # Per-sample theta (initialized to None; set after first EM run)
        self.theta_list = [None] * len(sample_dirs)

    def run(self) -> dict:
        """
        Run the outer GD loop.

        Each outer round:
          1. For each sample: run MAP EM (warm-start from previous theta).
          2. Stack theta_matrix and run GD step to update alpha.
          3. Check GD convergence.

        Returns:
            dict with keys:
              "theta_list"      : list[np.ndarray] -- per-sample theta at convergence
              "alpha"           : np.ndarray       -- final Dirichlet alpha [T]
              "gd_loss_history" : list[list[float]]-- GD losses per outer round
              "n_gd_rounds"     : int              -- number of GD rounds completed
              "converged"       : bool             -- whether GD converged
        """
        print(f"\n[MultiSampleJoliEM] Starting outer GD loop: "
              f"max_gd_rounds={self.max_gd_rounds}, "
              f"gd_convergence_tol={self.gd_convergence_tol}")

        gd_loss_history = []   # list of loss lists, one per outer round
        prev_gd_loss    = None
        gd_converged    = False

        # Initialize tracker
        tracker = TrainingTracker(self.sample_names)

        for gd_round in range(self.max_gd_rounds):
            print(f"\n[MultiSampleJoliEM] === GD round {gd_round + 1} / {self.max_gd_rounds} ===")

            # Retrieve current shared alpha for this round's MAP EM
            current_alpha = self.dirichlet_opt.get_alpha()   # shape (T,)

            # --- Inner loop: per-sample MAP EM ---
            em_rounds_this_round    = []
            em_converged_this_round = []

            for s_idx, (em, sname) in enumerate(zip(self.em_list, self.sample_names)):
                print(f"\n  [Sample {s_idx + 1}/{len(self.em_list)}: {sname}]")

                result = em.run(
                    max_em_rounds    = self.max_em_rounds,
                    min_rounds       = self.min_em_rounds,
                    convergence_mode = self.convergence_mode,
                    alpha_prior      = current_alpha,
                    init_theta       = self.theta_list[s_idx],   # None on round 0
                )

                # Store normalized theta (not raw alpha) for GD step
                # theta = alpha / alpha.sum() gives normalized abundance
                alpha_s = result.alpha
                total_s = alpha_s.sum()
                self.theta_list[s_idx] = alpha_s / total_s if total_s > 0 else alpha_s

                em_rounds_this_round.append(result.n_rounds)
                em_converged_this_round.append(result.converged)

            # --- GD step: update shared alpha ---
            # Stack per-sample normalized theta into (S, T) matrix
            theta_matrix = np.stack(self.theta_list)   # shape (S, T)

            alpha_updated, round_loss_history = self.dirichlet_opt.update(
                theta_matrix,
                max_iterations = self.gd_steps_per_round,
            )
            gd_loss_history.append(round_loss_history)

            # Check GD convergence on the last loss value of this round
            last_loss = round_loss_history[-1]
            if prev_gd_loss is not None:
                loss_change = abs(prev_gd_loss - last_loss)
                print(f"[MultiSampleJoliEM] GD round {gd_round + 1}: "
                      f"loss={last_loss:.4f}, |change|={loss_change:.2e}")
                if loss_change < self.gd_convergence_tol:
                    gd_converged = True
                    print(f"[MultiSampleJoliEM] GD converged at round {gd_round + 1}.")
            else:
                print(f"[MultiSampleJoliEM] GD round {gd_round + 1}: loss={last_loss:.4f}")

            prev_gd_loss = last_loss

            # --- Record metrics for this round ---
            tracker.record(
                gd_round          = gd_round,
                theta_list        = self.theta_list,
                alpha             = self.dirichlet_opt.get_alpha(),
                em_rounds_list    = em_rounds_this_round,
                em_converged_list = em_converged_this_round,
                gd_loss           = last_loss,
            )
            tracker.print_round_summary(gd_round)

            if gd_converged:
                break

        if not gd_converged:
            print(f"[MultiSampleJoliEM] Reached max_gd_rounds={self.max_gd_rounds} "
                  f"without GD convergence.")

        final_alpha = self.dirichlet_opt.get_alpha()
        print(f"\n[MultiSampleJoliEM] Done. "
              f"GD rounds={gd_round + 1}, converged={gd_converged}, "
              f"alpha: min={final_alpha.min():.4f}, max={final_alpha.max():.4f}, "
              f"mean={final_alpha.mean():.4f}")

        return {
            "theta_list":      self.theta_list,
            "alpha":           final_alpha,
            "gd_loss_history": gd_loss_history,
            "n_gd_rounds":     gd_round + 1,
            "converged":       gd_converged,
            "tracker":         tracker,
        }

    def write_results(self, output_dir: str, results: dict) -> None:
        """
        Write per-sample abundance.tsv, alpha_final.npy, and gd_loss_history.pkl.

        Args:
            output_dir : str  -- Base output directory (timestamped experiment folder).
            results    : dict -- Return value of run().
        """
        os.makedirs(output_dir, exist_ok=True)

        # Per-sample abundance.tsv
        for s_idx, sname in enumerate(self.sample_names):
            sample_out = os.path.join(output_dir, sname)
            os.makedirs(sample_out, exist_ok=True)

            # Convert normalized theta back to raw expected counts for output_writer
            theta_s    = results["theta_list"][s_idx]
            tcc_data   = self.tcc_data_list[s_idx]
            weight_data = self.weight_data_list[s_idx]

            # Reconstruct raw alpha: theta * total_multi_reads + single-tx counts
            # We use the alpha stored in the last EMResult via a final EM pass
            # instead of reconstructing manually — just re-run one MAP EM round
            # to get a proper EMResult for output_writer.
            em     = self.em_list[s_idx]
            result = em.run(
                max_em_rounds    = 1,
                min_rounds       = 1,
                convergence_mode = self.convergence_mode,
                alpha_prior      = results["alpha"],
                init_theta       = theta_s,
            )

            out_path = os.path.join(sample_out, "abundance.tsv")
            write_abundance(
                alpha            = result.alpha,
                eff_lens         = weight_data.eff_lens,
                transcript_names = tcc_data.transcript_names,
                output_path      = out_path,
            )
            print(f"[MultiSampleJoliEM] Written: {out_path}")

        # Shared alpha
        alpha_path = os.path.join(output_dir, "alpha_final.npy")
        np.save(alpha_path, results["alpha"])
        print(f"[MultiSampleJoliEM] Saved alpha: {alpha_path}")

        # GD loss history
        loss_path = os.path.join(output_dir, "gd_loss_history.pkl")
        with open(loss_path, "wb") as fh:
            pickle.dump(results["gd_loss_history"], fh)
        print(f"[MultiSampleJoliEM] Saved GD loss history: {loss_path}")

        # Training tracker
        tracker = results.get("tracker")
        if tracker is not None:
            stats_path = os.path.join(output_dir, "training_stats.pkl")
            tracker.save(stats_path)

            # Generate figures — imported here to keep matplotlib out of the
            # critical training path (only loaded when writing results)
            try:
                import sys
                sys.path.insert(0, os.path.join(
                    os.path.dirname(os.path.abspath(__file__)), "..", "analysis"
                ))
                from plot_training import generate_all_plots
                figures_dir = os.path.join(output_dir, "figures")
                generate_all_plots(tracker, figures_dir)
            except Exception as e:
                print(f"[MultiSampleJoliEM] Warning: figure generation failed: {e}")
