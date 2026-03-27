# JOLI-Kallisto: Multi-Sample MAP EM Implementation Plan
**Date:** 2026-03-26
**Last updated:** 2026-03-26

---

## Overview

Extend JOLI-Kallisto from single-sample plain EM (Phase 1) to multi-sample MAP EM with a
shared Dirichlet prior optimized via gradient descent (Phase 2). This mirrors the AT_code
experiment-4 approach but operates on EC-level TCC data (kallisto/bustools output) rather
than read-level BAM alignments.

**Both single-sample and multi-sample pipelines are available.** Single-sample
(`main_joli.py` / `run_joli_kallisto.sh`) is unchanged. Multi-sample is a new entry point
built on top of the same core modules.

**Both pipelines support long reads and short reads.** The EM math is identical; only the
upstream bustools call differs (`--long` for long reads, paired R1+R2 for short reads).
A `READ_TYPE` flag in the bash scripts switches the bustools command accordingly.

**Why multi-sample?** A shared Dirichlet prior across samples regularizes isoform abundance
estimates: transcripts consistently expressed across samples get stronger support, while
noise-driven assignments are suppressed. The prior (alpha) is learned from data via GD,
so it adapts to the dataset rather than being hand-tuned.

---

## Math

### Single-sample plain EM (Phase 1 — unchanged)

E-step + M-step per round:
```
n[t] = sum over multi-tx ECs { ec_count * theta[t] * w[t,ec] / sum_t'(theta[t'] * w[t',ec]) }
theta_new[t] = n[t] / sum_t(n[t])
```

### Multi-sample posterior-mean EM (Phase 2 — new)

The only change is the M-step: add the shared Dirichlet prior alpha to the numerator.
This is the **Dirichlet posterior mean** — the same update AT_code uses in
`update_theta_vectorized`:

```
theta_new[t] = (n[t] + alpha[t]) / sum_t(n[t] + alpha[t])
```

where `alpha` is the shared Dirichlet concentration vector (shape [T], one value per transcript).

- When alpha[t] = 1 for all t → reduces to plain EM (uniform prior = no regularization)
- When alpha[t] > 1 → pushes theta toward uniform (shrinks small estimates)
- When alpha[t] < 1 → sparsity-inducing (pulls small estimates further toward zero)
- alpha is shared across ALL samples; theta is per-sample

The posterior mean is always valid (numerator always ≥ 0 as long as alpha > 0), which is
why we use it over the strict MAP formula `(n + alpha - 1) / (N + sum(alpha) - T)` — the
strict MAP can produce negative numerators when alpha[t] < 1 and n[t] = 0.

### GD step (outer loop)

After each inner EM round over all samples, update alpha to maximize the Dirichlet
log-likelihood summed over all samples:
```
L(alpha) = sum_s [ log Gamma(sum_t alpha[t]) - sum_t log Gamma(alpha[t])
                   + sum_t (alpha[t] - 1) * log theta_s[t] ]
```
Optimize `log_alpha` (unconstrained) using Adam → `alpha = exp(log_alpha)`.
This keeps alpha strictly positive throughout training.

---

## Key Simplification vs AT_code

| Aspect | AT_code | JOLI multi-sample |
|--------|---------|-------------------|
| Data per sample | Read-level BAM (Yri matrix) | EC-level TCC (count.mtx) |
| Transcript space | Union of BAM-observed isoforms (runtime `np.searchsorted`) | Same `transcripts.txt` for all samples — no union step needed |
| theta vector | Sparse per sample (only observed isoforms) | Dense: length T (all transcripts in reference) |
| DirichletOptimizer alignment | `np.searchsorted` to map sample isoforms → union index | `np.stack(all_theta)` — O(1), no alignment needed |
| Inner EM | Read-to-transcript assignments (Phi_ri) | EC-to-transcript assignments (n from `_em_step`) |
| theta between GD rounds | Re-initialized each time | Carried over from previous round naturally |

The shared transcript index means the GD bottleneck (alignment/searchsorted) disappears entirely.
theta is simply carried over between GD rounds — no re-initialization needed.

---

## Architecture

### Outer loop

```
Initialize alpha (shape [T], all = alpha_initial)
Initialize theta_s = None for each sample s

For each GD round (outer_round = 0, 1, ..., max_gd_rounds):

    For each sample s (sequential loop):
        # Round 0: start from uniform; rounds 1+: warm-start from previous theta_s
        theta_s = run_posterior_mean_em(
            sample_s,
            alpha_prior    = alpha,
            init_theta     = theta_s,          # None on round 0 → uniform init
            max_em_rounds  = max_em_rounds,    # round 0: full convergence
        )

    theta_matrix = np.stack(all_theta_s)       # shape (S, T)
    alpha, gd_loss = dirichlet_optimizer.update(theta_matrix)

    log GD loss, check convergence
    if |gd_loss - prev_gd_loss| < gd_convergence_tol:
        break

Output: theta_s per sample → write abundance.tsv per sample
```

theta_s is naturally carried over between GD rounds (no special mechanism): after each
round, theta_s holds the converged value under the current alpha. On the next round,
it starts from there rather than uniform, so convergence is fast (typically < 100 rounds).

### No parallelization (for now — planned for future)

Inner EM per sample runs sequentially. Simple to implement and debug.
**Future:** add `multiprocessing.Pool` to parallelize the per-sample inner EM loop.
Each sample's EM is fully independent given a fixed alpha, so parallelization is
straightforward — just replace the sequential loop with `Pool.map()`. Add this once
the sequential version is validated and profiling confirms it is the bottleneck.

---

## Input / Output

### Input (per sample)
- **Long read**: one FASTQ file per sample
- **Short read**: paired R1 + R2 FASTQ files per sample
- Both go through the same bustools pipeline; `READ_TYPE` flag controls the command

### Output (per sample)
- `abundance.tsv` — per-transcript quantification (same format as single-sample)

### Experiment folder structure
```
files/results/exprmnt_{timestamp}/
├── experiment_description.log
├── running.log
├── runtime.txt
├── code_snapshot/
├── gd_loss_history.pkl         ← GD loss per outer round (for convergence plots)
├── alpha_final.npy             ← learned Dirichlet hyperparameters [T]
├── sample_1/
│   └── abundance.tsv
├── sample_2/
│   └── abundance.tsv
└── ...
```

---

## Files to Create / Modify

### New files

| File | Role |
|------|------|
| `multi_sample_em.py` | `MultiSampleJoliEM` class — outer GD loop, sequential inner MAP EM per sample |
| `dirichlet_optimizer.py` | `DirichletOptimizer` — simplified from AT_code; dense theta_matrix input (S × T), no alignment step |
| `main_multisample_joli.py` | CLI entry point with CONFIG: sample dirs, READ_TYPE, GD lr, alpha_initial, max_gd_rounds, convergence tol |
| `run_multisample_joli.sh` | Bash pipeline: bustools per sample → JOLI multi-sample EM; supports long + short read via READ_TYPE flag |
| `submit_multisample_joli.sh` | Thin SLURM submission wrapper for run_multisample_joli.sh |

### Modified files

| File | Change |
|------|--------|
| `em_algorithm.py` | Add `alpha_prior` (optional, shape [T]) and `init_theta` (optional, shape [T]) to `JoliEM.run()`. If `alpha_prior=None` → plain EM unchanged. If provided → posterior-mean M-step. If `init_theta` provided → start from that theta instead of uniform. |
| `run_joli_kallisto.sh` | Add `READ_TYPE` flag (long/short) to control bustools command. Currently long-read only. |
| `run_lr_kallisto.sh` | Same `READ_TYPE` flag for consistency. |

### Unchanged files

`load_tcc.py`, `weights.py`, `output_writer.py`, `main_joli.py` — no changes needed.

---

## `em_algorithm.py` changes (detail)

```python
def run(self,
        max_em_rounds: int      = 10000,
        min_rounds: int         = DEFAULT_MIN_ROUNDS,
        convergence_mode: str   = "kallisto",
        alpha_prior: np.ndarray = None,   # NEW: Dirichlet prior [T]; None → plain EM
        init_theta: np.ndarray  = None,   # NEW: starting theta [T]; None → uniform
) -> EMResult:

    # Initialization
    if init_theta is not None:
        theta = init_theta.copy()
        theta = theta / theta.sum()       # re-normalize for safety
    else:
        theta = np.full(self.n_transcripts, 1.0 / self.n_transcripts, dtype=np.float64)

    for round_num in range(max_em_rounds):
        n = self._em_step(theta)

        # M-step: posterior-mean EM (with prior) or plain EM (no prior)
        if alpha_prior is not None:
            numerator = n + alpha_prior
        else:
            numerator = n

        total = numerator.sum()
        theta_new = numerator / total if total > 0 else theta.copy()

        # ... convergence check, zeroing, finalRound — unchanged ...
```

---

## `dirichlet_optimizer.py` (detail)

Adapted from `AT_code/DirichletOptimizer_vector.py`. Key change: input is already a
dense (S × T) matrix — no searchsorted or alignment needed.

```python
class DirichletOptimizer:
    def __init__(self, n_transcripts: int, gd_lr: float, alpha_initial: float = 1.0):
        log_alpha_init = np.full(n_transcripts, np.log(alpha_initial), dtype=np.float32)
        self.log_alpha = torch.tensor(log_alpha_init, requires_grad=True)
        self.optimizer = torch.optim.Adam([self.log_alpha], lr=gd_lr)

    def update(self,
               theta_matrix: np.ndarray,  # shape (S, T) — dense, all samples aligned
               max_iterations: int = 10,
               tolerance: float = 1e-6,
    ) -> tuple[np.ndarray, list[float]]:
        """
        Run max_iterations Adam steps to maximize Dirichlet log-likelihood
        over all samples. Returns updated alpha (T,) and loss history.
        """
        data = torch.tensor(theta_matrix, dtype=torch.float32)   # (S, T)
        loss_history = []
        for _ in range(max_iterations):
            self.optimizer.zero_grad()
            alpha = torch.exp(self.log_alpha)                     # keep alpha > 0
            dirichlet = torch.distributions.Dirichlet(alpha)
            log_likelihood = dirichlet.log_prob(data).sum()       # sum over samples
            loss = -log_likelihood
            loss.backward()
            self.optimizer.step()
            loss_history.append(loss.item())
        alpha_np = torch.exp(self.log_alpha).detach().numpy()
        return alpha_np, loss_history
```

---

## Sample Input Options

Two mutually exclusive ways to specify samples in the CONFIG — the code resolves
whichever is set:

**Option A — explicit list** (use when samples are spread across different directories):
```python
SAMPLE_DIRS = [
    "/path/to/kallisto_output/sample_1/",
    "/path/to/kallisto_output/sample_2/",
    "/path/to/kallisto_output/sample_3/",
    # add as many samples as needed
]
SAMPLES_FOLDER = ""   # leave empty when using SAMPLE_DIRS
```

**Option B — folder scan** (use when all samples live under one parent directory):
```python
SAMPLE_DIRS    = []   # leave empty when using SAMPLES_FOLDER
SAMPLES_FOLDER = "/path/to/kallisto_output/"
# code automatically discovers all subdirectories under this folder
# each subdirectory is treated as one sample
```

Resolution logic (in `main_multisample_joli.py`):
```python
if SAMPLES_FOLDER:
    sample_dirs = sorted([
        os.path.join(SAMPLES_FOLDER, d)
        for d in os.listdir(SAMPLES_FOLDER)
        if os.path.isdir(os.path.join(SAMPLES_FOLDER, d))
    ])
else:
    sample_dirs = SAMPLE_DIRS
```

---

## `main_multisample_joli.py` CONFIG

```python
# ============================================================
# CONFIG — edit everything here; do not touch pipeline logic below
# ============================================================

# --- Sample input: use ONE of the two options below ---
SAMPLE_DIRS    = [                        # Option A: explicit list of sample dirs
    "/path/to/kallisto_output/sample_1/",
    "/path/to/kallisto_output/sample_2/",
]
SAMPLES_FOLDER = ""                       # Option B: parent folder; leave "" if using SAMPLE_DIRS
#                                         #   code will loop over all subdirs automatically

RESULTS_BASE       = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"
READ_TYPE          = "long"          # "long" (PacBio/ONT) | "short" (paired Illumina)
EFF_LEN_MODE       = "uniform"       # "uniform" | "kallisto"
MAX_EM_ROUNDS      = 10000           # max inner EM rounds per sample per GD round
MIN_EM_ROUNDS      = 50              # min inner EM rounds before convergence check
CONVERGENCE_MODE   = "joli"          # "joli" (normalized theta, faster for MAP)
MAX_GD_ROUNDS      = 500             # max outer GD iterations
GD_LR              = 0.01            # Adam learning rate for alpha
ALPHA_INITIAL      = 1.0             # initial Dirichlet concentration (uniform prior)
GD_CONVERGENCE_TOL = 1e-6            # stop outer loop when |loss change| < this
GD_STEPS_PER_ROUND = 10              # Adam steps per outer GD round
# ============================================================
```

---

## Implementation Order

```
Step 1 — Modify em_algorithm.py
    → Add alpha_prior and init_theta parameters to JoliEM.run()
    → Posterior-mean M-step when alpha_prior is not None
    → Plain EM unchanged when alpha_prior is None
    → Fully backward-compatible

Step 2 — Create dirichlet_optimizer.py
    → Adapt AT_code/DirichletOptimizer_vector.py
    → Remove alignment/searchsorted (dense input, no alignment needed)
    → Input: theta_matrix (S × T numpy array)
    → Output: alpha (T,), loss history list

Step 3 — Create multi_sample_em.py
    → MultiSampleJoliEM class
    → Loads list of sample dirs → TCCData + WeightData per sample
    → Sequential outer GD loop with per-sample inner EM
    → Saves gd_loss_history.pkl, alpha_final.npy, per-sample abundance.tsv

Step 4 — Create main_multisample_joli.py
    → CLI entry point with CONFIG section
    → Saves experiment outputs (experiment_description.log, runtime.txt, code_snapshot)

Step 5 — Add READ_TYPE to bash scripts
    → run_joli_kallisto.sh: add READ_TYPE flag to control bustools --long vs paired mode
    → run_multisample_joli.sh: new bash pipeline for multi-sample

Step 6 — Test
    → Run on 2 samples
    → Verify GD loss decreases monotonically
    → Verify per-sample correlation vs single-sample plain EM stays high (Spearman > 0.99)
    → Verify MAP shifts abundance for consistently-expressed transcripts
```

---

## Open Questions / Decisions

| Question | Decision | Notes |
|----------|----------|-------|
| Posterior mean vs strict MAP? | Posterior mean `(n+alpha)/(N+sum(alpha))` | Always valid; same as AT_code; strict MAP risks negative numerator when alpha < 1 |
| Parallelize inner EM? | No — sequential loop | Add later if profiling shows it is the bottleneck |
| theta between GD rounds? | Carry over naturally | No re-initialization; fast convergence in rounds 1+ |
| VI mode (Phase 4)? | Not in this phase | VI replaces point-estimate theta with variational distribution `q(theta) ~ Dirichlet(alpha')`. E-step uses `E[log theta] = psi(alpha') - psi(sum(alpha'))` (digamma) instead of theta directly. GD maximizes ELBO instead of Dirichlet log-likelihood. Requires per-sample alpha' parameter. Plan separately. |
| Parallelize inner EM? | Sequential for now | Future: replace sample loop with `multiprocessing.Pool.map()`. Each sample's EM is fully independent given alpha so the change is localized. Add after the sequential version is validated. |
| alpha floor? | exp(log_alpha) always > 0 | Parameterizing log_alpha ensures positivity throughout |

---

## Comparison to AT_code

| AT_code | JOLI multi-sample |
|---------|-------------------|
| `EM_VIorMAP_GD_vector.py` | `multi_sample_em.py` |
| `DirichletOptimizer_vector.py` | `dirichlet_optimizer.py` (simplified) |
| `main_EM_VIorMAP_GD_vector.py` | `main_multisample_joli.py` |
| BAM → Yri (read × transcript) | count.mtx → TCCData (EC × transcript) |
| Runtime union of observed isoforms + `np.searchsorted` | Shared `transcripts.txt` — `np.stack(theta_list)`, O(1) |
| `update_theta_vectorized` M-step | `JoliEM.run(alpha_prior=alpha)` M-step |
