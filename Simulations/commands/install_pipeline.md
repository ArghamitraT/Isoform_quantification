# /install-simulation — Install the RNA-seq Simulation Pipeline

You are implementing the full RNA-seq simulation pipeline from scratch inside `code/Simulations/`.
Before doing anything, read:
- `code/Simulations/plans/simulation_pipeline_plan.md` — authoritative spec for every file
- `code/Simulations/REFERENCE.md` — original pipeline reference

**NEVER touch `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation`.**
**NEVER modify the `lrgsp_simulation` conda environment.**

---

## What This Skill Does

Installs the complete pipeline implementation by writing every missing source file in the
correct order. Safe to re-run: checks each file before writing — skips files that already exist
and look complete (have a top-level docstring), overwrites stubs.

---

## Environment Facts

- **Conda env:** `lrgsp_simulation`
- **Packages in env:** numpy, pandas, scipy, pysam, biopython, matplotlib, seaborn, tqdm, nanocount, torch, pyro-ppl
- **No gffutils** — GTF parsing must use pandas only
- **Modules (loaded in SLURM scripts, NOT in Python):** minimap2/2.29, RSEM/1.3.3-foss-2022a, STAR/2.7.11b-GCC-13.2.0, SAMtools/1.21
- **NanoSim (read-only):** `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation/lrgasp-simulation/src/NanoSim/`
- **IsoSeqSim (read-only):** `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation/lrgasp-simulation/src/isoseqsim/`

---

## Step 0 — Pre-flight Check

Before writing any file, run:

```bash
find /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations/src -type f 2>/dev/null | sort
ls /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations/
conda env list | grep lrgsp_simulation
```

Report which files already exist. For each existing file, read the first 10 lines to check
for a top-level docstring — if it has one, skip it. If it is empty or a stub, overwrite it.

---

## Step 1 — `src/utility.py`

**Purpose:** Shared helpers used by every other module.

Must implement:
- `create_run_dir(base="files/results") -> Path` — creates `exprmnt_{timestamp}` folder, returns its Path
- `save_runtime(run_dir: Path, elapsed_seconds: float)` — writes `runtime.txt`
- `save_code_snapshot(run_dir: Path, src_root: Path)` — copies all `.py .sh .txt .yml .yaml .md` files from `src_root` into `run_dir/code_snapshot/`, preserving relative paths

Requirements:
- Top-level docstring: what it does, inputs, outputs
- Docstring on every function: description, Args, Returns
- Print checkpoint after each helper is called (e.g. `[utility] Run dir created: ...`)
- No CONFIG section needed (it's a library, not an entry point)

---

## Step 2 — `src/prepare_reference.py`

**Purpose:** Phase 1 — GTF + genome FASTA → simulation-ready reference files.

CLI (argparse):
```
--genome       path to GRCh38 FASTA
--gtf          path to GENCODE GTF
--output       output prefix (e.g. "human_sim")
--sqanti       optional SQANTI classification TSV for novel isoforms
--n_novel      int, how many novel isoforms to include (default 0)
--polya_len    int, poly-A tail length in bp (default 100)
```

Steps (with print checkpoints at each):
1. Parse GTF with pandas (no gffutils) — build transcript→exon table (chrom, start, end, strand, exon_number)
2. Extract transcript sequences from genome FASTA using pysam
3. Append poly-A tail (`"A" * polya_len`) to each sequence
4. If `--sqanti` provided: load classification TSV, filter out `not_in_catalog` and `novel_in_catalog` rows, add remaining sequences
5. Write `<output>.transcripts.fasta`, `<output>.annotation.gtf`, `<output>.genome.fasta`, `<output>.novel_isoforms.tsv`

Output summary: print total transcripts written, novel isoforms added (if any).

Requirements: top-level + per-function docstrings, CONFIG at top of `main()`, run instructions in docstring.

---

## Step 3 — `src/generate_abundances.py`

**Purpose:** Phase 2 — generate synthetic TPM abundance profile.

CLI (argparse):
```
--mode         "lognormal" | "custom" | "uniform"
--transcripts  path to .transcripts.fasta (to get transcript IDs)
--output       path for output TSV (transcript_id, TPM)
--seed         int (default 42)
--dropout      float 0–1, fraction with TPM=0 (lognormal only, default 0.3)
--mu           float, log-normal mean (default 2.0)
--sigma        float, log-normal std  (default 2.0)
--custom_tpm   path to TSV (transcript_id, TPM) for mode=custom
--cancer_genes path to file with one gene ID per line to up-regulate
--fold_change  float, up-regulation multiplier for cancer genes (default 5.0)
```

Mode implementations:
- **lognormal:** draw from LogNormal(mu, sigma), zero out `dropout` fraction randomly, normalize to sum=1e6 (TPM)
- **custom:** load provided TSV, optionally multiply cancer gene transcripts by `fold_change`, renormalize
- **uniform:** equal TPM to all transcripts

Output: TSV with columns `transcript_id`, `TPM`.
After writing: print total transcripts, expressed count (TPM > 0), top 10 by TPM.

Requirements: top-level + per-function docstrings, CONFIG at top of `main()`, run instructions in docstring.

---

## Step 4 — `src/simulate_illumina.py`

**Purpose:** Phase 3a — simulate paired-end Illumina reads using RSEM.

CLI (argparse):
```
--reference_prefix   path prefix from Phase 1
--abundances         path to abundances TSV from Phase 2
--output_dir         directory for output files
--threads            int (default 16)
--read_count         int, read pairs to simulate (default 100_000_000)
--seed               int (default 42)
--noise_reads        flag, include background noise
```

Steps:
1. Build RSEM reference from transcript FASTA (`rsem-prepare-reference`)
2. Convert abundances TSV → RSEM-format input
3. Run `rsem-simulate-reads` with paired-end 150 bp model
4. Rename outputs to `Illumina.simulated_1.fq`, `Illumina.simulated_2.fq`
5. Write `Illumina.isoform_counts.tsv` and `Illumina.read_to_isoform.tsv`

All external commands run via `subprocess.run(..., check=True)` with stdout/stderr captured and printed.

Requirements: top-level + per-function docstrings, CONFIG at top of `main()`, run instructions in docstring.

---

## Step 5 — `src/simulate_pacbio.py`

**Purpose:** Phase 3b — simulate PacBio CCS reads using IsoSeqSim.

IsoSeqSim path (read-only): `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation/lrgasp-simulation/src/isoseqsim/`

CLI (argparse):
```
--reference_prefix   path prefix from Phase 1
--abundances         path to abundances TSV from Phase 2
--output_dir         directory for output files
--threads            int (default 16)
--read_count         int, reads to simulate (default 10_000_000)
--seed               int (default 42)
--sub_rate           float, substitution error rate (default 0.004)
--ins_rate           float, insertion error rate    (default 0.006)
--del_rate           float, deletion error rate     (default 0.006)
```

Steps:
1. Convert abundances TSV → IsoSeqSim count input format
2. Run IsoSeqSim via `python <isoseqsim_path>/isoseqsim.py ...`
3. Rename output to `PacBio.simulated.fasta`
4. Write `PacBio.isoform_counts.tsv` and `PacBio.read_to_isoform.tsv`

All external commands via `subprocess.run(..., check=True)`.

Requirements: top-level + per-function docstrings, CONFIG at top of `main()`, run instructions in docstring.

---

## Step 6 — `src/simulate_ont.py`

**Purpose:** Phase 3c — simulate ONT reads using Trans-NanoSim.

NanoSim path (read-only): `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation/lrgasp-simulation/src/NanoSim/`

CLI (argparse):
```
--reference_prefix   path prefix from Phase 1
--abundances         path to abundances TSV from Phase 2
--output_dir         directory for output files
--threads            int (default 16)
--read_count         int, reads to simulate (default 30_000_000)
--seed               int (default 42)
--ont_type           "cDNA" | "dRNA" (default "cDNA")
--noise_reads        flag, include unaligned/noise reads
```

Steps:
1. Convert abundances TSV → NanoSim expression profile format
2. Run Trans-NanoSim via `python <nanosim_path>/simulator.py transcriptome ...`
3. Rename output to `ONT.simulated.fastq`
4. Write `ONT.isoform_counts.tsv` and `ONT.read_to_isoform.tsv`

All external commands via `subprocess.run(..., check=True)`.

Requirements: top-level + per-function docstrings, CONFIG at top of `main()`, run instructions in docstring.

---

## Step 7 — `main.py`

**Purpose:** Entry point tying all phases together with a single CONFIG section.

CONFIG block (at very top of `main()`, clearly delimited):
```python
# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================
reference_genome   = "/path/to/GRCh38.fa"
reference_gtf      = "/path/to/gencode.vXX.annotation.gtf"
output_prefix      = "human_sim"
abundance_mode     = "lognormal"    # "lognormal" | "custom" | "uniform"
custom_tpm_file    = ""
cancer_genes_file  = ""
cancer_fold_change = 5.0
dropout_fraction   = 0.3
illumina_count     = 100_000_000
pb_count           = 10_000_000
ont_count          = 30_000_000
ont_type           = "cDNA"
threads            = 16
seed               = 42
noise_reads        = False
test_mode          = False          # True = low read counts for quick testing
# ============================================================
```

Flow: Phase 1 → Phase 2 → Phase 3 (all three simulators). Uses `utility.py` for run dir,
runtime, and code snapshot. Writes `experiment_description.log` with full config dump.
First line of `running.log`: `Script: <path/to/main.py>`.

---

## Step 8 — SLURM Scripts

Write three scripts: `submit_illumina.sh`, `submit_pacbio.sh`, `submit_ont.sh`.

Each must:
- Have a CONFIG section at the top (paths, read counts, threads)
- Load the correct modules (`module load ...`)
- Activate `lrgsp_simulation` via conda
- Set `#SBATCH` directives (see resources below)
- Write first line of log: `echo "Script: $(realpath "$0")" > "${LOG}"`

| Script | CPUs | RAM | Modules |
|--------|------|-----|---------|
| submit_illumina.sh | 16 | 64G | RSEM/1.3.3-foss-2022a, SAMtools/1.21 |
| submit_pacbio.sh   | 16 | 64G | minimap2/2.29, SAMtools/1.21 |
| submit_ont.sh      | 16 | 128G | minimap2/2.29, SAMtools/1.21 |

---

## Step 9 — Verify

After all files are written:

```bash
find /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations/src -name "*.py" | sort
ls /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations/*.py \
   /gpfs/commons/home/atalukder/RNA_Splicing/code/Simulations/*.sh 2>/dev/null
```

Confirm every file in the target directory structure exists. Report any still missing.

---

## Reminders

- Every file: top-level docstring, per-function docstrings, print checkpoints
- CONFIG in entry points only — never scatter paths through script bodies
- `subprocess.run(..., check=True)` for all external tool calls
- Do NOT modify `lrgsp_simulation` env — flag missing packages to user instead
- Do NOT touch `data/Shree_stuff/Simulation/` — read-only reference only
