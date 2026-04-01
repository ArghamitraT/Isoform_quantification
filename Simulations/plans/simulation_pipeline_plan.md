# Simulation Pipeline — Implementation Plan

**Goal:** Clean re-implementation of the RNA-seq read simulation pipeline in `code/Simulations/`, supporting synthetic abundance generation and all three sequencing technologies.

**Reference:** Original pipeline documented in `code/Simulations/REFERENCE.md`.

---

## Environment

**Conda env:** `lrgsp_simulation` — do not modify it.

**Packages already available in `lrgsp_simulation`:**
- numpy, pandas, scipy, pysam, biopython
- matplotlib, seaborn, tqdm
- nanocount, torch, pyro-ppl

**Packages NOT in env — load via `module load` in SLURM scripts:**
| Tool | Module |
|------|--------|
| minimap2 | `minimap2/2.29` |
| RSEM | `RSEM/1.3.3-foss-2022a` |
| STAR | `STAR/2.7.11b-GCC-13.2.0` |
| SAMtools | `SAMtools/1.21` |

**⚠ Missing: `gffutils`** — not in `NanoCount_5_shree`. GTF parsing will use Biopython + pandas instead (no gffutils dependency).

**External simulator tools (read-only, do not modify):**
- NanoSim: `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation/lrgasp-simulation/src/NanoSim/`
- IsoSeqSim: `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation/lrgasp-simulation/src/isoseqsim/`

---

## Directory Structure (target — all inside `code/Simulations/`)

```
code/Simulations/
├── plans/
│   └── simulation_pipeline_plan.md   # this file
├── commands/
│   ├── simulation.md                  # Claude skill for guided workflow
│   └── install_pipeline.md            # Claude skill for installing all pipeline files
├── src/
│   ├── prepare_reference.py           # Phase 1: GTF + genome → simulation-ready refs
│   ├── generate_abundances.py         # Phase 2: synthetic TPM profile generation
│   ├── simulate_illumina.py           # Phase 3a: Illumina (wraps RSEM)
│   ├── simulate_pacbio.py             # Phase 3b: PacBio (wraps IsoSeqSim)
│   ├── simulate_ont.py                # Phase 3c: ONT (wraps Trans-NanoSim)
│   └── utility.py                     # Shared: run dir, runtime, code snapshot
├── main.py                            # Entry point with CONFIG section
├── submit_illumina.sh                 # SLURM job: Illumina only
├── submit_pacbio.sh                   # SLURM job: PacBio only
├── submit_ont.sh                      # SLURM job: ONT only
└── REFERENCE.md                       # Reference documentation (already created)
```

---

## Pipeline Phases

### Phase 1 — Reference Preparation (`src/prepare_reference.py`)

**Input:**
- Reference genome FASTA (human GRCh38)
- GTF annotation (e.g. GENCODE)
- Optional: SQANTI output for novel isoforms

**Steps:**
1. Parse GTF using pandas (no gffutils) → build transcript-to-exon table
2. Extract transcript sequences from genome using pysam/biopython
3. Add 100 bp poly-A tails to all transcripts
4. Optionally insert novel isoforms from SQANTI (skip `not_in_catalog`, `novel_in_catalog`)

**Output:** `<prefix>.annotation.gtf`, `<prefix>.transcripts.fasta`, `<prefix>.genome.fasta`, `<prefix>.novel_isoforms.tsv`

---

### Phase 2 — Synthetic Abundance Generation (`src/generate_abundances.py`)

Three modes (selected via CONFIG):

#### Mode A — Log-normal (realistic)
- Draw TPM values from a log-normal distribution
- Parameters: `mu` and `sigma` calibrated to match real RNA-seq data
- Applies a dropout fraction (fraction of transcripts with TPM = 0)
- Produces a realistic long tail: few highly expressed, many lowly expressed

#### Mode B — Custom / Cancer
- User provides a custom TPM table (TSV: `transcript_id`, `TPM`)
- Optionally: up-regulate a specified gene set (e.g. known oncogenes) by a fold-change factor
- Supports mixing: e.g. 80% from custom table + 20% log-normal noise
- Useful for simulating tumor-specific expression profiles

#### Mode C — Uniform (benchmarking)
- All transcripts assigned equal TPM
- Useful as a null/control condition

**Output:** `<prefix>.abundances.tsv` (columns: `transcript_id`, `TPM`)

---

### Phase 3 — Read Simulation (separate scripts per technology)

Each simulator is independent and submitted as its own SLURM job.

#### Phase 3a — Illumina (`src/simulate_illumina.py`)
- Wraps **RSEM** (loaded via `module load RSEM/1.3.3-foss-2022a`)
- Uses paired-end Illumina 150 bp model
- Optional noise fraction (default 10%)
- **Input:** reference prefix + abundances TSV
- **Output:** `Illumina.simulated_1.fq`, `Illumina.simulated_2.fq`, isoform results

#### Phase 3b — PacBio (`src/simulate_pacbio.py`)
- Wraps **IsoSeqSim** (called from its read-only path in Shree's folder)
- PacBio Sequel CCS error profile (sub 0.4%, ins 0.6%, del 0.6%)
- **Input:** reference prefix + abundances TSV
- **Output:** `PacBio.simulated.fasta`, `*.isoform_counts.tsv`, `*.read_to_isoform.tsv`

#### Phase 3c — ONT (`src/simulate_ont.py`)
- Wraps **Trans-NanoSim** (called from its read-only path in Shree's folder)
- Supports both dRNA and cDNA chemistries
- Optional noise/unaligned reads
- **Input:** reference prefix + abundances TSV
- **Output:** `ONT.simulated.fastq`, `*.isoform_counts.tsv`, `*.read_to_isoform.tsv`

---

## SLURM Job Strategy

Each technology runs as a **separate SLURM job** for flexibility (re-run one without re-running others):

| Script | Technology | Resources | Modules to load |
|--------|-----------|-----------|----------------|
| `submit_illumina.sh` | Illumina | 16 CPUs, 64G RAM | RSEM, SAMtools |
| `submit_pacbio.sh` | PacBio | 16 CPUs, 64G RAM | minimap2, SAMtools |
| `submit_ont.sh` | ONT | 16 CPUs, 128G RAM | minimap2, SAMtools |

Submit independently:
```bash
sbatch submit_illumina.sh
sbatch submit_pacbio.sh
sbatch submit_ont.sh
```

Each SLURM script activates `lrgsp_simulation` and loads the relevant modules.

---

## Output Layout (per run)

All results follow the project standard timestamped folder format under `files/results/`:

```
files/results/exprmnt_{timestamp}/
├── experiment_description.log   # config dump + what/why/expected outcome
├── running.log                  # stdout (first line: "Script: <path>")
├── abundances.tsv               # synthetic TPM profile used
├── Illumina.simulated_1.fq      # Illumina reads
├── Illumina.simulated_2.fq
├── PacBio.simulated.fasta       # PacBio reads
├── ONT.simulated.fastq          # ONT reads
├── *.isoform_counts.tsv         # per-technology isoform counts
├── *.read_to_isoform.tsv        # read → isoform mapping
├── results_summary.txt          # read counts, isoform coverage stats
├── runtime.txt                  # total wall-clock time
└── code_snapshot/               # exact copy of code/Simulations/ at run time
```

---

## CONFIG Section (in `main.py`)

```python
# ============================================================
# CONFIG — edit these variables before running; do not edit below
# ============================================================
reference_genome   = "/path/to/GRCh38.fa"
reference_gtf      = "/path/to/gencode.vXX.annotation.gtf"
output_prefix      = "human_sim"

# Abundance mode: "lognormal", "custom", or "uniform"
abundance_mode     = "lognormal"
custom_tpm_file    = ""              # path to TSV if mode="custom"
cancer_genes_file  = ""              # optional: genes to up-regulate (mode="custom")
cancer_fold_change = 5.0             # fold-change for cancer genes
dropout_fraction   = 0.3             # fraction of transcripts with TPM=0 (lognormal only)

# Simulation read counts
illumina_count     = 100_000_000     # read pairs
pb_count           = 10_000_000      # reads
ont_count          = 30_000_000      # reads
ont_type           = "cDNA"          # "cDNA" or "dRNA"

# Run settings
threads            = 16
seed               = 42
noise_reads        = False
keep_isoform_ids   = False
test_mode          = False           # True = low read counts for quick testing
# ============================================================
```

---

## Implementation Order

1. [ ] `src/utility.py` — run dir, runtime, code snapshot helpers
2. [ ] `src/prepare_reference.py` — Phase 1 (GTF parsing via pandas, no gffutils)
3. [ ] `src/generate_abundances.py` — Phase 2 (all three modes)
4. [ ] `src/simulate_illumina.py` — Phase 3a (wraps RSEM module)
5. [ ] `src/simulate_pacbio.py` — Phase 3b (wraps IsoSeqSim from Shree's read-only path)
6. [ ] `src/simulate_ont.py` — Phase 3c (wraps NanoSim from Shree's read-only path)
7. [ ] `main.py` — entry point tying all phases together with CONFIG
8. [ ] `submit_illumina.sh`, `submit_pacbio.sh`, `submit_ont.sh` — SLURM wrappers
9. [ ] Test each phase independently before wiring together
