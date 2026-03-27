# Simulations — Reference Documentation

> **DO NOT TOUCH:** `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation`
> This file documents it for reference only, to inform building the clean re-implementation here.

---

## Purpose

A 3-phase RNA-seq read simulation pipeline that generates synthetic sequencing reads for three technologies based on real reference genomes and expression profiles:
- **Illumina** (short-read, paired-end)
- **PacBio CCS** (long-read)
- **Oxford Nanopore (ONT)** (long-read, dRNA or cDNA)

---

## Original Directory Structure (`data/Shree_stuff/Simulation/`)

```
data/Shree_stuff/Simulation/
├── lrgasp-simulation/              # Main simulation pipeline (LRGASP project)
│   ├── simulate.py                 # Master simulation orchestrator
│   ├── prepare_reference_data.py   # GTF + FASTA reference preparation
│   ├── quantify.py                 # Transcript abundance from real reads
│   ├── run_sim.sh                  # Shell wrapper for full pipeline
│   ├── run_sim_job.sh              # SLURM job submission script
│   ├── run_simulate_test.sh        # Quick test run script
│   ├── src/
│   │   ├── simulate_illumina.py    # Illumina read simulator (wraps RSEM)
│   │   ├── simulate_pacbio.py      # PacBio read simulator (wraps IsoSeqSim)
│   │   ├── simulate_ont.py         # ONT read simulator (wraps Trans-NanoSim)
│   │   ├── isoseqsim/              # IsoSeqSim tool source
│   │   ├── NanoSim/                # Trans-NanoSim tool source
│   │   └── RSEM/                   # RSEM tool source
│   ├── data/
│   │   ├── test_data/              # Minimal test dataset
│   │   └── human_chr22/            # Human chr22 example data + SQANTI output
│   ├── sim_result*/                # Full simulation result folders
│   ├── similation_test*/           # Test run folders (note typo: "similation")
│   └── README.md                   # LRGASP pipeline documentation
│
├── STAR_gen_inds.sh                # STAR genome index generation
├── run_lrgasp.sh                   # LRGASP runner (hardcoded to Shree's paths)
├── run_aln.sh                      # PacBio alignment with minimap2
├── run_aln_STAR.sh                 # Illumina alignment with STAR
├── alignments/                     # Alignment output BAM files
├── reference/                      # STAR reference genome indices
│
├── sota/                           # State-of-the-art quantification tools
│   ├── lr_kallisto/                # Kallisto pseudo-alignment (~64 GB)
│   ├── nanocount/                  # Nanocount quantification
│   ├── oarfish/                    # Oarfish quantification
│   └── mpaqt/                      # MPAQT quantification
│
├── common_transcripts.py           # Find common isoforms across two methods
├── compare_ill_pb.py               # Correlate Illumina vs PacBio quantification
├── compare_thetas.py               # Compare estimated parameters
├── gen_tpm.py                      # Convert read counts to TPM
├── gen_read_counts.py              # Generate read count tables
├── get_stats.py                    # Basic statistics (isoform counts, reads)
│
├── PB_run*.csv / PB_*_tpms.csv    # PacBio quantification result CSVs
├── ONT_df.csv                      # ONT result dataframe
└── stdout_*.log                    # SLURM job output logs
```

---

## Pipeline Flow (3 Phases)

### Phase 1 — Reference Preparation (`prepare_reference_data.py`)
- Loads GTF annotation into a gffutils database
- Optionally integrates novel isoforms from SQANTI output
  - Skips "not_in_catalog" and "novel_in_catalog" isoform types
- Adds 100 bp poly-A tails to all transcript sequences
- **Outputs:** `<prefix>.annotation.gtf`, `<prefix>.transcripts.fasta`, `<prefix>.genome.fasta`, `<prefix>.novel_isoforms.tsv`

### Phase 2 — Abundance Quantification (`quantify.py`)
- Maps real long reads to reference transcripts using **minimap2**
- Counts primary alignments per transcript (excludes secondary/supplementary)
- Optionally assigns random high counts to mandatory (novel) isoforms
- Calculates TPM
- **Output:** TSV with columns `transcript_id`, `counts`, `TPM`

### Phase 3 — Read Simulation (`simulate.py`)
- Seeds all RNGs for reproducibility
- Calls three sub-simulators:

  | Branch | External Tool | Wrapper Script | Outputs |
  |--------|--------------|----------------|---------|
  | Illumina | RSEM | `simulate_illumina.py` | `Illumina.simulated_1.fq`, `_2.fq` |
  | PacBio | IsoSeqSim | `simulate_pacbio.py` | `PacBio.simulated.fasta` |
  | ONT | Trans-NanoSim | `simulate_ont.py` | `ONT.simulated.fastq` |

- Each branch also produces `*.isoform_counts.tsv` and `*.read_to_isoform.tsv`

### Post-Simulation — Alignment & Analysis
- `run_aln.sh`: minimap2 (`-ax splice -uf -C5`) → BAM (PacBio)
- `run_aln_STAR.sh`: STAR → BAM (Illumina paired-end)
- Root-level Python scripts compare quantification outputs across methods

---

## Key Parameters (`simulate.py`)

| Flag | Default | Meaning |
|------|---------|---------|
| `--reference_prefix` | — | Path prefix to Phase 1 outputs |
| `--counts` | — | TSV from Phase 2 (quantify.py output) |
| `--output` | — | Output directory |
| `--threads` `-t` | 16 | CPU threads |
| `--seed` `-s` | 11 | Random seed for reproducibility |
| `--illumina_count` | 100,000,000 | Illumina read pairs to simulate |
| `--pb_count` | 10,000,000 | PacBio reads to simulate |
| `--ont_count` | 30,000,000 | ONT reads to simulate |
| `--ont_type` | `cDNA` | `dRNA` or `cDNA` chemistry |
| `--noise_reads` | off | Include background noise reads |
| `--keep_isoform_ids` | off | Preserve isoform IDs in read names |
| `--test_mode` | off | Low read counts for quick testing |

---

## How to Run (original pipeline, for reference)

```bash
cd /gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation/lrgasp-simulation

# Phase 1: Prepare reference
python prepare_reference_data.py \
  --reference_annotation data/human_chr22/gencode.v36.annotation.chr22.gtf \
  --reference_transcripts data/human_chr22/gencode.v36.transcripts.chr22.fa \
  --reference_genome data/human_chr22/GRCh38.chr22.fa \
  --sqanti_prefix data/human_chr22/rat_human_chr22 \
  --n_random_isoforms 50 \
  --output reference_data_chr22/human.chr22

# Phase 2: Quantify transcripts from real reads
python quantify.py \
  --fastq data/human_chr22/Human.PacBio.ENCFF.chr22.fq \
  -t 16 \
  --reference_transcripts reference_data_chr22/human.chr22.transcripts.fasta \
  --mandatory reference_data_chr22/human.chr22.novel_isoforms.tsv \
  --output reference_data_chr22/human.chr22.counts.tsv

# Phase 3: Simulate reads
python simulate.py \
  --reference_prefix reference_data_chr22/human.chr22 \
  --counts reference_data_chr22/human.chr22.counts.tsv \
  -t 16 -s 22 \
  --output chr22_simulated/
```

---

## Known Issues in Shree's Code (do NOT replicate here)

1. **Hardcoded paths** — All shell scripts reference `/gpfs/commons/home/sraghavendra/...`; fail for any other user
2. **Hardcoded email** — SLURM scripts mail to `sraghavendra@nygenome.org`
3. **Typo in result dir names** — `similation_test*` (missing "u" in "simulation")
4. **Stale backup file** — `src/simulate_pacbio copy.py` (unused copy)
5. **Commented-out code** — `run_sim.sh` has ~30% commented lines with no explanation
6. **No CONFIG section** — Parameters scattered throughout scripts, not centralized
7. **No docstrings** — No top-level or function-level documentation
8. **No runtime logging or code snapshots** — Cannot reproduce exact runs
9. **Large output files in repo** — `ONT.simulated.read_to_isoform.tsv` (1.3 GB), `PacBio.simulated.read_to_isoform.tsv` (487 MB)
10. **No versioning metadata** — Multiple result folders with no record of what parameters produced them

---

## Goals for Clean Re-implementation Here

- Single CONFIG section at top of every entry point
- Top-level and per-function docstrings (description, Args, Returns)
- Modular structure split by responsibility (one file per simulator, shared utilities)
- Checkpoints/print statements at each stage
- Runtime logging (`runtime.txt`), code snapshots (`code_snapshot/`), timestamped result folders
- No hardcoded paths or user-specific values
- Conda environment in `env/` (e.g. `sim_env.yml`)
