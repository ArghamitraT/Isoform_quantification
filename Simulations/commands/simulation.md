# /simulation — Guided RNA-seq Simulation Workflow

You are helping the user run the RNA-seq simulation pipeline in `code/Simulations/`.
Always read `code/Simulations/plans/simulation_pipeline_plan.md` and `code/Simulations/REFERENCE.md` before doing anything.

**NEVER touch `/gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff/Simulation`.**
**NEVER modify the `NanoCount_5_shree` conda environment.**

---

## Environment Facts (verified)

- **Conda env:** `NanoCount_5_shree`
- **Packages in env:** numpy, pandas, scipy, pysam, biopython, matplotlib, seaborn, tqdm, nanocount, torch, pyro-ppl
- **Missing from env — load as modules in SLURM scripts:**
  - `minimap2/2.29`
  - `RSEM/1.3.3-foss-2022a`
  - `STAR/2.7.11b-GCC-13.2.0`
  - `SAMtools/1.21`
- **No gffutils** — GTF parsing uses pandas instead
- **NanoSim path (read-only):** `.../data/Shree_stuff/Simulation/lrgasp-simulation/src/NanoSim/`
- **IsoSeqSim path (read-only):** `.../data/Shree_stuff/Simulation/lrgasp-simulation/src/isoseqsim/`

---

## What This Skill Does

Guides the user through a complete RNA-seq simulation workflow:

1. **Setup** — validate paths, CONFIG, and environment
2. **Reference prep** — Phase 1: GTF + genome → simulation-ready references
3. **Abundance generation** — Phase 2: synthetic TPM profile (log-normal, custom/cancer, or uniform)
4. **Read simulation** — Phase 3: Illumina, PacBio, ONT as separate SLURM jobs
5. **Analysis** — check outputs, basic QC, isoform coverage stats

---

## Step 1 — Setup & Validation

When invoked:

1. Read `code/Simulations/main.py` CONFIG section. If the file does not exist yet, tell the user to run the implementation steps in the plan first.

2. Verify key inputs exist:
   ```bash
   ls $reference_genome $reference_gtf
   conda env list | grep NanoCount_5_shree
   module avail 2>&1 | grep -E "RSEM|minimap2|STAR|SAMtools"
   ```

3. Report any missing files or modules before proceeding.

---

## Step 2 — Reference Preparation

```bash
conda activate NanoCount_5_shree
python code/Simulations/src/prepare_reference.py \
    --genome $reference_genome \
    --gtf $reference_gtf \
    --output $output_prefix
```

Expected outputs: `<prefix>.annotation.gtf`, `<prefix>.transcripts.fasta`, `<prefix>.genome.fasta`

Verify all three exist before moving to Phase 2.

---

## Step 3 — Synthetic Abundance Generation

Ask the user which mode if not set in CONFIG, then run:

### Log-normal (realistic)
```bash
python code/Simulations/src/generate_abundances.py \
    --mode lognormal \
    --transcripts <prefix>.transcripts.fasta \
    --dropout $dropout_fraction \
    --seed $seed \
    --output <prefix>.abundances.tsv
```

### Custom / Cancer
```bash
python code/Simulations/src/generate_abundances.py \
    --mode custom \
    --transcripts <prefix>.transcripts.fasta \
    --custom_tpm $custom_tpm_file \
    --cancer_genes $cancer_genes_file \
    --fold_change $cancer_fold_change \
    --seed $seed \
    --output <prefix>.abundances.tsv
```
If `cancer_genes_file` is empty, ask the user before proceeding.

### Uniform (benchmarking)
```bash
python code/Simulations/src/generate_abundances.py \
    --mode uniform \
    --transcripts <prefix>.transcripts.fasta \
    --output <prefix>.abundances.tsv
```

After generation, print: total transcripts, expressed (TPM > 0), top 10 by TPM.

---

## Step 4 — Read Simulation (Separate SLURM Jobs)

Each technology is a separate SLURM job. Show the user all three commands:

```bash
sbatch code/Simulations/submit_illumina.sh   # 16 CPUs, 64G  — loads RSEM, SAMtools
sbatch code/Simulations/submit_pacbio.sh     # 16 CPUs, 64G  — loads minimap2, SAMtools
sbatch code/Simulations/submit_ont.sh        # 16 CPUs, 128G — loads minimap2, SAMtools
```

Ask which to submit. Only submit what the user confirms.

After submission, show job IDs and how to monitor:
```bash
squeue -u $USER
```

---

## Step 5 — Analysis & QC

Once jobs complete:
```bash
conda activate NanoCount_5_shree
python code/Simulations/src/analyze_outputs.py --run_dir $run_dir
```

Report:
- Total reads simulated per technology
- Isoform coverage: % of expressed transcripts with ≥1 read
- Top 20 most simulated isoforms vs input TPM (spot-check correlation)
- Isoforms with TPM > 0 but 0 simulated reads (flag as potential issues)

---

## Reminders for Every Run

- All scripts live in `code/Simulations/` — nothing outside this folder
- CONFIG lives only in `main.py` — never hardcode paths elsewhere
- Every run saves to `files/results/exprmnt_{timestamp}/`
- `runtime.txt` and `code_snapshot/` are mandatory in every result folder
- First line of `running.log` must be `Script: <path/to/script>`
- Do not modify `NanoCount_5_shree` — if a package is missing, flag it to the user
