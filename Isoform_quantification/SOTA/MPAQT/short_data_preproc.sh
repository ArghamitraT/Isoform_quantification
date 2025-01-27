#!/bin/bash
#SBATCH --job-name=sim1_mpaqt_short_data_preproc      # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

#cd /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT
~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate mpaqt


/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto bus \
	--index /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/references/gencode.v46.index.kallisto \
	--num \
	--paired \
	--technology bulk \
	--output-dir /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/nc/day0_rep1 \
	/gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_01_R1.fastq \
 /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_01_R2.fastq

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto bus \
	--index /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/references/gencode.v46.index.kallisto \
	--num \
	--paired \
	--technology bulk \
	--output-dir /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/nc/day0_rep2 \
	/gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_02_R1.fastq \
 /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_02_R2.fastq

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto bus \
	--index /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/references/gencode.v46.index.kallisto \
	--num \
	--paired \
	--technology bulk \
	--output-dir /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/kallisto/day0_rep1 \
	/gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_01_R1.fastq \
 /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_01_R2.fastq

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto bus \
	--index /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/references/gencode.v46.index.kallisto \
	--num \
	--paired \
	--technology bulk \
	--output-dir /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/kallisto/day0_rep2 \
	/gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_02_R1.fastq \
 /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_02_R2.fastq

# kallisto bus \
# --index /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/references/gencode.v46.index.kallisto \
# --num \
# --paired \
# --technology bulk \
# --output-dir /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/sim1 \
# /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/human_simulated_job_correct/Illumina.simulated_1.fq \
# /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/human_simulated_job_correct/Illumina.simulated_2.fq

# kallisto bus \
# --index /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/references/gencode.v46.index.kallisto \
# --num \
# --paired \
# --technology bulk \
# --output-dir /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/sanity_check/day0_rep1 \
# /gpfs/commons/home/sraghavendra/PacBio/reads/short/downsampled/ds_01_R1.fastq \
# /gpfs/commons/home/sraghavendra/PacBio/reads/short/downsampled/ds_01_R2.fastq

# kallisto bus --index /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/references/gencode.v46.index.kallisto \
# --num \
# --paired \
# --technology bulk \
# --output-dir /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/d5r1 \
# /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_51_R1.fastq /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_51_R2.fastq

mpaqt prepare short-read \
	--project nc \
	--sample day0_rep1 \
	--bus /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/nc/day0_rep1/output.bus \
	--matrix_ec /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/nc/day0_rep1/matrix.ec

mpaqt prepare short-read \
	--project nc \
	--sample day0_rep2 \
	--bus /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/nc/day0_rep2/output.bus \
	--matrix_ec /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/nc/day0_rep2/matrix.ec

mpaqt prepare short-read \
	--project kallisto \
	--sample day0_rep1 \
	--bus /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/kallisto/day0_rep1/output.bus \
	--matrix_ec /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/kallisto/day0_rep1/matrix.ec

mpaqt prepare short-read \
	--project kallisto \
	--sample day0_rep2 \
	--bus /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/kallisto/day0_rep2/output.bus \
	--matrix_ec /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/bus/kallisto/day0_rep2/matrix.ec
