#!/bin/bash
#SBATCH --job-name=kallisto_index                      # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                   # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

#cd /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT
#~/miniconda3/bin/conda init bash
#source /gpfs/commons/home/sraghavendra/.bashrc
#conda activate mpaqt

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto index \
# 	--index /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/index \
# 	/gpfs/commons/home/sraghavendra/PacBio/reference/sota/transcriptome.fna

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto index \
# 	--index /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/sim_index \
# 	/gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/ref_data/human.transcripts.fasta

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto index \
# 	--index /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/gen_index \
# 	/gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto index \
	--index /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/sim_gen_index \
	/gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/ref_data/human.genome.fasta
