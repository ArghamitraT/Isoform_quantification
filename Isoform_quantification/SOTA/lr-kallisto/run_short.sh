#!/bin/bash
#SBATCH --job-name=short_kallisto_gen               # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

#cd /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT
#~/miniconda3/bin/conda init bash
#source /gpfs/commons/home/sraghavendra/.bashrc
#conda activate mpaqt

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/gen_index \
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_01 \
#  /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_01_R1.fastq /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_01_R2.fastq

#  /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/gen_index \
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_02 \
#  /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_02_R1.fastq /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_02_R2.fastq

#  /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/gen_index \
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_51 \
#  /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_51_R1.fastq /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_51_R2.fastq

#  /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/gen_index \
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_52 \
#  /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_52_R1.fastq /gpfs/commons/home/sraghavendra/PacBio/reads/short/ill_52_R2.fastq

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/sim_gen_index \
 -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_sim1 \
 /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/human_simulated_job_correct/Illumina.simulated_1.fq \
 /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/human_simulated_job_correct/Illumina.simulated_2.fq

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/sim_gen_index \
 -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_sim2 \
 /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/Illumina.simulated_1.fq \
 /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/Illumina.simulated_2.fq
