#!/bin/bash
#SBATCH --job-name=lr_kallisto_gen                  # Job name
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
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_01 --long /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_01.fastq

#  /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/gen_index \
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_02 --long /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_02.fastq

#  /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/gen_index \
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_51 --long /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_51.fastq

#  /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/gen_index \
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_52 --long /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_52.fastq

 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/sim_gen_index \
 -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_sim1 --long /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/human_simulated_job_correct/PacBio.simulated.fasta

  /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/kallisto quant -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/sim_gen_index \
 -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_sim2 --long /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/PacBio.simulated.fasta