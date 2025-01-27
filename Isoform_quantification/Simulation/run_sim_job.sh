#!/bin/bash
#SBATCH --job-name=run_sim_2                        # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate mpaqt

python /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/simulate.py \
  --reference_prefix /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/ref_data/human \
  --counts /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/ref_data/human.counts.tsv \
  -t 16 -s 22 \
  --illumina_count 10000000 \
  --pb_count 1000000 \
  --output /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/
