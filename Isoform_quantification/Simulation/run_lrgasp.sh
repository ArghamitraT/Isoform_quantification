#!/bin/bash
#SBATCH --job-name=run_lrgasp_sim                  # Set the job name
#SBATCH --partition=pe2                            # Partition Name
#SBATCH --mail-type=END,FAIL                       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=50G                                  # Request 20G
#SBATCH --mail-user=sraghavendra@nygenome.org      # Where to send mail
#SBATCH --output=stdout_%j.log                     # Standard output and error log

/gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/run_simulate_test.sh
