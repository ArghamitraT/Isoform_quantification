#!/bin/bash
#SBATCH --job-name=run_mpaqt                        # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate mpaqt

mpaqt quant \
        --project nc \
        --sample day0_rep1

mpaqt quant \
        --project nc \
        --sample day0_rep2

mpaqt quant \
        --project kallisto \
        --sample day0_rep1

mpaqt quant \
        --project kallisto \
        --sample day0_rep2