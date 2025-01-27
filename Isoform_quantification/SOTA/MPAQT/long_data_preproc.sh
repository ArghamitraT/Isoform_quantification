#!/bin/bash
#SBATCH --job-name=mpaqt_long_data_preproc         # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                   # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate mpaqt

mpaqt prepare long-read \
        --project nc \
        --sample day0_rep1 \
        --counts /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/lr/tx_counts_NC_for_MPAQT_01.tsv

mpaqt prepare long-read \
        --project nc \
        --sample day0_rep2 \
        --counts /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/lr/tx_counts_NC_for_MPAQT_01.tsv

mpaqt prepare long-read \
        --project kallisto \
        --sample day0_rep1 \
        --counts /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/lr/tx_counts_Kallisto_for_MPAQT_01.tsv

mpaqt prepare long-read \
        --project kallisto \
        --sample day0_rep2 \
        --counts /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT/data/lr/tx_counts_Kallisto_for_MPAQT_01.tsv