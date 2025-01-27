#!/bin/bash
#SBATCH --job-name=build_salmon                     # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

#cd /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT
~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate mpaqt

make -C /gpfs/commons/home/sraghavendra/SOTA/salmon/salmon/build/
