#!/bin/bash
#SBATCH --job-name=isoquant_run                     # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

#cd /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT
~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate isoquant

isoquant.py --reference /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna \
--genedb /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf \
--fastq /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_02.fastq \
--data_type pacbio_ccs -o /gpfs/commons/home/sraghavendra/SOTA/IsoQuant/output/

