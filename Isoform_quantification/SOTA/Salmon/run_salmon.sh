#!/bin/bash
#SBATCH --job-name=run_salmon                     # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=50G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load salmon

salmon quant -t /gpfs/commons/home/sraghavendra/PacBio/reference/sota/transcriptome.fna -l A -a /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transcriptome/aln_52_shortAligned.sortedByCoord.out.bam -o salmon_quant52
