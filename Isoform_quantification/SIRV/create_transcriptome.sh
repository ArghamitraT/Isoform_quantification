#!/bin/bash
#SBATCH --job-name=create_transcriptome             # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load gffread/0.12.8

gffread -w /gpfs/commons/home/sraghavendra/PacBio/reference/transcriptome.fa -g /gpfs/commons/home/sraghavendra/PacBio/reference/ncbi/GCF_000001405.40_GRCh38.p14_genomic.fna /gpfs/commons/home/sraghavendra/PacBio/reference/ncbi/genomic.gtf
