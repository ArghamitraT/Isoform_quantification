#!/bin/bash
#SBATCH --job-name=downsample_01_1                  # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=40G                                   # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load minimap2/2.17
module load samtools/1.9

samtools view -b -s 0.35 /gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/aln_01.bam > /gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/downsampled/aln_01_1.bam