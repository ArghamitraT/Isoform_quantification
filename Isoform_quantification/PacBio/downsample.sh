#!/bin/bash
#SBATCH --job-name=10_2_downsample_long                # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=100G                                  # Job memory request
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load minimap2/2.17
module load samtools/1.9


samtools view -b -s 1.10 /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/orig/aln_01_long.bam > /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/ds_10_num2_aln_01_long.bam

samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/ds_10_num2_aln_01_long.bam > /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/report_ds_10_num2_01.txt

samtools view -b -s 1.10 /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/orig/aln_02_long.bam > /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/ds_10_num2_aln_02_long.bam

samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/ds_10_num2_aln_02_long.bam > /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/report_ds_10_num2_02.txt
