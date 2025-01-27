#!/bin/bash
#SBATCH --job-name=h1_pb_aln                        # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=20G                                   # Job memory request
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load minimap2/2.17
module load samtools/1.9

minimap2 -ax splice /gpfs/commons/home/sraghavendra/SOTA/h1_cell/ref/lrgasp_grch38_sirvs.fasta /gpfs/commons/home/sraghavendra/SOTA/h1_cell/ont.fastq > /gpfs/commons/home/sraghavendra/SOTA/h1_cell/aln_ont_long.sam
samtools view -bS /gpfs/commons/home/sraghavendra/SOTA/h1_cell/aln_ont_long.sam > /gpfs/commons/home/sraghavendra/SOTA/h1_cell/aln_ont_long.bam
samtools flagstat /gpfs/commons/home/sraghavendra/SOTA/h1_cell/aln_ont_long.bam > /gpfs/commons/home/sraghavendra/SOTA/h1_cell/report_ont_long.txt

