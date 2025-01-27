#!/bin/bash
#SBATCH --job-name=h1_ill_aln                      # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=100G                                   # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load star/2.5.2a
module load samtools/1.9

STAR --genomeDir /gpfs/commons/home/sraghavendra/SOTA/h1_cell/ref/ \
--runThreadN 6 \
--readFilesIn /gpfs/commons/home/sraghavendra/SOTA/h1_cell/illumina.fastq \
--outFileNamePrefix /gpfs/commons/home/sraghavendra/SOTA/h1_cell/aln_ill_short_ \
--outSAMtype BAM SortedByCoordinate \

samtools flagstat /gpfs/commons/home/sraghavendra/SOTA/h1_cell/aln_ill_short_Aligned.sortedByCoord.out.bam > /gpfs/commons/home/sraghavendra/SOTA/h1_cell/report_ill_short.txt
