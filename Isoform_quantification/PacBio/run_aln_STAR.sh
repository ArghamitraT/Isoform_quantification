#!/bin/bash
#SBATCH --job-name=STAR_trans_aln_02                # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=100G                                   # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load star/2.5.2a
module load samtools/1.9

STAR --genomeDir /gpfs/commons/home/sraghavendra/PacBio/reference/transc/ \
--runThreadN 6 \
--readFilesIn /gpfs/commons/home/sraghavendra/PacBio/reads/short/day0_rep2_R1.fastq /gpfs/commons/home/sraghavendra/PacBio/reads/short/day0_rep2_R2.fastq \
--outFileNamePrefix /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transcriptome/aln_02_short_ \
--outSAMtype BAM SortedByCoordinate \

samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transcriptome/aln_01_short_Aligned.sortedByCoord.out.bam > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transcriptome/report_01_short.txt
samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transcriptome/aln_02_short_Aligned.sortedByCoord.out.bam > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transcriptome/report_02_short.txt

#--quantMode TranscriptomeSAM \

