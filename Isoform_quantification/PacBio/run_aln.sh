#!/bin/bash
#SBATCH --job-name=ill_alns_trans                    # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=20G                                   # Job memory request
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load minimap2/2.17
module load samtools/1.9

minimap2 -ax sr -t4 /gpfs/commons/home/sraghavendra/PacBio/reference/transcriptome.fna /gpfs/commons/home/sraghavendra/PacBio/reads/short/day0_rep2_R1.fastq /gpfs/commons/home/sraghavendra/PacBio/reads/short/day0_rep2_R2.fastq > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02_short.sam
#minimap2 -ax sr -t4 /gpfs/commons/home/sraghavendra/PacBio/reference/GRCh38_latest_genomic.fna /gpfs/commons/home/sraghavendra/PacBio/reads/short/transc_minimap/day0_rep2_R1.fastq.gz /gpfs/commons/home/sraghavendra/PacBio/reads/short/transc_minimap/day0_rep2_R2.fastq.gz > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02.sam
samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02_short.sam > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/report_02_short.txt
samtools view -bS /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02_short.sam > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02_short.bam

#minimap2 -ax splice -uf -C5 /gpfs/commons/home/sraghavendra/PacBio/reference/transcriptome.fna /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_02.fastq > /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/aln_02_long.sam
#samtools view -bS /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/aln_02_long.sam > /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/aln_02_long.bam
#samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/aln_02_long.bam > /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/report_02_long.txt

