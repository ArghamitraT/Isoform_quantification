#!/bin/bash
#SBATCH --job-name=flnc_bam_to_fq                   # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=8G                                    # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=4:00:00                              # Time limit 4 hours
#SBATCH --output=stdout_%j.log                      # Standard output and error log

#wget https://downloads.pacbcloud.com/public/dataset/MAS-Seq-bulk-2023-WTC11/day0-rep1/2-FLNC/flnc.bam

#wget https://downloads.pacbcloud.com/public/dataset/MAS-Seq-bulk-2023-WTC11/day0-rep2/2-FLNC/flnc.bam

module load samtools/1.9

#samtools view /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_01.bam | less > /gpfs/commons/home/sraghavendra/PacBio/reads/long/view.txt

#samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_01.bam > /gpfs/commons/home/sraghavendra/PacBio/reads/long/report_01.txt
#samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_02.bam	> /gpfs/commons/home/sraghavendra/PacBio/reads/long/report_02.txt

samtools fastq /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_01.bam > /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_01.fastq
samtools fastq /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_02.bam > /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_02.fastq
