#!/bin/bash
#SBATCH --job-name=SIRV_aln_E0_ST                   # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=4G                                    # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --time=4:00:00                              # Time limit 4 hours
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load minimap2/2.17
module load samtools/1.9

#minimap2 -ax splice --splice-flank=no -k14  /gpfs/commons/home/sraghavendra/SIRV/reference/sirv_transcriptome.fa /gpfs/commons/home/sraghavendra/SIRV/reads/SRR6058584_E0.fastq > /gpfs/commons/home/sraghavendra/SIRV/alignments/aln_E0.sam

#minimap2 -ax splice --splice-flank=no -k14  /gpfs/commons/home/sraghavendra/SIRV/reference/sirv_transcriptome.fa /gpfs/commons/home/sraghavendra/SIRV/reads/SRR6058583_E2.fastq > /gpfs/commons/home/sraghavendra/SIRV/alignments/aln_E2.sam

#samtools flagstat /gpfs/commons/home/sraghavendra/SIRV/alignments/aln_E2.sam > /gpfs/commons/home/sraghavendra/SIRV/alignments/report_E2.txt

samtools view -bS /gpfs/commons/home/sraghavendra/SIRV/alignments/aln_E0.sam > /gpfs/commons/home/sraghavendra/SIRV/alignments/aln_E0.bam

samtools view -bS /gpfs/commons/home/sraghavendra/SIRV/alignments/aln_E2.sam > /gpfs/commons/home/sraghavendra/SIRV/alignments/aln_E2.bam
