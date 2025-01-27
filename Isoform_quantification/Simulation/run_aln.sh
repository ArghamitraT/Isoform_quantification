#!/bin/bash
#SBATCH --job-name=pb_aln_sim                       # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=20G                                   # Job memory request
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load minimap2/2.17
module load samtools/1.9

# minimap2 -ax sr -t4 /gpfs/commons/home/sraghavendra/PacBio/reference/transcriptome.fna /gpfs/commons/home/sraghavendra/PacBio/reads/short/day0_rep2_R1.fastq /gpfs/commons/home/sraghavendra/PacBio/reads/short/day0_rep2_R2.fastq > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02_short.sam
# #minimap2 -ax sr -t4 /gpfs/commons/home/sraghavendra/PacBio/reference/GRCh38_latest_genomic.fna /gpfs/commons/home/sraghavendra/PacBio/reads/short/transc_minimap/day0_rep2_R1.fastq.gz /gpfs/commons/home/sraghavendra/PacBio/reads/short/transc_minimap/day0_rep2_R2.fastq.gz > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02.sam
# samtools flagstat /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02_short.sam > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/report_02_short.txt
# samtools view -bS /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02_short.sam > /gpfs/commons/home/sraghavendra/PacBio/alignments/short/transc_minimap/aln_02_short.bam

minimap2 -ax splice -uf -C5 /gpfs/commons/home/sraghavendra/PacBio/reference/sota/transcriptome.fna /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/PacBio.simulated.fasta > /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_pb_2.sam
samtools view -bS /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_pb_2.sam > /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_pb_2.bam
samtools flagstat /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_pb_2.bam > /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_pb_2_report.txt

