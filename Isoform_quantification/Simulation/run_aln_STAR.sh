#!/bin/bash
#SBATCH --job-name=ill_aln_sim                      # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load star/2.5.2a
module load samtools/1.9

STAR --genomeDir /gpfs/commons/home/sraghavendra/Simulation/reference/ \
--runThreadN 6 \
--readFilesIn /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/Illumina.simulated_1.fq /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/Illumina.simulated_2.fq \
--outFileNamePrefix /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_ill_2_ \
--outSAMtype BAM SortedByCoordinate \

samtools flagstat /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_ill_2_Aligned.sortedByCoord.out.bam > /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_ill_2_report.txt

