#!/bin/bash
#SBATCH --job-name=nc                               # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=100G                                   # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load python/3.9.7
source /gpfs/commons/home/sraghavendra/nanocount_env/bin/activate

# NanoCount -i /gpfs/commons/home/sraghavendra/PacBio/alignments/long/transcriptome/long_aln_51.bam -b /gpfs/commons/home/sraghavendra/SOTA/NanoCount/aligned_reads_selected_51.bam --extra_tx_info -o /gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_51.tsv


# NanoCount -i /gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/subfolder/ds_10_num1_aln_52_long.bam -b /gpfs/commons/home/sraghavendra/SOTA/NanoCount/aligned_reads_selected_52.bam --extra_tx_info -o /gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_ds_52.tsv

# NanoCount -i /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_pb_1.bam -b /gpfs/commons/home/sraghavendra/SOTA/NanoCount/aligned_reads_selected_sim_1.bam --extra_tx_info -o /gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_ds_sim_1.tsv

NanoCount -i /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_pb_2.bam -b /gpfs/commons/home/sraghavendra/SOTA/NanoCount/aligned_reads_selected_sim_2.bam --extra_tx_info -o /gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_ds_sim_2.tsv
