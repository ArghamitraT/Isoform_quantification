#!/bin/bash
#SBATCH --job-name=ds_distr                    # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=200G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load python/3.9.7
source /gpfs/commons/home/sraghavendra/nanocount_env/bin/activate

#python3 get_read_distribution.py Encode/Alignments/short_aln.bam Encode/enc_short.png
#python3 AlignmentStatistics_shree.py Encode/Alignments/short_aln.bam > Encode/report_short.txt
#python3 get_read_distribution.py PacBio/alignments/short/transcriptome/aln_01_short_Aligned.sortedByCoord.out.bam PacBio/ill_short_trans_STAR_01.png
#python3 AlignmentStatistics_shree.py PacBio/alignments/short/transc_minimap/aln_01_short.bam > PacBio/report_short_trans_min_01.txt
#python3 get_read_distribution.py PacBio/alignments/short/transcriptome/aln_02_short_Aligned.sortedByCoord.out.bam PacBio/ill_short_trans_STAR_02.png
#python3 AlignmentStatistics_shree.py PacBio/alignments/short/transcriptome/aln_02_short_Aligned.sortedByCoord.out.bam > PacBio/report_short_trans_STAR_02.txt
#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/aln_02_long.bam PacBio/pb_long_trans_02.png
#python3 AlignmentStatistics_shree.py PacBio/alignments/long/transcriptome/aln_01_long.bam > PacBio/report_long_trans_01.txt
#python3 AlignmentStatistics_shree.py PacBio/alignments/long/transcriptome/aln_02_long.bam > PacBio/report_long_trans_02.txt
#python3 isoform_stats.py PacBio/alignments/short/transc_minimap/aln_02_short.bam PacBio/alignments/long/transcriptome/aln_01_long.bam > PacBio/reps/common_report_trans_min_01.txt
python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_01_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/aln_02_long.bam > PacBio/reps_cross/long_02_short_01.txt$

#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/ds_10_aln_01_long.bam PacBio/imgs/ds_10_01.png

#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/ds_10_aln_02_long.bam PacBio/imgs/ds_10_02.png

#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/ds_8_aln_01_long.bam PacBio/imgs/ds_8_01.png

#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/ds_8_aln_02_long.bam PacBio/imgs/ds_8_02.png

#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/ds_5_aln_01_long.bam PacBio/imgs/ds_5_01.png

#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/ds_5_aln_02_long.bam PacBio/imgs/ds_5_02.png

#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/ds_2_aln_01_long.bam PacBio/imgs/ds_2_01.png

#python3 get_read_distribution.py PacBio/alignments/long/transcriptome/ds_2_aln_02_long.bam PacBio/imgs/ds_2_02.png


#python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_02_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/ds_10_aln_01_long.bam > PacBio/reps_cross/common_ds_10_01.txt

#python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_02_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/ds_8_aln_01_long.bam > PacBio/reps_cross/common_ds_8_01.txt

#python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_02_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/ds_5_aln_01_long.bam > PacBio/reps_cross/common_ds_5_01.txt

#python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_02_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/ds_2_aln_01_long.bam > PacBio/reps_cross/common_ds_2_01.txt

#python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_01_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/ds_10_aln_02_long.bam > PacBio/reps_cross/common_ds_10_02.txt

#python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_01_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/ds_8_aln_02_long.bam > PacBio/reps_cross/common_ds_8_02.txt

#python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_01_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/ds_5_aln_02_long.bam > PacBio/reps_cross/common_ds_5_02.txt

#python3 isoform_stats.py PacBio/alignments/short/transcriptome/aln_01_short_Aligned.sortedByCoord.out.bam PacBio/alignments/long/transcriptome/ds_2_aln_02_long.bam > PacBio/reps_cross/common_ds_2_02.txt

