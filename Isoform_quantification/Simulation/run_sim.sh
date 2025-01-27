module load samtools
module load minimap2
module load rsem


mkdir -p sim_result_2

tar -xzf data/human_chr22.tar.gz -C sim_result_2/

python prepare_reference_data.py \
  --reference_annotation /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf \
  --reference_transcripts /gpfs/commons/home/sraghavendra/PacBio/reference/sota/transcriptome.fna \
  --reference_genome /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna \
  --n_random_isoforms 50 \
  --output sim_result_2/ref_data/human

# Decided to simulate reads based on reference transcripts only
#--sqanti_prefix sim_result_2/human_chr22/sqanti \ (prep ref data)
#--mandatory sim_result_2/ref_data/human.novel_isoforms.tsv \ (quantify)


#Will simulate using our real PacBio data

python quantify.py \
  --fastq /gpfs/commons/home/sraghavendra/PacBio/reads/long/flnc_51.fastq \
  -t 16 \
  --reference_transcripts sim_result_2/ref_data/human.transcripts.fasta \
  --output sim_result_2/ref_data/human.counts.tsv

# python simulate.py \
#   --reference_prefix sim_result_2/ref_data/human \
#   --counts sim_result_2/ref_data/human.counts.tsv \
#   -t 16 -s 22 \
#   --output sim_result_2/human_simulated_10_reads/
