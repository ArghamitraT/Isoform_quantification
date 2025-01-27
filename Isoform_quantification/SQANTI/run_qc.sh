#!/bin/bash
#SBATCH --job-name=sqanti_qc                        # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=200G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

#cd /gpfs/commons/home/sraghavendra/SOTA/MPAQT/MPAQT
~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate SQANTI3.env

#python SQANTI3-5.2.1/sqanti3_qc.py /gpfs/commons/home/sraghavendra/bam2gff/aln_01_long_ds_2.gff /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna

python SQANTI3/sqanti3_qc.py --fasta /gpfs/commons/home/sraghavendra/SQANTI/ds_10_aln_01_long.fastq /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna --force_id_ignore --skipORF
