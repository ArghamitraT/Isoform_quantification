#!/bin/bash
#SBATCH --job-name=STAR_gen_inds                    # Job name
#SBATCH --partition=pe2                             # Partition Name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=150G                                   # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load star/2.5.2a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /gpfs/commons/home/sraghavendra/Simulation/reference/ \
--limitGenomeGenerateRAM 150000000000 \
--genomeFastaFiles /gpfs/commons/home/sraghavendra/PacBio/reference/sota/transcriptome.fna
#--sjdbGTFfile /gpfs/commons/home/sraghavendra/PacBio/reference/transc/genomic.gtf \
#--sjdbOverhang 99
