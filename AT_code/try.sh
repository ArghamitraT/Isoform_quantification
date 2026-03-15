#!/bin/bash
##ENVIRONMENT SETTINGS; REPLACE WITH CAUTION
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=Simulation      #Set the job name too "JobExample1"
#SBATCH --time=15:45:00              #Set the wall clock limit to 1hr and 30min,takes 100min/EM iteration **CHANGE (AT)**
#SBATCH --mem=30G              
#SBATCH --cpus-per-task=8                   
#SBATCH --mail-type=END,FAIL    
#SBATCH --output=/gpfs/commons/home/atalukder/RNA_Splicing/files/results/simulations/files/output_files/simulation.%j      #Send stdout/err to
#SBATCH --mail-user=atalukder@nygenome.org 



# cp -r /gpfs/commons/home/spark/knowles_lab/Argha /gpfs/commons/home/atalukder/RNA_Splicing/data/knowles_lab
# cp -r /gpfs/commons/home/sraghavendra/PacBio/ /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_fastq

cp -r /gpfs/commons/home/sraghavendra/Simulation /gpfs/commons/groups/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/simulations_copiedFromShree



# set -e
# cd $HOME
# source ~/.bashrc
# conda activate NanoCount_5
# # python /gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code/try.py
# # cp -r /gpfs/commons/home/sraghavendra/SOTA /gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff
# rsync -ah --progress /gpfs/commons/home/sraghavendra/PacBio /gpfs/commons/home/atalukder/RNA_Splicing/data/Shree_stuff
