#!/bin/bash
##ENVIRONMENT SETTINGS; REPLACE WITH CAUTION
##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=Simulation      #Set the job name too "JobExample1"
#SBATCH --time=8:45:00              #Set the wall clock limit to 1hr and 30min,takes 100min/EM iteration **CHANGE (AT)**
#SBATCH --mem=256G              
#SBATCH --cpus-per-task=8                   
#SBATCH --mail-type=END,FAIL    
#SBATCH --output=/gpfs/commons/home/atalukder/RNA_Splicing/files/results/simulations/files/output_files/simulation.%j      #Send stdout/err to
#SBATCH --mail-user=atalukder@nygenome.org 

set -e
cd $HOME
source ~/.bashrc
conda activate NanoCount_5
python /gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code/simulation.py