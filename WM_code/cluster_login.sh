
# Testing JOLI in the Cluster
    # Wayne Monical
    # 4/2/2025

# This shell script outlines all of the steps I took to get JOLI running in the NYGC cluster
# Do not just run this script straight through! Run it line by line, and make sure each step works


# get to right directory
cd "/mnt/c/Users/Wayne Monical/OneDrive/Documents/JOLI"

# login to cluster 
# it will ask for password
    # ssh wmonical@10.4.30.43 # old logins
    # ssh wmonical@pe1-login # old logins
ssh wmonical@ne1-login ## 10.4.51.122

# load necessary modules
module avail
module load miniconda3

# get the code from git
git clone https://github.com/ArghamitraT/JOLI.git

# change to the folder
cd JOLI

# in order to use conda, YOU NEED TO INITIALIZE IT, but only the first time
conda init

# create a virtual environment to test JOLI. It has to be python version 3.10!
conda create --name joli_test python=3.10
conda activate joli_test

# install the required packages for JOLI
pip install -r Environments/JOLI_pip_requirements.txt

# exit the ssh
exit


## getting data to the cluster
cd 'mnt/c/Users/Wayne Monical/OneDrive/Documents/JOLI'

# copy the zipped data file from your local machine to the cluster
# it will ask for your password
scp example_data.zip wmonical@10.4.30.43:/gpfs/commons/home/wmonical/JOLI

## log back into ssh, then unzip
ssh wmonical@10.4.30.43
cd JOLI
unzip example_data.zip 

# running JOLI
conda activate joli_test
python main.py --data_folder 'Data/' \
               --output_path 'results/' \
               --sample1 'aln_E0' \
               --sample2 'aln_E2' \
               --experiment_num 4\
               --alpha_initial 10000\
               --max_em_rounds 30



# srun -p ne1 --mem 20G --cpus-per-task=4 --pty -u bash


