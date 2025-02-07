"""
This files is used to submit files in the slurm
"""
import os
import time
import random

trimester = time.strftime("_%Y_%m_%d__%H_%M_%S")
def create_job_dir(dir="", fold_name = ""):
    if dir:
        job_path = os.path.join(dir, fold_name)
        os.mkdir(job_path)

    else:
        job_path = os.path.join(os.getcwd(), fold_name)
        if not os.path.exists(job_path):
            os.mkdir(job_path)

    return job_path


#job_path = create_job_dir(fold_name="job")

### it calls the .py file
def create_prg_file(python_file_path, prg_file_path, output_file_path, input_file_names, alpha_initial, GD_lr, EM_round, load_filename, load, experiment_num):
   
    if old_prg_file:
        if experiment_num == 5:
            header = f"#!/bin/bash\n" + \
            "set -e\n" + \
            "cd $HOME\n" + \
            "source ~/.bashrc\n" + \
            "conda activate NanoCount_5\n" + \
            f"python {python_file_path} --data_folder {input_data_folder} --output_path {output_file_path} --sample1 {input_file_names[0][0]} {input_file_names[0][1]} --sample2 {input_file_names[1][0]} {input_file_names[1][1]} --alpha_initial {alpha_initial} --GD_lr {GD_lr} --max_em_rounds {EM_round} --load {load} --load_filename {load_filename} --experiment_num {experiment_num}"

        else:
            header = f"#!/bin/bash\n" + \
            "set -e\n" + \
            "cd $HOME\n" + \
            "source ~/.bashrc\n" + \
            "conda activate NanoCount_5\n" + \
            f"python {python_file_path} --data_folder {input_data_folder} --output_path {output_file_path} --sample1 {input_file_names[0]}  --sample2 {input_file_names[1]} --alpha_initial {alpha_initial} --GD_lr {GD_lr} --max_em_rounds {EM_round} --load {load} --load_filename {load_filename} --experiment_num {experiment_num} --dirichlet_builtin {dirichlet_builtin}"
        
        with open(prg_file_path, "w") as f:
            f.write(header)
    else:
        if experiment_num == 5:
            header = f"#!/bin/bash\n" + \
            "set -e\n" + \
            "cd $HOME\n" + \
            "source ~/.bashrc\n" + \
            "conda activate NanoCount_5\n" + \
            f"python {python_file_path} --data_folder {input_data_folder} --output_path {output_file_path} --sample1 {input_file_names[0]} {input_file_names[1]} --sample2 {input_file_names[2]} {input_file_names[3]} --alpha_initial {alpha_initial} --GD_lr {GD_lr} --max_em_rounds {EM_round} --load {load} --load_filename {load_filename} --experiment_num {experiment_num} --dirichlet_builtin {dirichlet_builtin} --EM_type {EM_type} --dirichlet_process {dirichlet_process} --process_bam_required {process_bam_required}"

        else:
            header = f"#!/bin/bash\n" + \
            "set -e\n" + \
            "cd $HOME\n" + \
            "source ~/.bashrc\n" + \
            "conda activate NanoCount_5\n" + \
            f"python {python_file_path} --data_folder {input_data_folder} --output_path {output_file_path} --sample1 {input_file_names[0]}  --sample2 {input_file_names[1]} --alpha_initial {alpha_initial} --GD_lr {GD_lr} --max_em_rounds {EM_round} --load {load} --load_filename {load_filename} --experiment_num {experiment_num} --dirichlet_builtin {dirichlet_builtin} --EM_type {EM_type} --dirichlet_process {dirichlet_process} --process_bam_required {process_bam_required}"
        
        with open(prg_file_path, "w") as f:
            f.write(header)
    return prg_file_path
    
    
def copy_weights(desired_weight_path, to_be_saved_path ):
    for file_name in os.listdir(desired_weight_path):
        dir = os.path.join(desired_weight_path, file_name)
        os.system(f"cp {dir} {to_be_saved_path}")


def create_slurm_file(prg_file_path, job_name, slurm_file_path):

    show_name = '_'.join(job_name.split('_')[1:])
    show_name = f"{slurm_file_name}_{show_name}"


    header = f"#!/bin/bash\n" + \
    "##ENVIRONMENT SETTINGS; REPLACE WITH CAUTION\n" + \
    "##NECESSARY JOB SPECIFICATIONS\n" + \
    f"#SBATCH --job-name={show_name}      #Set the job name to \"JobExample1\"\n" + \
    f"#SBATCH --time={hour}:45:00              #Set the wall clock limit \n" + \
    f"#SBATCH --mem={memory}G              \n" + \
    f"#SBATCH --cpus-per-task={nthred}                   \n" + \
    "#SBATCH --mail-type=END,FAIL    \n" + \
    f"#SBATCH --output={output_dir}/out_{job_name}.%j      #Send stdout/err to\n" + \
    "#SBATCH --mail-user=atalukder@nygenome.org                    \n" + \
    f"{prg_file_path}"

    with open (slurm_file_path, "w") as f:
        f.write(header)
    return slurm_file_path


def get_file_name(kind, l0=0, l1=0, l2=0, l3=0, ext=True):

    file_name = f"{kind}_{trimester}"
    if ext:
        file_name = f"{file_name}.sh"
    return file_name



main_data_dir = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results"
job_path = "/gpfs/commons/home/atalukder/RNA_Splicing/files/cluster_job_submission_files"
code_dir = "/gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code"
where_to_save = "/gpfs/commons/home/atalukder/RNA_Splicing/files/results/saved_weights_toResume"

data_dir_0   = create_job_dir(dir= main_data_dir, fold_name= "exprmnt"+trimester)
data_dir   = create_job_dir(dir= data_dir_0, fold_name= "files")
weight_dir = create_job_dir(dir= data_dir_0, fold_name="weights")
output_dir = create_job_dir(dir= data_dir, fold_name="output_files")

""" Parameters: **CHANGE (AT)** """
# **** the first element should be LR, the second should be SR ****
samples_file_names = [['ds_100_num1_aln_01_long', 'ds_100_num1_aln_21_short'],
                      ['ds_100_num1_aln_21_long', 'ds_100_num1_aln_01_short']]
# samples_file_names = [['ds_100_num1_aln_01_long', 'ds_100_num1_aln_01_short'],
#                       ['ds_100_num1_aln_21_long', 'ds_100_num1_aln_21_short']]
# samples_file_names = [['ds_100_num1_aln_01_long', 'NA'],
#                         ['ds_100_num1_aln_01_short', 'NA'],
#                         ['ds_100_num1_aln_21_long', 'NA'], 
#                        ['ds_100_num1_aln_21_short', 'NA']]
# samples_file_names = [['ds_10_num1_aln_52_long', 'ds_100_num1_aln_02_short']]
# samples_file_names = [['ds_10_num1_aln_52_long', 'NA'],
#                        ['ds_10_num1_aln_51_long', 'NA']]
other_script_names = ['EM_VIorMAP_GD_vector.py', 'DirichletOptimizer_vector.py', 'generate_bash.py', 'process_bam_files.py']
name_arr = ['main_EM_VIorMAP_GD_vector.py']
from_where_to_copy = "exprmnt_2024_08_10__02_05_36"
last_EM_round = 25
copy_needed = 0 #(AT)
alpha_val_arr = [10000]
GDlr_val_arr = [0.01]
EM_round_arr = [30] 
experiment_num = 4      #"Different experiment setup, 1: for 1 sample, 2 for merged, 4 for multisample, 5 for merged multisample"
dirichlet_builtin = 0
simulation = 1
old_prg_file = 0
hour=7
memory=100 # GB
nthred = 8 # number of CPU
EM_type = 'MAP'  # "Inference canbe through VI or MAP"
dirichlet_process = 'theta' # "Dirichlet optimization canbe 'expectation_log_theta' or 'theta'"
process_bam_required = 0
readme_comment = f"conference experiment, exp4, simulation, MAP, reduced delta theta to 0.001, newly aligned LR .bam files with new minimap version 2.28"
if simulation:
    input_data_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/pklfiles/new_minimap_file/'
    output_file_name = "output_Simulation_VIGD_token_"
else:
    # input_data_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/subfolder/'
    input_data_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln_pklfiles/'
    output_file_name = "output_PacIllu_VIGD_token_"
""" Parameters: **CHANGE (AT)** """ 


slurm_file_name = 'exp4SimMAP'

def create_readme():
    name = os.path.join(data_dir, "readme")
    readme = open(name, "a")
    comment = readme_comment
    readme.write(comment)
    readme.close()


def create_load_file_path(file_names_list, GD_lr, alpha_initial):
    default_load_filepath = os.path.join(main_data_dir, from_where_to_copy, 'weights')
    fileName = "allWeights"     # Initialize an empty string to store the result

    # Iterate through the list of file names
    for index, file_path in enumerate(file_names_list, start=1):
        file_name = file_path.split('/')[-1]
        file_identifier = ''.join(file_name.split('_')).split('.')[0]
        fileName += f"_file{index}_{file_identifier}"
    fileName = f"{fileName}_GDlr_{GD_lr}_AlphaInitial_{alpha_initial}.0_EMround_{last_EM_round}"
    file_name = fileName.strip()
    
    # Search for files that match this pattern in the specified directory
    file_list = []
    for file in os.listdir(default_load_filepath):
        if file.startswith(file_name) and file.endswith('.pkl'):
            file_list.append(os.path.join(default_load_filepath, file))
 
    if not file_list:
        print(f"No file found matching the format {fileName}'")
        return
    elif len(file_list)>1:
        print(f"More than 1 file found matching the format {fileName}'")
        return
    else:
        return file_list[0]

def gen_combination():
    create_readme()
    
    for name in name_arr:
        for set in samples_file_names:                
            for alpha_val in alpha_val_arr:   
                for GDlr_val in GDlr_val_arr:    
                    for EM_round in EM_round_arr:  

                        kind = name.split('_')[0]
                        python_file_path = os.path.join(code_dir, name)

                        hash_obj = random.getrandbits(25)
                        
                        prg_file_path = os.path.join(job_path, get_file_name(kind= f"prg_{kind}_{hash_obj}"))
                        slurm_file_path = os.path.join(job_path, get_file_name(kind= f"slurm_{kind}_{hash_obj}"))
                        output_file_path = os.path.join(output_dir, f"{output_file_name}{hash_obj}")
                        if copy_needed:
                            load_filename = create_load_file_path(file_names_list=set, GD_lr=GDlr_val, alpha_initial=alpha_val)
                        else:
                            load_filename = 'generic'
                        
                        # create_prg_file(python_file_path=python_file_path, prg_file_path=prg_file_path, output_file_path=output_file_path, input_file_names=set, alpha_initial=alpha_val)
                        create_prg_file(python_file_path=os.path.join(data_dir, name), 
                                        prg_file_path=prg_file_path, 
                                        output_file_path=output_file_path, 
                                        input_file_names=set, 
                                        alpha_initial=alpha_val, 
                                        GD_lr= GDlr_val,
                                        EM_round = EM_round, 
                                        load_filename=load_filename, 
                                        load=copy_needed,
                                        experiment_num=experiment_num)
                        
                        
                        create_slurm_file(prg_file_path=prg_file_path, 
                                        job_name=get_file_name(kind=f"{kind}_{hash_obj}", ext=False), 
                                        slurm_file_path=slurm_file_path)

                        os.system(f"cp {python_file_path} {data_dir}")
                        for script in other_script_names:
                            script_dir = os.path.join(code_dir, script)
                            os.system(f"cp {script_dir} {data_dir}")

                        # #os.system(f"cp {utility_file_path} {data_dir}")
                        ## (AT)
                        os.system(f"chmod u+x {prg_file_path}")
                        os.system(f"chmod u+x {slurm_file_path}")
                        os.system(f"sbatch {slurm_file_path}")
                    

def main():
    gen_combination()

if __name__ == "__main__":
    main()
