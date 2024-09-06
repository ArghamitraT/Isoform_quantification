import os
import time
import random

trimester = time.strftime("_%Y_%m_%d__%H_%M_%S")


def create_job_dir(dir="", fold_name=""):
    if dir:
        job_path = os.path.join(dir, fold_name)
        os.mkdir(job_path)

    else:
        job_path = os.path.join(os.getcwd(), fold_name)
        if not os.path.exists(job_path):
            os.mkdir(job_path)

    return job_path


job_path = create_job_dir(dir="/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/files/NanoCount_output/",
                          fold_name="job")


### it calls the .py file
def create_prg_file(python_file_path, l0, l1, l2, l3, prg_file_path):
    header = f"#!/bin/bash\n" + \
             "module purge\n" + \
             "module load Anaconda3/2021.05\n" + \
             "source activate xmodality2\n" + \
             f"python {python_file_path} --stp_sz {l0} --l1 {l1} --GATlayer {l2} --dropout {l3} --weight_dir {weight_dir}  --epoch 450 --resume 1 --batch_size 2 --data_prlal 1"

    """ **CHANGE** THE PREVIOUS LINE """

    with open(prg_file_path, "w") as f:
        f.write(header)
    return prg_file_path


def copy_weights(desired_weight_path, to_be_saved_path):
    for file_name in os.listdir(desired_weight_path):
        dir = os.path.join(desired_weight_path, file_name)
        os.system(f"cp {dir} {to_be_saved_path}")


def gen_dir():
    dir = "terra_output_afterMLSB_ASA/" + path_weights_to_copy + "/weights"
    desired_weight_path = os.path.join(os.getcwd(), dir)
    to_be_saved_path = os.path.join(os.getcwd(), "weights")
    copy_weights(desired_weight_path, to_be_saved_path)


def create_slurm_file(prg_file_path, job_name, slurm_file_path):
    header = f"#!/bin/bash\n" + \
             "##ENVIRONMENT SETTINGS; REPLACE WITH CAUTION\n" + \
             "#SBATCH --export=NONE                #Do not propagate environment\n" + \
             "#SBATCH --get-user-env=L             #Replicate login environment\n" + \
             "##NECESSARY JOB SPECIFICATIONS\n" + \
             f"#SBATCH --job-name={job_name}      #Set the job name to \"JobExample1\"\n" + \
             "#SBATCH --time=65:00:00              #Set the wall clock limit to 1hr and 30min\n" + \
             "##SBATCH --time=45:00              #Set the wall clock limit to 1hr and 30min **CHANGE**\n" + \
             "#SBATCH --ntasks=48                   #Request 1 task\n" + \
             "#SBATCH --mem=180000M                  #Request 2560MB (2.5GB) per node **CHANGE**\n" + \
             f"#SBATCH --output={output_dir}/out_{job_name}.%j      #Send stdout/err to\n" + \
             "#SBATCH --gres=gpu:a100:2                    #Request 2 GPU per node can be 1 or 2 \n" + \
             "##OPTIONAL JOB SPECIFICATIONS\n" + \
             "#SBATCH --account=132755309533             #Set billing account to 123456\n" + \
             "##SBATCH --mail-type=ALL              #Send email on all job events\n" + \
             f"{prg_file_path}"

    with open(slurm_file_path, "w") as f:
        f.write(header)
    return slurm_file_path


def get_file_name(kind, l0, l1, l2, l3, ext=True):
    """ **CHANGE** THE FILE NAME """
    file_name = f"{kind}_stpsz_{l0}_l1_{l1}_dropout_{l3}_GATlayer_{l2}"
    if ext:
        file_name = f"{file_name}.sh"
    return file_name


""" **CHANGE** THE MAIN FOLDER NAME """
main_data_dir = create_job_dir(fold_name="terra_output_afterMLSB_ASA")

data_dir_0 = create_job_dir(dir=main_data_dir, fold_name="exprmnt" + trimester)
data_dir = create_job_dir(dir=data_dir_0, fold_name="files")
weight_dir = create_job_dir(dir=data_dir_0, fold_name="weights")
output_dir = create_job_dir(dir=data_dir, fold_name="output_files")

""" **CHANGE** THE DIRECTORY NAME FROM WHERE WEIGHTS NEEDS TO BE COPIED INTO ./WEIGHTS FOLDER(THE UNIVERSAL WEIGHT FOLDER)"""
path_weights_to_copy = "exprmnt_2022_09_02__23_39_10"
""" **CHANGE** IF WE NEED TO COPY THE WEIGHTS (LIKE FOR RESUME OF EVALUATION), THEN MAKE IT 1 OR MAKE IT 0"""
copy_needed = 0

if copy_needed:
    gen_dir()


def create_readme():
    name = os.path.join(data_dir, "readme")
    readme = open(name, "a")

    """ **CHANGE** WRITE THE COMMENT"""
    comment = f"Final PPI model, BERT frozen; training for 150 more epochs: total 450 epochs"

    readme.write(comment)
    readme.close()


def gen_combination():
    create_readme()

    """ **CHANGE** X0, X1, X2, NAME """
    for name in ['pretrain_ASA.py']:

        for x0 in [0.0001]:  # x0, l0 -> is step size
            for x1 in [0.0001]:  # x1, l1 -> is l1 or l2
                for x2 in [1, 2, 3, 4]:  # x2, l2 -> is GATlayer
                    for x3 in [0.1, 0.2]:  # x3, l3 -> is dropout
                        kind = name.split('_')[0]
                        python_file_path = os.path.join(os.getcwd(), name)
                        utility_file_path = os.path.join(os.getcwd(), "pretrain_utils_1.py")

                        hash_obj = random.getrandbits(25)
                        prg_file_path = os.path.join(job_path,
                                                     get_file_name(kind=f"prg_{kind}_{hash_obj}", l0=x0, l1=x1, l2=x2,
                                                                   l3=x3))
                        slurm_file_path = os.path.join(job_path,
                                                       get_file_name(kind=f"slurm_{kind}_{hash_obj}", l0=x0, l1=x1,
                                                                     l2=x2, l3=x3))

                        create_prg_file(python_file_path=python_file_path, prg_file_path=prg_file_path, l0=x0, l1=x1,
                                        l2=x2, l3=x3)
                        create_slurm_file(prg_file_path=prg_file_path,
                                          job_name=get_file_name(kind=kind, l0=x0, l1=x1, l2=x2, l3=x3, ext=False),
                                          slurm_file_path=slurm_file_path)

                        os.system(f"cp {python_file_path} {data_dir}")
                        os.system(f"cp {utility_file_path} {data_dir}")

                        # os.system(f"chmod u+x {prg_file_path}")
                        # os.system(f"chmod u+x {slurm_file_path}")
                        # os.system(f"sbatch {slurm_file_path}")


def main():
    gen_combination()


if __name__ == "__main__":
    main()

"""
MAY NEED LATER
#main_data_dir = os.path.join(pro_dir, "./terra_output_task2_WrmS/")
python_file_path = "/scratch/user/arghamitra.talukder/WORK/PPI/task2_iSqnc_oRr.py"
"""