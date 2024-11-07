from EM_VI_GD_vector import Expec_Max
import os
import argparse


def create_load_file_path():
    default_load_filepath = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_2024_08_10__02_05_36/weights'
    # Initialize an empty string to store the result
    fileName = "allWeights"

    # Iterate through the list of file names
    for index, file_path in enumerate(file_names_list, start=1):
        # Split the file path by '/' and take the last part (the file name)
        file_name = file_path.split('/')[-1]
        # Extract a specific part of the file name if necessary (e.g., removing extension)
        file_identifier = ''.join(file_name.split('_')).split('.')[0]
        # Construct the string
        fileName += f"_file{index}_{file_identifier}"
    fileName = f"{fileName}_GDlr_{GD_lr}_AlphaInitial_{alpha_initial}.0_EMround_{last_EM_round}"
    file_name = fileName.strip()
    
    # Search for files that match this pattern in the specified directory
    file_list = []
    for file in os.listdir(default_load_filepath):
        if file.startswith(file_name) and file.endswith('.pkl'):
            file_list.append(os.path.join(default_load_filepath, file))
    
    #file_list = ['/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_2024_08_10__02_05_36/weights/allWeights_file1_ds2num1aln02long_file2_ds100num1aln02short_GDlr_0.01_AlphaInitial_10000.0_EMround_25_token_29675364_2024_8_10_02_24_57.pkl']
 
    if not file_list:
        print(f"No file found matching the format {fileName}'")
        return
    elif len(file_list)>1:
        print(f"More than 1 file found matching the format {fileName}'")
        return
    else:
        return file_list[0]


main_folder_default = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln_pklfiles/'
output_file_default = os.path.join(os.getcwd(), 'results/exprmntSingleRun_2024_00_00__00_00_00/files/output_files/outputTRIAL_PacIllu_VIGD_token_00000')

# Set up argparse to handle command-line arguments (AT) TAKE A CLOSE LOOK AT THE DEFAULT VALUES
parser = argparse.ArgumentParser(description="Process BAM files and output results.")
parser.add_argument("--data_folder", type=str, default=main_folder_default,
                    help="Path for the output file. File name should be in this format 'outputTRIAL_PacIllu_VIGD_token_00000'. Default is files/results/exprmntSingleRun_2024_00_00__00_00_00/files/output_files/outputTRIAL_PacIllu_VIGD_token_00000")

parser.add_argument("--output_path", type=str, default=output_file_default,
                    help="Path for the output file. File name should be in this format 'outputTRIAL_PacIllu_VIGD_token_00000'. Default is files/results/exprmntSingleRun_2024_00_00__00_00_00/files/output_files/outputTRIAL_PacIllu_VIGD_token_00000")
parser.add_argument("--sample1", type=str, default='aln_E0', help="Sample1 (LR) file name.")
parser.add_argument("--sample2", type=str, default='aln_E2', help="Sample2 (SR) file name.")
parser.add_argument("--GD_lr", type=float, default=0.01, help="Learning rate for dirichlet gradient descent.")
parser.add_argument("--alpha_initial", type=float, default=10000, help="The fixed sum value of alpha")
parser.add_argument("--max_em_rounds", type=int, default=2, help="The maximum EM iterations")
parser.add_argument("--load", type=int, default=0, help="0 or 1 If load is 1, we will load the old weights and resume training") ## (AT)
parser.add_argument("--load_filename", type=str, default='generic', help="if load==1 (the model needs to resumed training), need a path to load the weights")
parser.add_argument("--experiment_num", type=int, default=4, help="Different experiment setup, 1: for 1 sample, 2 for merged, 4 for multisample, 5 for merged multisample")
parser.add_argument("--dirichlet_builtin", type=int, default=0, help="if 1 then uses torch builtin dirichlet function, if 0 then implements the manual calculation")

# Parse the arguments
args = parser.parse_args()
# Assign the user-provided or default path to the output_file variable
input_folder = args.data_folder
output_file = args.output_path
experiment_num = args.experiment_num

if experiment_num == 5:
    sample1 = args.sample1
    sample2 = args.sample2
else:
    sample1 = input_folder + args.sample1
    sample2 = input_folder + args.sample2

GD_lr = args.GD_lr
alpha_initial = args.alpha_initial
max_em_rounds = args.max_em_rounds
load = args.load
load_filename = args.load_filename
dirichlet_builtin = args.dirichlet_builtin


## (AT) SIRV
# sample1 = input_folder  + 'aln_E0_short'
# sample2 = input_folder + 'aln_E2_short'

# Print all the parameters
last_EM_round = 25
if experiment_num == 1:
    file_names_list = [sample1]
else:
    file_names_list = [sample1, sample2]

print("Output_file_path", output_file)
print("GD_lr", GD_lr)
print("alpha_initial", alpha_initial)
print("max_em_rounds", max_em_rounds)
print("load", load)
print("dirichlet_builtin", dirichlet_builtin)
print("experiment_num", experiment_num)
if experiment_num == 1:
    print("Single sample, no gradient decent")
elif experiment_num == 2:
    print("2 samples merged as one, no gradient decent")
elif experiment_num == 4:
    print(print("Multi-sample with gradient decent"))
elif experiment_num == 5:
    print(print("Merged multi-sample with gradient decent"))


if load:
    print("load_filename", load_filename)
    if load_filename == "generic":
        load_filename= create_load_file_path()

Expec_Max (file_names=file_names_list, 
           count_file=output_file, 
           GD_lr=GD_lr, 
           alpha_initial=alpha_initial, 
           max_em_rounds=max_em_rounds,
           load=load,
           load_filename=load_filename,
           experiment_num = experiment_num,
           dirichlet_builtin = dirichlet_builtin)

print("#########END###########")#


"""
/gpfs/commons/home/atalukder/miniconda3/envs/NanoCount_5/bin/python /gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code/main_EM_VI_GD.py > /gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmntSingleRun_2024_00_00__00_00_00/files/out_main_000000__2024_07_13__22_52_00.txt 2>&1

"""
"""
git DirichletOptimizer.py EM_VI_GD.py find_theta_corr.py main_EM_VI_GD.py main_EM_VI_GD_simulation.py plt_experiment_stats.py simulation.py try2.py corrWGroundTruth.py find_theta_alphaInitial_relation.py vector_main.py vector_version.py 
"""