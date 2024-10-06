from EM_VI_GD import Expec_Max
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
    
    # (AT)
    #file_list = ['/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_2024_08_10__02_05_36/weights/allWeights_file1_ds2num1aln02long_file2_ds100num1aln02short_GDlr_0.01_AlphaInitial_10000.0_EMround_25_token_29675364_2024_8_10_02_24_57.pkl']
 
    if not file_list:
        print(f"No file found matching the format {fileName}'")
        return
    elif len(file_list)>1:
        print(f"More than 1 file found matching the format {fileName}'")
        return
    else:
        return file_list[0]


##main_folder = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/data/'

main_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/'
output_file_default = os.path.join(os.getcwd(), '../../files/results/exprmntSingleRun_2024_00_00__00_00_00/files/output_files/outputTRIAL_PacIllu_VIGD_token_00000')

# Set up argparse to handle command-line arguments (AT) TAKE A CLOSE LOOK AT THE DEFAULT VALUES
parser = argparse.ArgumentParser(description="Process BAM files and output results.")
parser.add_argument("--output_path", type=str, default=output_file_default,
                    help="Path for the output file. Default is '../../files/results/single_run/files/output_files/output_generic_VIGD_00000'.")
parser.add_argument("--sample1", type=str, default='ds_5_aln_02_long.bam', help="Sample1 (LR) file name.")
parser.add_argument("--sample2", type=str, default='aln_01_short.bam', help="Sample2 (SR) file name.")
parser.add_argument("--GD_lr", type=float, default=0.01, help="Learning rate for dirichlet gradient descent.")
parser.add_argument("--alpha_initial", type=float, default=10000, help="The fixed sum value of alpha")
parser.add_argument("--max_em_rounds", type=int, default=2, help="The maximum EM iterations")
parser.add_argument("--load", type=int, default=0, help="0 or 1 If load is 1, we will load the old weights and resume training") ## (AT)
parser.add_argument("--load_filename", type=str, default='generic', help="if load==1 (the model needs to resumed training), need a path to load the weights")


# Parse the arguments
args = parser.parse_args()
# Assign the user-provided or default path to the output_file variable
output_file = args.output_path
sample1 = main_folder + args.sample1
sample2 = main_folder + args.sample2
GD_lr = args.GD_lr
alpha_initial = args.alpha_initial
max_em_rounds = args.max_em_rounds
load = args.load
load_filename = args.load_filename


## (AT) SIRV
main_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/SIRV/'
# sample1 = main_folder + 'aln_E0_short.bam'
# sample2 = main_folder + 'aln_E2_short.bam'
sample1 = main_folder + 'aln_E0.bam'
sample2 = main_folder + 'aln_E2.bam'
# output_file = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmntSingleRun_2024_00_00__00_00_00/output_SIRV_VIGD_token_000000'
# sample1 = main_folder + 'ds_2_num1_aln_02_long.bam'
# sample2 = main_folder + 'ds_100_num1_aln_02_short.bam'
# output_file = os.path.join(os.getcwd(), '../../files/results/exprmnt_2024_08_28__01_31_50/files/outputTRIAL_PacIllu_VIGD_token_000000')
# main_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/'
# sample1 = main_folder + 'ds_100_num1_aln_01_short.bam'
# sample2 = main_folder + 'ds_5_num1_aln_51_long.bam'

# Print all the parameters
last_EM_round = 25
file_names_list = [sample1, sample2] ## (AT)
print("Output_file_path", output_file)
print("GD_lr", GD_lr)
print("alpha_initial", alpha_initial)
print("max_em_rounds", max_em_rounds)
print("load", load)

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
           load_filename=load_filename)

# takes 100min/EM iteration

print("#########END###########")#


"""
/gpfs/commons/home/atalukder/miniconda3/envs/NanoCount_5/bin/python /gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code/main_EM_VI_GD.py > /gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmntSingleRun_2024_00_00__00_00_00/files/out_main_000000__2024_07_13__22_52_00.txt 2>&1

"""
