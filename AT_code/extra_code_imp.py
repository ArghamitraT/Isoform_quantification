import os
import shutil
import fnmatch


x_25 = -357186225.84963673
x_26 = -357161032.7245745

print(x_25-x_26)

x_49 = -356868260.53346413
x_50 = -356863429.0546106

print(x_49-x_50)

def delet_weight_folder(destination_dir):
    for folder in os.listdir(destination_dir):
        folder_path = os.path.join(destination_dir, folder)

        for subfolder in os.listdir(folder_path):
            subfolder_path = os.path.join(folder_path, subfolder)

            if os.path.isdir(subfolder_path):
                if subfolder == 'weights':
                    shutil.rmtree(subfolder_path)


def delet_output_folder(destination_dir):
    for folder in os.listdir(destination_dir):
        folder_path = os.path.join(destination_dir, folder)

        for subfolder in os.listdir(folder_path):
            subfolder_path = os.path.join(folder_path, subfolder)
            
            for minifolder in os.listdir(subfolder_path):
                minifolder_path = os.path.join(subfolder_path, minifolder)
                if os.path.isdir(minifolder_path):
                    if minifolder == 'output_files':
                        shutil.rmtree(minifolder_path)



# this functions copies the result folder without heavy self weights
def copy_except_weights_and_files(source_dir, destination_dir):
    # Step 1: Go inside the source directory
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Step 2: Read all folders in the directory
    for folder_name in os.listdir(source_dir):
        folder_path = os.path.join(source_dir, folder_name)

        # Step 3: If the folder is 'weights', handle it separately
        for sub_folder_name in os.listdir(folder_path):
            sub_folder_path = os.path.join(folder_path, sub_folder_name)
            if os.path.isdir(sub_folder_path):
                if sub_folder_name == 'weights':
                    # Handle copying of files in the 'weights' folder, skipping 'allWeights_*.pkl'
                    for root, dirs, files in os.walk(sub_folder_path):
                        for file in files:
                            if not fnmatch.fnmatch(file, 'allWeights_*.pkl'):
                                # Copy the file if it does not match 'allWeights_*.pkl'
                                source_file = os.path.join(root, file)
                                dest_dir = root.replace(source_dir, destination_dir, 1)
                                if not os.path.exists(dest_dir):
                                    os.makedirs(dest_dir)
                                shutil.copy2(source_file, dest_dir)
                else:
                    # Copy everything else in the folder to the destination directory
                    dest_folder_path = os.path.join(destination_dir, folder_name, sub_folder_name)
                    shutil.copytree(sub_folder_path, dest_folder_path, dirs_exist_ok=True)



# this functions copies the result folder without heavy self weights
def copy_except_weights_and_files_expFolder(source_dir, destination_dir):
    # Step 1: Go inside the source directory
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Step 2: Read all folders in the directory
    #for folder_name in os.listdir(source_dir):
        #folder_path = os.path.join(source_dir, folder_name)
    
    folder_name = source_dir.split('/')[-1]
    folder_path = source_dir
    # Step 3: If the folder is 'weights', handle it separately
    for sub_folder_name in os.listdir(folder_path):
        sub_folder_path = os.path.join(folder_path, sub_folder_name)
        #sub_folder_path = folder_path
        print(sub_folder_path)
        if os.path.isdir(sub_folder_path):
            if sub_folder_name == 'weights':
                # Handle copying of files in the 'weights' folder, skipping 'allWeights_*.pkl'
                for root, dirs, files in os.walk(sub_folder_path):
                    for file in files:
                        if not fnmatch.fnmatch(file, 'allWeights_*.pkl'):
                            # Copy the file if it does not match 'allWeights_*.pkl'
                            source_file = os.path.join(root, file)
                            dest_dir = root.replace(source_dir, destination_dir, 1)
                            if not os.path.exists(dest_dir):
                                os.makedirs(dest_dir)
                            shutil.copy2(source_file, dest_dir)
            else:
                # Copy everything else in the folder to the destination directory
                dest_folder_path = os.path.join(destination_dir, folder_name, sub_folder_name)
                shutil.copytree(sub_folder_path, dest_folder_path, dirs_exist_ok=True)




def copy_scripts(script_name, destination):
    shutil.copy(script_name, destination)


# destination_folder =  'exprmnt_2024_09_29__22_02_56'
# script_names_list = ['generate_result_stat.py', 'plt_experiment_stats.py']
# script_main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code'
# destination_main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results' 

# for script in script_names_list:
#     destination_final_dir = os.path.join(destination_main_dir, destination_folder, 'files')
#     script_final_dir = os.path.join(script_main_dir, script)
#     copy_scripts(script_final_dir, destination_final_dir)


import os
import argparse
parser = argparse.ArgumentParser(description="Process BAM files and output results.")
parser.add_argument("--input_folder", type=str,
                    help="Path for the output file. File name should be in this format 'outputTRIAL_PacIllu_VIGD_token_00000'. Default is files/results/exprmntSingleRun_2024_00_00__00_00_00/files/output_files/outputTRIAL_PacIllu_VIGD_token_00000")

args = parser.parse_args()
# Assign the user-provided or default path to the output_file variable
input_folder = args.input_folder

# Usage
experiment_folder = input_folder
source_directory = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/'+experiment_folder
destination_directory = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/Recomb_25_seq_experiments/with_pnm_Noweight/'+experiment_folder
# #delet_weight_folder('/gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code/extra_codes/new_simulation_results2')
#copy_except_weights_and_files(source_directory, destination_directory)
copy_except_weights_and_files_expFolder(source_directory, destination_directory)
#delet_output_folder('/gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code/extra_codes/new_simulation_results2')


