"""
Script Description:
====================
This script processes experiment data files, compares them with ground truth values, 
and calculates Spearman and Pearson correlation coefficients for results using specific parameters.

Inputs:
-------
1. `main_result_dir` (str): Main directory containing the experimental results.
2. `experiment_file` (str): Subdirectory under the main directory where the experiment results are stored.
3. `simulation` (int): Flag indicating whether the data is simulated (1) or real (0).
4. `experiment` (int): Experiment type identifier (e.g., 1, 2, 4, 5).
   - 1: Single sample
   - 2: Merged samples
   - 4: Multi-sample
   - 5: Merged multi-sample
5. `groundTruth_main_dir` (str): Directory containing ground truth data files.
6. `groundTruth_file` (dict): Mapping of sample identifiers to ground truth file names.

Outputs:
--------
1. Logs of extracted information and the computed correlation coefficients for each processed file.
2. Correlation calculations are handled by the `result_process.spearman_pearson_corr_generic` function.

Script Flow:
------------
1. **File Processing**:
   - Files are iterated over in `file_arr` (list of file names from the output directory).
   - For each file, the script extracts relevant metadata using `ExperimentFileProcessor`.
   - The metadata includes information such as sample type, length type (long/short), and replica details.

2. **Ground Truth Association**:
   - Based on the extracted metadata, the corresponding ground truth file is identified using `groundTruth_file`.

3. **Correlation Calculation**:
   - The script compares the predicted result (from the processed file) and the ground truth file.
   - Calculates Spearman and Pearson correlation coefficients using `result_process.spearman_pearson_corr_generic`.

4. **Dynamic File Handling**:
   - Handles both simulation and real data dynamically by adjusting parameters and file paths.

Notes:
------
- Ensure that the ground truth file paths are correctly specified in `groundTruth_file`.
- The directory structure for results and ground truth files must be consistent with the script's assumptions.

"""

from util import ExperimentFileProcessor
import os
import generate_result_stat as result_process

def process_experiment_data(directory, experiment, data_type, file):
    # Initialize the file processor with the directory
    processor = ExperimentFileProcessor(directory)
    
    # Process files for the specified experiment and data type
    extracted_data = processor.process_files(experiment, data_type, file)
    
    # Print the extracted information
    print(f"Extracted data for {experiment} ({data_type}):")
    return extracted_data


def process_file_arr_allgroundTruth(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file):
    
    log_file = directory+'/corr.txt'
    with open(log_file, 'a') as f:
        for file in file_arr:
            ground_truth_arr = []
            extracted_data = process_experiment_data(directory, f"exp{experiment}", exp_data, file)
            
            if extracted_data:
                # if extracted_data['result_aln_replica'][0] == '0':
                #         ground_truth_arr.append(groundTruth_file["PB_sample1"])
                #         ground_truth_arr.append(groundTruth_file["ill_sample1"])
                #         ground_truth_arr.append(groundTruth_file["Ground_truth1"])
                # if extracted_data['result_aln_replica'][0] == '2':
                #         ground_truth_arr.append(groundTruth_file["PB_sample2"])
                #         ground_truth_arr.append(groundTruth_file["ill_sample2"])
                #         ground_truth_arr.append(groundTruth_file["Ground_truth2"])
                if extracted_data['result_aln_replica'][0] == '0':
                        
                        ground_truth_arr.append(groundTruth_file["Ground_truth1"])
                if extracted_data['result_aln_replica'][0] == '2':
                        
                        ground_truth_arr.append(groundTruth_file["Ground_truth2"])
                for ground_truth in ground_truth_arr:
                    file_name = ground_truth
                    
                    ground_truth = os.path.join(groundTruth_main_dir, ground_truth)
                    predicted_result = os.path.join(directory, file)
                    if file_name.rsplit('.', 1)[0].rsplit('_')[0] == 'PB':
                        spearman_corr, pearson_corr =result_process.spearman_pearson_corr_generic(file_path1=predicted_result, file_path2=ground_truth, sep_by=",")
                    else:
                        spearman_corr, pearson_corr = result_process.spearman_pearson_corr_generic(file_path1=predicted_result, file_path2=ground_truth)
                    print(f"ground_truth {file_name}")

                    part1, part2 = result_process.format_file_name(predicted_result, ground_truth)
                    formatted_output = (

                        f"{part1} and {part2}.\nSpearman correlation: {spearman_corr:.3f}\nPearson correlation: {pearson_corr:.3f}\nground_truth {file_name}"
                    )
                    f.write(formatted_output + '\n')
            else:
                continue
            
            
        
def process_file_arr(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file):
     
     for file in file_arr:
        extracted_data = process_experiment_data(directory, f"exp{experiment}", exp_data, file)
        if extracted_data:
            if extracted_data['result_length'] == 'long':
                if extracted_data['result_aln_replica'][0] == '0':
                    ground_truth = groundTruth_file["PB_sample1"]
                if extracted_data['result_aln_replica'][0] == '2':
                    ground_truth = groundTruth_file["PB_sample2"]
            elif extracted_data['result_length'] == 'short':
                if extracted_data['result_aln_replica'][0] == '0':
                    ground_truth = groundTruth_file["ill_sample1"]
                if extracted_data['result_aln_replica'][0] == '2':
                    ground_truth = groundTruth_file["ill_sample2"]
            
            ground_truth = os.path.join(groundTruth_main_dir, ground_truth)
            predicted_result = os.path.join(directory, file)
            if extracted_data['result_length'] == 'long':
                result_process.spearman_pearson_corr_generic(file_path1=predicted_result, file_path2=ground_truth, sep_by=",")
            else:
                result_process.spearman_pearson_corr_generic(file_path1=predicted_result, file_path2=ground_truth)
        else:
            continue
        
        
        
   
if __name__ == "__main__":
    
    ######## parameters ##########
    main_result_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results'
    experiment_file = 'exprmntSingleRun_2024_00_00__00_00_00'
    simulation = 1
    experiment = 1  # "Different experiment setup, 1: for 1 sample, 2 for merged, 4 for multisample, 5 for merged multisample"
    # file_arr  = ['output_Simulation_VIGD_token_10896541_sample1_file1_ds100num1aln01long_GDlr_0.01_AlphaInitial_1.0_EMround_25_2024_12_8_22_22_08.tsv']
    directory = os.path.join(main_result_dir, experiment_file, 'files/output_files/')
    file_arr  = os.listdir(directory)
    groundTruth_main_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths'
#     groundTruth_file = {
#     "ill_sample1": "ill_sample1_gt.tsv",
#     "ill_sample2": "ill_sample2_gt.tsv",
#     "PB_sample1": "PB_sample1_gt.tsv",
#     "PB_sample2": "PB_sample1_gt.tsv",  
#     "Ground_truth1": "sample1_gt.tsv",
#     "Ground_truth2": "sample2_gt.tsv",
# }
    groundTruth_file = {  
    "Ground_truth1": "sample1_gt.tsv",
    "Ground_truth2": "sample2_gt.tsv",
}
    ######## parameters ##########

    print(f"experiment {experiment}")
    # Call the function with experiment name and data type
    # directory = os.path.join(main_result_dir, experiment_file, 'files/output_files/')
    
    if simulation:
        exp_data = 'simulation'
    else:
        exp_data = 'real'

    process_file_arr_allgroundTruth(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file)
    # process_file_arr(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file)