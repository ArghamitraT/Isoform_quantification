import pandas as pd
import os
import generate_result_stat as result_process
from util import ExperimentFileProcessor
import matplotlib.pyplot as plt

def plot_predictedVSgroundtruth_theta(common_isoforms, x_axis, y_axis, figure_path):

    predicted_theta = common_isoforms['tpm_quant']
    ground_truth_theta = common_isoforms['tpm_truth']

    # Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(predicted_theta, ground_truth_theta, alpha=0.7, label='Isoforms')
    plt.plot([0, 1], [0, 1], color='red', linestyle='--', label='Ideal Agreement')  # Add diagonal line
    plt.title('Isoform Estimate vs Isoform Truth')
    plt.xlabel(f'{x_axis}(Predicted)')
    plt.ylabel(f'{y_axis}(True)')
    plt.legend()
    plt.grid(True)

    # Save or show plot
    os.makedirs(figure_path, exist_ok=True)
    plt.savefig(os.path.join(figure_path, f'{x_axis}AND{y_axis}isoform_estimate_vs_truth.png'))  # Save as a file
    print('saved')
    # plt.show()



def process_experiment_data(directory, experiment, data_type, file):
    # Initialize the file processor with the directory
    processor = ExperimentFileProcessor(directory)
    
    # Process files for the specified experiment and data type
    extracted_data = processor.process_files(experiment, data_type, file)
    
    # Print the extracted information
    print(f"Extracted data for {experiment} ({data_type}):")
    return extracted_data

def process_file_arr(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file):
    
     for file in file_arr:
        ground_truth_arr = []
        extracted_data = process_experiment_data(directory, f"exp{experiment}", exp_data, file)
        if extracted_data:
            if extracted_data['result_aln_replica'][0] == '0':
                    ground_truth_arr.append(groundTruth_file["PB_sample1"])
                    ground_truth_arr.append(groundTruth_file["ill_sample1"])
                    ground_truth_arr.append(groundTruth_file["Ground_truth1"])
            if extracted_data['result_aln_replica'][0] == '2':
                    ground_truth_arr.append(groundTruth_file["PB_sample2"])
                    ground_truth_arr.append(groundTruth_file["ill_sample2"])
                    ground_truth_arr.append(groundTruth_file["Ground_truth2"])
            
            for ground_truth in ground_truth_arr:
                file_name = ground_truth    
                ground_truth = os.path.join(groundTruth_main_dir, ground_truth)
                predicted_result = os.path.join(directory, file)
                if file_name.rsplit('.', 1)[0].rsplit('_')[0] == 'PB':
                    common_isoforms = result_process.csv_tpm_processing(file_path1=predicted_result, file_path2=ground_truth, sep_by=",")
                else:
                    common_isoforms =result_process.csv_tpm_processing(file_path1=predicted_result, file_path2=ground_truth)
                
                if 'tpm_quant' in common_isoforms and 'tpm_truth' in common_isoforms:
                    common_isoforms['tpm_truth'] = common_isoforms['tpm_truth'].apply(result_process.fraction_to_float_gen)

                figure_path = '/'.join(directory.rstrip('/').split('/')[:-2])+'/figures/'
                plot_predictedVSgroundtruth_theta(common_isoforms, file, file_name, figure_path)
                print(f"ground_truth {file_name}")
        else:
            continue

def main():

    ######## parameters ##########
    main_result_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results'
    experiment_file = 'exprmnt_2024_12_08__22_18_02'
    simulation = 1
    experiment = 1  # "Different experiment setup, 1: for 1 sample, 2 for merged, 4 for multisample, 5 for merged multisample"
    groundTruth_main_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths'
    groundTruth_file = {
    "ill_sample1": "ill_sample1_gt.tsv",
    "ill_sample2": "ill_sample2_gt.tsv",
    "PB_sample1": "PB_sample1_gt.tsv",
    "PB_sample2": "PB_sample1_gt.tsv",  
    "Ground_truth1": "sample1_gt.tsv",
    "Ground_truth2": "sample2_gt.tsv",
}
    ######## parameters ##########


    predicted_theta_path_full = os.path.join(main_result_dir, experiment_file, 'files/output_files/')
    file_arr  = os.listdir(predicted_theta_path_full)
    # ground_truth_path_full = os.path.join(main_dir, ground_truth_path)

    if simulation:
        exp_data = 'simulation'
    else:
        exp_data = 'real'
    
    process_file_arr(file_arr, predicted_theta_path_full, experiment, exp_data, groundTruth_main_dir, groundTruth_file)
    # process_file_arr(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file)


if __name__ == "__main__":
    main()
