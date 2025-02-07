"""
This script compares the simulation results with ground truth
"""
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr
import pandas as pd
import datetime
import time
import os
import re
from collections import defaultdict
from scipy.interpolate import UnivariateSpline
# import seaborn as sns
from scipy.integrate import trapz
import sys
import csv
import pickle


def find_correlation_W_grndTurth_exp2(directory, ground_truth_file, theta_Ground_truth, file_pattern):

    for file in os.listdir(directory):
            match = file_pattern.search(file)
            if match:
                # Parse the matched groups
                #if sample == sample_ground_truth:
                our_quant = pd.read_csv(directory+file, sep="\t")
                quant_values = our_quant.set_index('transcript_name')['raw']*1000000
                # Get the theta values for the corresponding transcript names in theta_Ground_truth
                theta_values = pd.Series(theta_Ground_truth[f'sample_1'])
                aligned_quant, aligned_theta = quant_values.align(theta_values, join='inner')
                spearman_corr, _ = spearmanr(aligned_quant, aligned_theta)
                pearson_corr, p_value = pearsonr(aligned_quant, aligned_theta)
                # print(f"predicted_theta {file}\nground_truth {ground_truth_file}\nSpearman_corr {spearman_corr}\nPearson_corr {pearson_corr} ")
                print(f"Spearman_corr {spearman_corr:.3f}\nPearson_corr {pearson_corr:.3f} ")


                print()

def find_correlation_W_grndTurth_exp5(directory, ground_truth_file, theta_Ground_truth, file_pattern, sample_ground_truth):

    theta_file_list = []
    for file in os.listdir(directory):
            match = file_pattern.search(file)
            if match:
                # Parse the matched groups
                (sample) = match.groups()

                if sample[0] == sample_ground_truth:
                    our_quant = pd.read_csv(directory+file, sep="\t")
                    quant_values = our_quant.set_index('transcript_name')['raw']*1000000
                    # Get the theta values for the corresponding transcript names in theta_Ground_truth
                    theta_values = pd.Series(theta_Ground_truth[f'sample_1'])
                    aligned_quant, aligned_theta = quant_values.align(theta_values, join='inner')
                    spearman_corr, _ = spearmanr(aligned_quant, aligned_theta)
                    pearson_corr, p_value = pearsonr(aligned_quant, aligned_theta)
                    print(f"predicted_theta {file}\nground_truth {ground_truth_file}\nSpearman_corr {spearman_corr}\nPearson_corr {pearson_corr} ")


                    print()

def compute_groupwise_correlations_exp2(directory, ground_truth_file, theta_Ground_truth, file_pattern):
    """
    Computes the Spearman and Pearson correlations between predicted theta values and ground truth theta values
    for three different groups (top, middle, bottom) based on the sorted ground truth values.
    
    Parameters:
    - directory: Directory containing the predicted theta files.
    - theta_Ground_truth: DataFrame containing the ground truth theta values with a column matching `sample_key`.
    - file_pattern: Compiled regex pattern to filter files in the directory.
    - sample_key: The key for the sample column in `theta_Ground_truth` to use for ground truth values.
    
    Returns:
    None. Prints the correlation results for each group.
    """
    
    sample_key='sample_1'
    # Extract ground truth theta values and sort them
    theta_values = pd.Series(theta_Ground_truth[sample_key])
    sorted_theta = dict(sorted(theta_values.items(), key=lambda item: item[1]))
    # Convert to a list of items for slicing
    sorted_theta_items = list(sorted_theta.items())
    # Determine indices to split into three equal groups
    split_indices = [int(len(sorted_theta) / 3), int(2 * len(sorted_theta) / 3)]
    # Slice the list into three groups
    groups = {
        'top_1_3': sorted_theta_items[:split_indices[0]],
        'mid_2_3': sorted_theta_items[split_indices[0]:split_indices[1]],
        'bottom_3_3': sorted_theta_items[split_indices[1]:]
    }

    
    # Iterate over files in the directory and compute correlations for each group
    for file in os.listdir(directory):
        match = file_pattern.search(file)
        if match:
            # Read the predicted theta values from the file
            our_quant = pd.read_csv(os.path.join(directory, file), sep="\t")
            quant_values = our_quant.set_index('transcript_name')['raw'] * 1000000

            print(f"File: {file}")
            print(f"ground_truth: {ground_truth_file}")

            # Iterate over each group to calculate correlation
            for group_name, group_indices in groups.items():
                # Convert the list of tuples into a pandas Series
                group_indices_series = pd.Series(
                    data=[value for _, value in group_indices],
                    index=[isoform_name for isoform_name, _ in group_indices])

                # Align quant_values with group_indices
                aligned_quant, aligned_theta = quant_values.align(group_indices_series, join='inner')

                # # Align predicted theta with the ground truth for the current group
                # aligned_quant, aligned_theta = quant_values.align(group_indices, join='inner')

                # Calculate correlations
                spearman_corr, _ = spearmanr(aligned_quant, aligned_theta)
                pearson_corr, _ = pearsonr(aligned_quant, aligned_theta)

                # Print results for each group
                
                print(f"Ground Truth Group (ascending): {group_name}")
                print(f"Spearman Correlation: {spearman_corr:.3f}")
                print(f"Pearson Correlation: {pearson_corr:.3f}")
                print()


def main():
    
    file_pattern = re.compile(
        r'output_Simulation_VIGD_token_5851266_sample(\d+)_.*'
    ) 
    ground_truth_file = 'SRspCorr_0.99_LRspCorr_0.84_day_0_isoformAbundance_groundTruth.pkl'
    experiment_file = 'exprmnt_2024_10_12__00_43_44'
    experiment = 2
    if experiment ==5:
         sample_ground_truth = '1'

    
    
    main_directory = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/round11/'
    # sample_ground_truth = '2'
    with open(main_directory+ground_truth_file, 'rb') as f:
        theta_Ground_truth = pickle.load(f)
    experiment_main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/new_simulation_ActualResults3'
    directory = os.path.join(experiment_main_dir, experiment_file, 'files/output_files/')
    
    if experiment == 2:
        compute_groupwise_correlations_exp2(directory, ground_truth_file, theta_Ground_truth, file_pattern)
        find_correlation_W_grndTurth_exp2(directory, ground_truth_file, theta_Ground_truth, file_pattern)
    elif experiment == 5:
         find_correlation_W_grndTurth_exp5(directory, ground_truth_file, theta_Ground_truth, file_pattern, sample_ground_truth)

if __name__ == "__main__":
    main()


########## experiment 5 simulation comparison with ground truth #################

# main_directory = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/round11/'
# ground_truth_file = 'SRspCorr_0.99_LRspCorr_0.84_day_1_isoformAbundance_groundTruth.pkl'
# sample_ground_truth = '2'
# with open(main_directory+ground_truth_file, 'rb') as f:
#     theta_Ground_truth = pickle.load(f)

# experiment_file = 'exprmnt_2024_09_29__21_53_05'
# main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/new_simulation_ActualResults2'
# directory = os.path.join(main_dir, experiment_file, 'files/output_files/')

# file_pattern = re.compile(
#     r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_.*'
# )

# theta_file_list = []
# for file in os.listdir(directory):
#         match = file_pattern.search(file)
#         if match:
#             # Parse the matched groups
#             (token, sample) = match.groups()

#             if sample == sample_ground_truth:
#                   our_quant = pd.read_csv(directory+file, sep="\t")
#                   quant_values = our_quant.set_index('transcript_name')['raw']*1000000
#                   # Get the theta values for the corresponding transcript names in theta_Ground_truth
#                   theta_values = pd.Series(theta_Ground_truth[f'sample_1'])
#                   aligned_quant, aligned_theta = quant_values.align(theta_values, join='inner')
#                   spearman_corr, _ = spearmanr(aligned_quant, aligned_theta)
#                   pearson_corr, p_value = pearsonr(aligned_quant, aligned_theta)
#                   print(f"predicted_theta {file}\nground_truth {ground_truth_file}\nSpearman_corr {spearman_corr}\nPearson_corr {pearson_corr} ")


#                   print()