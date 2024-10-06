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


main_directory = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/round11/'
ground_truth_file = 'SRspCorr_0.99_LRspCorr_0.84_day_0_isoformAbundance_groundTruth.pkl'
# sample_ground_truth = '2'
with open(main_directory+ground_truth_file, 'rb') as f:
    theta_Ground_truth = pickle.load(f)

experiment_file = 'exprmnt_2024_09_29__22_02_56'
main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/new_simulation_ActualResults2'
directory = os.path.join(main_dir, experiment_file, 'files/output_files/')

file_pattern = re.compile(
    r'output_Simulation_VIGD_token_20428844_sample(\d+)_.*'
)

theta_file_list = []
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
            print(f"predicted_theta {file}\nground_truth {ground_truth_file}\nSpearman_corr {spearman_corr}\nPearson_corr {pearson_corr} ")


            print()

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