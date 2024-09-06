import os
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
import pandas as pd
import datetime
import time
import os
import re
from collections import defaultdict
from scipy.interpolate import UnivariateSpline
#import seaborn as sns
from scipy.integrate import trapz
import sys
import generate_stat as stat_formula
import csv

# Function to find log file paths
def find_log_file_paths(main_dir, experiment_file):
    final_result_dir = os.path.join(main_dir, experiment_file, 'files/output_files/') 
    #final_result_dir = os.path.join(main_dir, experiment_file, 'files/problem/') #(AT)
    log_file_paths = []
    for file in os.listdir(final_result_dir):
        if file.startswith('out_main_') and not file.endswith('.csv'):
            log_file_paths.append(os.path.join(final_result_dir, file))
    return log_file_paths

def parse_log_file(log_file_path, multi_occurrence_vars, single_occurrence_vars):
    with open(log_file_path, 'r') as file:
        lines = file.readlines()

    # Initialize dictionaries to store regular expressions and corresponding data lists
    data = {var: [] for var in multi_occurrence_vars}
    single_occurrence_data = {var: None for var in single_occurrence_vars}

    # Update the regex to handle numbers or the word 'inf'
    regex_patterns_multi = {var: re.compile(f'{var} (-?\\d*\\.?\\d+|\\d+|inf|-inf)') for var in multi_occurrence_vars}
    regex_patterns_single = {var: re.compile(f'{var} (-?\\d*\\.?\\d+|\\d+|inf|-inf|/[^\\s]+)') for var in single_occurrence_vars}

    # Parse the file line by line
    for line in lines:
        for var, regex in regex_patterns_multi.items():
            match = regex.search(line)
            if match:
                value = match.group(1)
                # Convert the value to float or inf as necessary
                data[var].append(float(value) if value not in ['inf', '-inf'] else float(value))

        for var, regex in regex_patterns_single.items():
            match = regex.search(line)
            if match and single_occurrence_data[var] is None:
                value = match.group(1)
                # Check if the value is a file path
                if '/' in value:
                    # Extract the file name from the path
                    single_occurrence_data[var] = value.split('/')[-1]
                else:
                    # Convert the value to float if it's not a file path
                    single_occurrence_data[var] = float(value) if value not in ['inf', '-inf'] else float(value)


    return single_occurrence_data


def find_tsv_files(directory, paired_codes):
    paired_tsv_files = []
    for code1, code2 in paired_codes:
        file1_long = None
        file2_long = None
        file1_short = None
        file2_short = None
        
        for file in os.listdir(directory):
            if file.endswith('.tsv'):
                if code1 in file and 'long' in file:
                    file1_long = file
                elif code2 in file and 'long' in file:
                    file2_long = file
                elif code1 in file and 'short' in file:
                    file1_short = file
                elif code2 in file and 'short' in file:
                    file2_short = file
                
                # Pair long files
                if file1_long and file2_long:
                    paired_tsv_files.append((file1_long, file2_long))
                    file1_long = file2_long = None  # Reset after pairing
                    
                # Pair short files
                if file1_short and file2_short:
                    paired_tsv_files.append((file1_short, file2_short))
                    file1_short = file2_short = None  # Reset after pairing

    return paired_tsv_files

 
def main():
    experiment_file = 'exprmnt_2024_08_10__02_05_36'
    main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/'
    directory = os.path.join(main_dir, experiment_file, 'files/output_files/')

    # Find all log file paths
    log_file_paths = find_log_file_paths(main_dir, experiment_file)
    if not log_file_paths:
        print("No file found matching the format 'out_main_'")
        return

    # Dictionary to store the codes based on alpha, ds, length, and replica
    file_info = defaultdict(list)

    # Parse each log file and store the relevant data
    # Parse each log file and store the relevant data
    for log_file_path in log_file_paths:
        single_occurrence_vars = ['alpha_initial', 'sample_1', 'sample_2']
        single_occurrence_data = parse_log_file(log_file_path, [], single_occurrence_vars)
        
        # Extract the code from the log file name
        log_code = re.search(r'out_main_(\d+)', log_file_path).group(1)

        def extract_info(sample, log_code, alpha):
            ds_match = re.search(r'ds_(\d+)', sample)
            aln_match = re.search(r'aln_(\d{2})', sample)
            length_match = re.search(r'(long|short)', sample)
            if ds_match and aln_match and length_match:
                ds = ds_match.group(1)
                aln = aln_match.group(1)
                length = length_match.group(0)
                if aln == '03' or length == 'short':
                    return
                key = (alpha, ds, length)
                file_info[key].append((log_code, aln))
        
        alpha_initial = single_occurrence_data['alpha_initial']
        try:
            extract_info(single_occurrence_data['sample_1'], log_code, alpha_initial)
            extract_info(single_occurrence_data['sample_2'], log_code, alpha_initial)
        except:
            print(f"problem: {log_file_path}")

    
    paired_codes = []
    for key, codes in file_info.items():
        # Sort codes by replica (aln) to ensure correct pairing
        codes.sort(key=lambda x: x[1])
        
        # Create all possible pairs within the same key
        for i in range(len(codes)):
            for j in range(i + 1, len(codes)):
                if codes[i][1] != codes[j][1]:  # Ensure different replicas
                    paired_codes.append((key, codes[i][0], codes[j][0]))

    # Find and pair the corresponding .tsv files
    paired_files = find_tsv_files(directory, [(pair[1], pair[2]) for pair in paired_codes])


    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    # log_file = directory + f'corr_results{timestamp}.txt'
    # sys.stdout = open(log_file, 'w')

    # for idx, pair in enumerate(paired_files):
    #     # Extract the corresponding key (which contains alpha, ds, length)
    #     for key, codes in file_info.items():
    #         for code, _ in codes:
    #             if code == pair[0].split('_')[3]:
    #                 alpha, ds, length = key
        

    #     print(f"Alpha: {alpha}, DS: {ds}, Length: {length}, Paired TSV Files: {pair[0]} and {pair[1]}")

    #     """ ## CALL 'spearman_corr_generic' FUNC FIRST ## """
    #     stat_formula.spearman_corr_generic(directory + pair[0], directory + pair[1], log_file)
    #     stat_formula.calculate_im_acvc(pair[0], pair[1], directory, log_file)

    # sys.stdout.close()
    # Prepare the CSV file for writing
    csv_file = directory + f'corr_results{timestamp}.csv'

    # Open the CSV file for writing
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # Write the header row
        writer.writerow(['Alpha', 'DS', 'Length', 'File1', 'File2', 'Code1', 'Code2', 'Spearman_Corr', 'ACVC', 'IM'])
        
        # Iterate through the paired files
        for idx, pair in enumerate(paired_files):
            # Extract the corresponding key (which contains alpha, ds, length)
            alpha, ds, length = None, None, None
            code1, code2 = pair[0].split('_')[3], pair[1].split('_')[3]
            
            for key, codes in file_info.items():
                for code, _ in codes:
                    if code == code1:
                        alpha, ds, length = key
                        break
                if alpha:  # Break outer loop if key is found
                    break
            
            # Perform calculations
            spearman_corr = stat_formula.spearman_corr_generic(directory + pair[0], directory + pair[1])
            ACVC, IM = stat_formula.calculate_im_acvc(pair[0], pair[1], directory)
            
            # Write the row to the CSV file
            writer.writerow([alpha, ds, length, "_".join(pair[0].split('_')[3:11]), "_".join(pair[1].split('_')[3:11]), code1, code2, spearman_corr, ACVC, IM])

     
# Run the main function
if __name__ == "__main__":
    main()