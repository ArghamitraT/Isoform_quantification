import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
import pandas as pd
import datetime
import time
import os
import re
from collections import defaultdict


def create_image_name(name, format=".png"):

    crnt_tm = datetime.datetime.now()
    image_name = (name+"_" + str(crnt_tm.year) + "_" + str(crnt_tm.month) + "_" + str(crnt_tm.day) + "_"
                  + time.strftime("%H_%M_%S") + format)
    return image_name

def plot_EM_results(alpha_history, convergence_history, theta_history):
    # Plot convergence over iterations
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.plot(convergence_history)
    plt.title('Convergence Over Iterations')
    plt.xlabel('Iteration')
    plt.ylabel('Convergence')

    # Plot alpha values over iterations
    plt.subplot(1, 3, 2)
    plt.plot(alpha_history)
    plt.title('Alpha Values Over Iterations')
    plt.xlabel('Iteration')
    plt.ylabel('Alpha Value')
    plt.legend()
    plt.savefig("figures/"+create_image_name("convergence_alpha"))


    for sample_key in theta_history:
        # Determine the consistent order of isoforms across all iterations
        isoforms = list(theta_history[sample_key][0].keys())  # Assuming the first Counter has all isoforms

        # Convert the Counters to a NumPy array
        theta_values = np.array([
            [counter[isoform] for isoform in isoforms]
            for counter in theta_history['sample1']
        ])

        # Now you can plot each column of theta_values across iterations
        plt.figure()
        for i, isoform in enumerate(isoforms):
            plt.plot(theta_values[:, i], label=isoform)
        plt.legend()
        plt.title(f'Isoform Values Across Iterations for {sample_key}')
        plt.xlabel('Iteration')
        plt.ylabel('Value')
        plt.savefig("figures/" + create_image_name("theta_"+sample_key))
        plt.show()


def fraction_to_float(fraction_str):
    # Strip the input to remove any leading/trailing whitespace
    fraction_str = fraction_str.strip()

    # Try to convert the string directly to a float (works if there is no fraction)
    try:
        return float(fraction_str)
    except ValueError:
        # Split the string on the '/'
        fraction_parts = fraction_str.split('/')
        if len(fraction_parts) == 2:  # It's a fraction
            numerator = float(fraction_parts[0])
            denominator = float(fraction_parts[1])
            return numerator / denominator
        else:
            raise ValueError(f"Invalid fraction string: {fraction_str}")


def spearman_corr_SIRV(file_path1, sample):

    ground_truth = pd.read_csv('../../data/SIRV/SIRV_Shree/E2_molarity.csv')
    # Reset the index of our_quant so that the isoform names become a column
    our_quant = pd.read_csv(file_path1, sep="\t")
    our_quant_reset = our_quant.reset_index()

    # Rename the columns accordingly
    #our_quant_reset.columns = ['transcript_name', 'other_column1', 'est_count', 'tpm']

    # Clean the isoform names in our_quant to match the naming convention in df2
    our_quant_reset['cleaned_name'] = our_quant_reset['transcript_name'].str.replace(r'\(\+\)|\(\-\)', '', regex=True)

    # Initialize lists to store matched TPM and molarity values
    matched_tpm = []
    matched_molarity = []

    # Iterate over our_quant to find matching molarity values in ground_truth
    for index, row in our_quant_reset.iterrows():
        # Extract the cleaned isoform name and tpm value
        cleaned_name = row['cleaned_name']
        tpm_value = row['tpm']

        # Find the corresponding molarity value in ground_truth (assuming 'E2' column has the molarity)
        molarity_value = ground_truth.loc[ground_truth['Name'] == cleaned_name, 'E2'].values

        # If a matching isoform is found in ground_truth
        if len(molarity_value) > 0:
            # Append the tpm and molarity to the respective lists
            matched_tpm.append(tpm_value)
            matched_molarity.append(fraction_to_float(molarity_value[0]))  # Take the first match in case of multiple

    # Calculate Spearman's correlation using the matched lists
    correlation, p_value = spearmanr(matched_tpm, matched_molarity)

    # Output the results
    print(f'Spearman correlation coefficient {sample}: {correlation}')
    print(f'P-value {sample}: {p_value}')

def fraction_to_float_gen(value):
    try:
        return float(value)
    except ValueError:
        return None

def split_transcript_name(transcript):
    return transcript.split('|')[0]

def spearman_corr_generic(file_path1, file_path2, directory):
    def split_transcript_name(name):
        # Add your implementation for split_transcript_name function
        return name

    def fraction_to_float_gen(value):
        # Add your implementation for fraction_to_float_gen function
        return value

    # Load the datasets
    our_quant = pd.read_csv(file_path1, sep="\t")
    ground_truth = pd.read_csv(file_path2, sep="\t")

    # Apply the function to the 'transcript_name' column of your DataFrame
    our_quant['transcript_name'] = our_quant['transcript_name'].apply(split_transcript_name)
    ground_truth['transcript_name'] = ground_truth['transcript_name'].apply(split_transcript_name)

    # Find common isoforms
    common_isoforms = pd.merge(our_quant, ground_truth, on='transcript_name', suffixes=('_quant', '_truth'))

    # Ensure both datasets have TPM values in the respective columns named accurately
    if 'tpm_quant' in common_isoforms and 'tpm_truth' in common_isoforms:
        # Convert any fractional notation in TPM values from ground_truth to float if necessary
        common_isoforms['tpm_truth'] = common_isoforms['tpm_truth'].apply(fraction_to_float_gen)

        # Calculate Spearman's correlation using the matched TPM lists
        correlation, p_value = spearmanr(common_isoforms['tpm_quant'], common_isoforms['tpm_truth'])

        # Extract the relevant parts from file names
        file_name1 = file_path1.split('/')[-1]
        file_name2 = file_path2.split('/')[-1]
        part1 = "_".join(file_name1.split('_')[5:10])
        part2 = "_".join(file_name2.split('_')[5:8])

        formatted_output = (
            f"Spearman correlation: {part1} and {part2} is {correlation:.6f}"
        )
        # Output the results
        timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
        output_file = directory+f'corr_results{timestamp}.txt'
        with open(output_file, 'a') as f:
            f.write(formatted_output + '\n')

        print(formatted_output)
        print(f"P-value for correlation: {p_value}")
    else:
        print("TPM columns missing or incorrectly named in one of the datasets.")


# SIRV_output = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_2024_07_15__14_27_36/files/output_files/output_SIRV_VIGD_15947722_sample2_aln_E2_2024_7_15_14_30_45.tsv'
# spearman_corr_SIRV(SIRV_output, sample='2')

# output_file = os.path.join(os.getcwd(), '../../files/results/single_run/')
# SIRV_output = output_file+'output_SIRV_VIGD_sample2_aln_E2_2024_7_15_15_17_56.tsv'
# spearman_corr_SIRV(SIRV_output, sample='2')

def pair_files(directory, criteria='replica'):
    # Dictionary to store the files based on pairing criteria
    file_pairs = defaultdict(list)

    # Regular expression to match the files
    file_pattern = re.compile(r'(.*)_(\d{2})_(long|short)(.*)')

    # List all files in the directory
    for file in os.listdir(directory):
        match = file_pattern.match(file)
        if match:
            prefix, day_replica, length, suffix = match.groups()
            if criteria == 'day':
                key = day_replica[:1]  # Use day (first digit)
            elif criteria == 'replica':
                key = day_replica[1]  # Use replica (second digit)
            elif criteria == 'alignment':
                key = prefix + suffix  # Use prefix and suffix (alignment)
            else:
                raise ValueError("Criteria must be 'day', 'replica', or 'alignment'")
            file_pairs[key].append(file)

    # Create the pairs
    paired_files = []
    for key, files in file_pairs.items():
        long_files = [f for f in files if 'long' in f]
        short_files = [f for f in files if 'short' in f]
        for long_file in long_files:
            for short_file in short_files:
                paired_files.append((long_file, short_file))

    return paired_files


# Example usage
if __name__ == "__main__":
    experiment_file = 'exprmnt_2024_07_15__15_38_55'
    main_dir = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/files/results/'
    directory = os.path.join(main_dir, experiment_file, 'files/output_files/')
    criteria = 'replica'  # Change to 'day', 'replica', or 'alignment' as needed
    paired_files = pair_files(directory, criteria)

    for pair in paired_files:
        #print(pair)
        spearman_corr_generic(directory+pair[0], directory+pair[1], directory)