

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
            matched_tpm.append(np.log(tpm_value+1))
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


def csv_tpm_processing(file_path1, file_path2, suffixes=('_quant', '_truth')):
    # Load the datasets
    our_quant = pd.read_csv(file_path1, sep="\t")
    ground_truth = pd.read_csv(file_path2, sep="\t")

    # Apply the function to the 'transcript_name' column of your DataFrame
    our_quant['transcript_name'] = our_quant['transcript_name'].apply(split_transcript_name)
    ground_truth['transcript_name'] = ground_truth['transcript_name'].apply(split_transcript_name)

    # Find common isoforms
    common_isoforms = pd.merge(our_quant, ground_truth, on='transcript_name', suffixes=suffixes)

    return common_isoforms

def spearman_corr_generic(file_path1, file_path2, directory):

    common_isoforms = csv_tpm_processing(file_path1, file_path2)

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

output_file = os.path.join(os.getcwd(), '../../files/results/single_run/files/output_files/')
#output_file = os.path.join(os.getcwd(), '../../files/results/exprmnt_2024_07_28__17_38_26/files/output_files/')
SIRV_output = output_file+'output_SIRV_VIGD_00000_sample2_aln_E2_2024_8_5_22_39_32.tsv'
spearman_corr_SIRV(SIRV_output, sample='2')
print()

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
            elif criteria == 'length':
                key = length  # Use prefix and suffix (alignment)
            else:
                raise ValueError("Criteria must be 'day', 'replica', or 'alignment'")
            file_pairs[key].append(file)

    # Create the pairs
    paired_files = []
    if criteria == 'replica':
        for key, files in file_pairs.items():
            long_files = [f for f in files if 'long' in f]
            short_files = [f for f in files if 'short' in f]
            for long_file in long_files:
                for short_file in short_files:
                    paired_files.append((long_file, short_file))
    elif criteria == 'length':
        for key, files in file_pairs.items():
            paired_files.append((files[0], files[1]))
        print()

    return paired_files

# Function to calculate CV
def calculate_cv(data, type):
    if type == 'long':
        colm_name = 'long'
    else:
        colm_name = 'short'
    data_log = np.log(data[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']] + 1)
    data['mean_abundance'] = data_log[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']].mean(axis=1)
    data['std_abundance'] = data_log[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']].std(axis=1)
    CV_ig  = data['std_abundance'] / data['mean_abundance']
    data['CV'] = CV_ig

    CV_ig_squared = CV_ig ** 2
    IM = np.sqrt(CV_ig_squared.mean())

    # Calculate ACVC (assuming equal intervals in abundance for simplicity)
    sorted_u_ig = data['mean_abundance'].sort_values()
    sorted_CV_ig = data['CV'].loc[sorted_u_ig.index]

    ACVC = np.trapz(sorted_CV_ig, x=sorted_u_ig)
    print(f"IM {IM}, ACVC {ACVC}")

    return data

def calculate_acvc(cv_values, abundance):
    return trapz(cv_values, abundance)

def calculate_im_acvc(rep1, rep2, directory, type='long'):

    if type == 'long':
        rep1_data = csv_tpm_processing(directory+rep1, directory+rep2, suffixes=('_longRep1', '_longRep2'))
    else:
        rep1_data = csv_tpm_processing(directory+rep1, directory+rep2, suffixes=('_shortRep1', '_shortRep2'))

    # Calculate CV for each replicate
    rep1_data = calculate_cv(rep1_data, type=type)


    # Sort data by mean abundance for plotting
    rep1_data_sorted = rep1_data.sort_values('mean_abundance')

    # Box Plot for IM
    plt.figure(figsize=(8, 6))
    sns.boxplot(y='CV', data=rep1_data_sorted, palette="Set2", boxprops=dict(alpha=0.5))
    #sns.stripplot(y='CV', data=rep1_data_sorted, jitter=True, linewidth=1, palette="Set2")
    plt.title('Box Plot with Scatter Overlay for IM')
    plt.ylabel('IM')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.ylim(0, 1.5)  # Adjust the y-axis limit based on your data range
    plt.show()
    # Plot
    plt.plot(rep1_data['mean_abundance'], rep1_data['CV'], 'o')
    plt.xlabel('Transcript abundance (log2(TPM+1))')
    plt.ylabel('CV')
    plt.title('CV vs Transcript Abundance')
    plt.ylim(0, 1.5)
    plt.xlim(0, 10)
    plt.legend()
    plt.show()
    print()

    # Log2-transform abundance estimates
    data_log2 = np.log2(rep1_data_sorted[[f'tpm_{type}Rep1', f'tpm_{type}Rep2']] - 1)

    # Calculate mean and standard deviation
    mean_abundance = data_log2.mean(axis=1)
    std_abundance = data_log2.std(axis=1)

    # Calculate CV
    CV = std_abundance / mean_abundance

    # Plot
    plt.plot(mean_abundance, CV, 'o')
    plt.xlabel('Transcript abundance (log2(TPM-1))')
    plt.ylabel('CV')
    plt.title('CV vs Transcript Abundance')
    plt.ylim(0, 1.5)
    plt.xlim(0, 10)
    plt.legend()
    plt.show()
    print()

    # # Plot
    # plt.plot(mean_abundance, rep1_data['CV'], 'o')
    # plt.xlabel('Transcript abundance (log2(TPM-1))')
    # plt.ylabel('CV')
    # plt.title('CV vs Transcript Abundance')
    # plt.ylim(0, 1.5)
    # plt.xlim(0, 10)
    # plt.legend()
    # plt.show()
    # print()

    # # Log2-transform abundance estimates
    # data_log2 = np.log2(rep1_data_sorted[[f'tpm_{type}Rep1', f'tpm_{type}Rep2']] - 1)
    #
    # # Calculate mean and standard deviation
    # mean_abundance = data_log2.mean(axis=1)
    # std_abundance = data_log2.std(axis=1)
    #
    # # Calculate CV
    # CV = std_abundance / mean_abundance



    # # Smooth the CV curves using a spline
    # rep1_spline = UnivariateSpline(np.log(rep1_data_sorted['mean_abundance']), rep1_data_sorted['CV'], s=1)
    #
    # # Generate smooth CV values for plotting
    # rep1_smooth_CV = rep1_spline(np.log(rep1_data_sorted['mean_abundance']))
    #
    # # Plot the CV curves
    # plt.figure(figsize=(10, 6))
    # plt.plot(rep1_data_sorted['mean_abundance'], rep1_smooth_CV, label='LR Rep 1 and 2 CV', color='blue')
    # #plt.plot(rep2_data_sorted['mean_abundance'], rep2_smooth_CV, label='Replicate 2 CV', color='red')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel('Transcript Abundance (TPM)')
    # plt.ylabel('Coefficient of Variation (CV)')
    # plt.title('CV vs. Transcript Abundance for Long Read and Short Read Combined')
    # plt.legend()
    # plt.grid(True)
    # plt.show()


def merge_csv_files(file1, file2, output_dir):
    # Load the first CSV file
    df1 = pd.read_csv(file1, sep="\t")

    # Load the second CSV file
    df2 = pd.read_csv(file2, sep="\t")

    # Identify common transcripts
    # Find common isoforms
    # common_isoforms = pd.merge(df1, df2, on='transcript_name')

    common_transcripts = pd.merge(df1, df2, left_on='transcript_name', right_on='transcript_name',
                                  suffixes=('_sample1', '_sample2'))

    # Select relevant columns
    merged_df = common_transcripts[['transcript_name', 'tpm_sample1', 'tpm_sample2']]

    # Rename columns to match the desired format
    merged_df.columns = ['ID', 'sample1', 'sample2']

    # Save the merged dataframe to a new CSV file
    output_file = output_dir+'merged'+time.strftime("_%Y_%m_%d__%H_%M_%S")+'.tsv'
    merged_df.to_csv(output_file, index=False)

    print(f'Merged file saved as {output_file}')


def main():
    experiment_file = 'exprmnt_2024_09_02__17_26_30' #the name of the folder where you have your .csv files for isoform quantifications
    main_dir = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/files/results/'
    directory = os.path.join(main_dir, experiment_file, 'files/output_files/')
    criteria = 'length'  # Change to 'day', 'replica', or 'length' as needed
    paired_files = pair_files(directory, criteria)

    # calculate_im_acvc(paired_files[0], paired_files[1], directory, type='short')
    # spearman_corr_generic(directory + paired_files[0][0], directory + paired_files[1][0], directory)
    # spearman_corr_generic(directory + paired_files[0][1], directory + paired_files[1][1], directory)

    i=0
    for pair in paired_files:
        #merge_csv_files(directory + pair[0], directory + pair[1], directory)

        print(pair)
        i+=1

        if i ==1:
            calculate_im_acvc(pair[0], pair[1], directory, type='long')
            spearman_corr_generic(directory + pair[0], directory + pair[1], directory)
        elif i ==2:
            calculate_im_acvc(pair[0], pair[1], directory, type='short')
            spearman_corr_generic(directory+pair[0], directory+pair[1], directory)


# Example usage
if __name__ == "__main__":
    main()
