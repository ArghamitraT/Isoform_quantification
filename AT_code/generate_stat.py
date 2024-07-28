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
import seaborn as sns
from scipy.integrate import trapz
import sys


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

def format_file_name(file_path1, file_path2, type):
    file_name1 = file_path1.split('/')[-1]
    file_name2 = file_path2.split('/')[-1]
    if type == 'long':
        part1 = "_".join(file_name1.split('_')[3:10])
        part2 = "_".join(file_name2.split('_')[3:10])
    else:
        part1 = "_".join(file_name1.split('_')[3:8])
        part2 = "_".join(file_name2.split('_')[3:8])

    print()
    return part1, part2
def spearman_corr_generic(file_path1, file_path2, log_file, type):
    common_isoforms = csv_tpm_processing(file_path1, file_path2)

    if 'tpm_quant' in common_isoforms and 'tpm_truth' in common_isoforms:
        common_isoforms['tpm_truth'] = common_isoforms['tpm_truth'].apply(fraction_to_float_gen)
        correlation, p_value = spearmanr(common_isoforms['tpm_quant'], common_isoforms['tpm_truth'])

        part1, part2 = format_file_name(file_path1, file_path2, type)

        formatted_output = (

            f"{part1} and {part2}.\nSpearman correlation: {correlation:.6f}"
        )
        with open(log_file, 'a') as f:
            f.write(formatted_output + '\n')

        print(formatted_output)
        print(f"P-value for correlation: {p_value}")
    else:
        print("TPM columns missing or incorrectly named in one of the datasets.")


def pair_files(directory, type):
    # Dictionary to store the files based on downsampling percentage and length type
    long_file_pairs = defaultdict(list)
    short_file_pairs = defaultdict(list)

    # Regular expression to match the files
    file_pattern = re.compile(r'(.*)_ds_(\d+)_aln_(\d{2})_(long|short)(.*)|(.*)_aln_(\d{2})_(long|short)(.*)')

    # This will give you LR and SR file names with same replica, eg: lr_01_replica1+lr_01_replica2,
    if type == 'replica':
        # List all files in the directory
        for file in os.listdir(directory):
            match = file_pattern.match(file)
            if match:
                if match.groups()[1]:  # Matches pattern with downsampling
                    prefix, ds_percentage, aln_replica, length, suffix = match.groups()[:5]
                    key = (ds_percentage, length)
                    long_file_pairs[key].append((file, aln_replica))
                else:  # Matches pattern without downsampling
                    prefix, aln_replica, length, suffix = match.groups()[5:]
                    token = re.search(r'VIGD_(\d+)', prefix).group(1)
                    key = (token)
                    short_file_pairs[key].append((file, aln_replica))

        # Create long read pairs
        paired_files = []
        valid_tokens = set()
        for key, files in long_file_pairs.items():
            paired_files.append((files[0][0], files[1][0]))
            # (AT)
            sr_file1 = short_file_pairs[files[0][0].split('_')[3]][0][0]
            sr_file2 = short_file_pairs[files[1][0].split('_')[3]][0][0]
            paired_files.append((sr_file1, sr_file2))

    # This will give you LR and SR file names those were trained together, eg: lr_01_replica1+sr_01_replica2,
    elif type =='within_trainee':
        # List all files in the directory
        for file in os.listdir(directory):
            match = file_pattern.match(file)
            if match:
                if match.groups()[1]:  # Matches pattern with downsampling
                    prefix, ds_percentage, aln_replica, length, suffix = match.groups()[:5]
                    token = re.search(r'VIGD_(\d+)', prefix).group(1)
                    key = (token)
                    long_file_pairs[key].append((file, aln_replica))
                else:  # Matches pattern without downsampling
                    prefix, aln_replica, length, suffix = match.groups()[5:]
                    token = re.search(r'VIGD_(\d+)', prefix).group(1)
                    key = (token)
                    short_file_pairs[key].append((file, aln_replica))

        # Create long read pairs
        paired_files = []
        valid_tokens = set()
        for key, files in long_file_pairs.items():
            paired_files.append((files[0][0], short_file_pairs[files[0][0].split('_')[3]][0][0]))

    # (AT)
    # x = list(short_file_pairs.items())
    # paired_files.append((x[0][1][0][0], x[1][1][0][0]))


    return paired_files



# Function to calculate CV
def calculate_cv(data, type, log_file):
    colm_name = 'long' if type == 'long' else 'short'
    data_log = np.log(data[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']] + 1)
    data['mean_abundance'] = data_log[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']].mean(axis=1)
    data['std_abundance'] = data_log[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']].std(axis=1)
    CV_ig = data['std_abundance'] / data['mean_abundance']
    data['CV'] = CV_ig

    CV_ig_squared = CV_ig ** 2
    IM = np.sqrt(CV_ig_squared.mean())

    sorted_u_ig = data['mean_abundance'].sort_values()
    sorted_CV_ig = data['CV'].loc[sorted_u_ig.index]

    ACVC = np.trapz(sorted_CV_ig, x=sorted_u_ig)
    result = f"IM {IM}, ACVC {ACVC}\n"
    with open(log_file, 'a') as f:
        f.write(result)
    print(result)

    return data

def calculate_acvc(cv_values, abundance):
    return trapz(cv_values, abundance)

def calculate_im_acvc(rep1, rep2, directory, log_file, type='long'):

    if type == 'long':
        rep1_data = csv_tpm_processing(directory+rep1, directory+rep2, suffixes=('_longRep1', '_longRep2'))
    else:
        rep1_data = csv_tpm_processing(directory+rep1, directory+rep2, suffixes=('_shortRep1', '_shortRep2'))

    rep1_data = calculate_cv(rep1_data, type=type, log_file=log_file)
    rep1_data_sorted = rep1_data.sort_values('mean_abundance')

    part1, part2 = format_file_name(rep1, rep2, type)

    # Split the path into parts
    path_parts = directory.split('/')
    # Remove the last two parts ('files' and 'output_files')
    new_path_parts = path_parts[:-3]
    # Append 'figures' to the path
    new_path_parts.append('figures')
    # Join the parts to form the new path
    figure_dir = '/'.join(new_path_parts)

    fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    sns.boxplot(y='CV', data=rep1_data_sorted, palette="Set2", ax=axes[0], boxprops=dict(alpha=0.5))
    axes[0].set_title('Box Plot for CV')
    axes[0].set_ylabel('CV')
    axes[0].grid(True, which='both', linestyle='--', linewidth=0.5)
    axes[0].set_ylim(0, 1.5)

    axes[1].plot(rep1_data_sorted['mean_abundance'], rep1_data_sorted['CV'], 'o')
    axes[1].set_xlabel('Transcript abundance (log2(TPM+1))')
    axes[1].set_ylabel('CV')
    axes[1].set_title('CV vs Transcript Abundance')
    axes[1].set_ylim(0, 1.5)
    axes[1].set_xlim(0, 10)
    plt.tight_layout()
    plt.title(f"stats for {part1} and {part2}")
    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    # (AT)
    # plt.savefig(os.path.join(figure_dir, create_image_name(f"IM_CV_{type}" + timestamp + '.png')))
    plt.show()


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
    experiment_file = 'exprmnt_2024_07_16__12_46_32'
    main_dir = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/files/results/'
    directory = os.path.join(main_dir, experiment_file, 'files/output_files/')
    # paired_files, 2 types
    # type == 'replica': This will give you LR and SR file names with same replica, eg: lr_01_replica1+lr_01_replica2,
    # type == 'within_trainee': This will give you LR and SR file names those were trained together, eg: lr_01_replica1+sr_01_replica2,
    paired_files = pair_files(directory, type='replica')
    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    log_file = directory + f'corr_results{timestamp}.txt'
    sys.stdout = open(log_file, 'w')

    for pair in paired_files:
        print(pair)

        # Determine the type based on file names
        # if 'long' in pair[0] and 'long' in pair[1]:
        #     data_type = 'long'
        # elif 'short' in pair[0] and 'short' in pair[1]:
        #     data_type = 'short'
        # else:
        #     raise ValueError("File pair does not match expected 'long' or 'short' types")

        """ ## CALL 'spearman_corr_generic' FUNC FIRST ## """
        data_type = 'long'
        spearman_corr_generic(directory + pair[0], directory + pair[1], log_file, type=data_type)
        calculate_im_acvc(pair[0], pair[1], directory, log_file, type=data_type)

    sys.stdout.close()


# Example usage
if __name__ == "__main__":
    main()
