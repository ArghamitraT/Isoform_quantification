import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
import pandas as pd
import datetime
import time
import os

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


def spearman_corr(our_quant, sample):

    ground_truth = pd.read_csv('../../data/SIRV/E2_molarity.csv')
    # Reset the index of our_quant so that the isoform names become a column
    our_quant_reset = our_quant.reset_index()

    # Rename the columns accordingly
    our_quant_reset.columns = ['transcript_name', 'other_column1', 'est_count', 'tpm']

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

def spearman_corr_generic(file_path1, file_path2):
    # Load the datasets
    our_quant = pd.read_csv(file_path1, sep="\t")
    ground_truth = pd.read_csv(file_path2, sep="\t")

    # Apply the function to the 'transcript_name' column of your DataFrame
    our_quant['transcript_name'] = our_quant['transcript_name'].apply(split_transcript_name)
    ground_truth['transcript_name'] = ground_truth['transcript_name'].apply(split_transcript_name)

    # Reset the index in our_quant so that the isoform names become a column
    # our_quant_reset = our_quant.reset_index()

    # Clean the isoform names in both datasets to match the naming convention
    # our_quant_reset['cleaned_name'] = our_quant_reset['transcript_name'].str.replace(r'\(\+\)|\(\-\)', '', regex=True)
    # ground_truth['cleaned_name'] = ground_truth['Name'].str.replace(r'\(\+\)|\(\-\)', '', regex=True)

    # Find common isoforms
    common_isoforms = pd.merge(our_quant, ground_truth, on='transcript_name', suffixes=('_quant', '_truth'))

    # Ensure both datasets have TPM values in the respective columns named accurately
    if 'tpm_quant' in common_isoforms and 'tpm_truth' in common_isoforms:
        # Convert any fractional notation in TPM values from ground_truth to float if necessary
        common_isoforms['tpm_truth'] = common_isoforms['tpm_truth'].apply(fraction_to_float_gen)

        # Calculate Spearman's correlation using the matched TPM lists
        correlation, p_value = spearmanr(common_isoforms['tpm_quant'], common_isoforms['tpm_truth'])

        # Output the results
        print(
            f"Spearman correlation coefficient between {file_path1.split('/')[-1]} & {file_path2.split('/')[-1]}: {correlation}")
        print(f"P-value for correlation: {p_value}")
    else:
        print("TPM columns missing or incorrectly named in one of the datasets.")


output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/')

# REP 2, big
# file_path1 = output_file+'output_day0rep2_sample2_Day0_Pacbio_2024_4_27_15_14_44.tsv'
# file_path2 = output_file+'output_day0rep2_sample1_Day0_illumina_2024_4_27_15_14_13.tsv'

# REP 1, big
# file_path1 = output_file+'output_day0rep1_sample1_Day0_illumina_2024_4_27_14_54_52.tsv'
# file_path2 = output_file+'output_day0rep1_sample2_Day0_Pacbio_2024_4_27_14_55_49.tsv'

# REP 2, small
# file_path1 = output_file+'output_day0rep2_sample2_Day0_2_Pacbio_2024_4_26_20_47_51.tsv'
# file_path2 = output_file+'output_day0rep2_sample1_Day0_1_illumina_2024_4_26_20_47_11.tsv'

# REP 1, small
# file_path1 = output_file+'output_day0rep1_sample2_Day0_2_Pacbio_2024_4_26_20_11_23.tsv'
# file_path2 = output_file+'output_day0rep1_sample1_Day0_1_illumina_2024_4_26_20_11_17.tsv'

####### cross replicate #########
# REP 1 small
# file_path1 = output_file+'output_IlluRep1_PacRep2_sample1_Day0_1_illumina_2024_4_27_20_48_19.tsv'
# file_path2 = output_file+'output_IlluRep2_PacRep1_sample2_Day0_2_Pacbio_2024_4_27_20_43_04.tsv'

# REP 2 small
# file_path1 = output_file+'output_IlluRep2_PacRep1_sample1_Day0_1_illumina_2024_4_27_20_43_03.tsv'
# file_path2 = output_file+'output_IlluRep1_PacRep2_sample2_Day0_2_Pacbio_2024_4_27_20_48_19.tsv'


# REP 1 big
# file_path1 = output_file+'output_IlluRep1_PacRep2_BIG_sample1_Day0_1_illumina_2024_4_27_21_04_37.tsv'
# file_path2 = output_file+'output_IlluRep2_PacRep1_BIG_sample2_Day0_2_Pacbio_2024_4_27_21_18_20.tsv'

# REP 2 small
# file_path1 = output_file+'output_IlluRep2_PacRep1_BIG_sample1_Day0_1_illumina_2024_4_27_21_18_19.tsv'
# file_path2 = output_file+'output_IlluRep1_PacRep2_BIG_sample2_Day0_2_Pacbio_2024_4_27_21_04_38.tsv'

""" NanoCount """
# REP 2, big
# file_path1 = output_file+'Illu_day0rep2_NC.tsv'
# file_path2 = output_file+'Pac_bigDown_day0rep2_NC.tsv'

# # REP 1, big
# file_path1 = output_file+'Illu_day0rep1_NC.tsv'
# file_path2 = output_file+'Pac_bigDown_day0rep1_NC.tsv'

# REP 1, small
# file_path1 = output_file+'Illu_day0rep1_NC.tsv'
# file_path2 = output_file+'Pac_smallDown_day0rep1_NC.tsv'

# REP 2, small
# file_path1 = output_file+'Illu_day0rep2_NC.tsv'
# file_path2 = output_file+'Pac_smallDown_day0rep2_NC.tsv'

""" StringTie2 """
# REP 1, small
file_path1 = ('/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/CU_courses/Spring_24/CBMF4761/Project/RNA_Splicing/data/StringTie_Illumina_Output/'
              'stringTie_Day0_1_illumina_output.tsv')
file_path2 = output_file+'Pac_smallDown_day0rep1_NC.tsv'

# REP 2, small
file_path1 = ('/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/CU_courses/Spring_24/CBMF4761/Project/RNA_Splicing/data/StringTie_Illumina_Output/'
              'stringTie_Day0_2_illumina_output.tsv')
file_path2 = output_file+'Pac_smallDown_day0rep2_NC.tsv'

spearman_corr_generic(file_path1, file_path2)