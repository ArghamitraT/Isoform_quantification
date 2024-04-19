import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
import pandas as pd


def create_image_name(name, format=".png"):
    import datetime
    import time
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