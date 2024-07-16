""" This script generates violin plot comparing with SOTA """

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
import generate_corr as gen_img

"""
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

# Function to read data and return log-transformed TPM
def getCV(file_path):
    df = pd.read_csv(file_path, sep="\t")
    tpm_data = df["tpm"]
    # Calculate the Coefficient of Variation (CV)
    cv = tpm_data.std() / tpm_data.mean()
    print(f"CV for {file_path.split('/')[-1]}: {cv}")

def get_spearsman_corr(file_path):
    our_quant = pd.read_csv(file_path, sep="\t")
    #our_quant = df["tpm"]
    ground_truth = pd.read_csv('../../data/SIRV/E2_molarity.csv')
    # Reset the index of our_quant so that the isoform names become a column
    our_quant_reset = our_quant.reset_index()

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
    print(f"Spearman correlation coefficient {file_path.split('/')[-1]}: {correlation}")
    print(f"P-value {file_path.split('/')[-1]}: {p_value}")

# Function to read data and return log-transformed TPM
def getLogTPM(file_path):
    df = pd.read_csv(file_path, sep="\t")
    tpm_data = df["tpm"]
    log_tpm = np.log(tpm_data + 1)  # Adding 1 to avoid log(0)
    return log_tpm

# Function to process specific E2 condition data from a merged DataFrame
def getLogTPM_E2(merged_df, condition):
    filtered_df = merged_df[merged_df["E2"].astype(str).str.strip() == condition]
    tpm = filtered_df["tpm"]
    log_tpm = np.log(tpm + 1)  # Adding 1 to avoid log(0)
    return log_tpm

# Function to process E2 data and handle multiple conditions
def getE2data(file_path_data, xlsx_df):
    tsv_df = pd.read_csv(file_path_data, sep="\t")
    tsv_df["transcript_name"] = tsv_df["transcript_name"].str[:-3]
    merged_df = pd.merge(tsv_df, xlsx_df, left_on="transcript_name", right_on="Name", how="inner")

    conditions = ["1/32", "1/4", "1", "4"]
    log_tpms = {f"E2 {c}": getLogTPM_E2(merged_df, c) for c in conditions}
    return log_tpms

# Specify the base path for data files
main_dir = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/CU_courses/Spring_24/CBMF4761/Project/RNA_Splicing/data/SIRV/Anvit_SOTA_Output/'
xlsx_df = pd.read_csv(main_dir + "E2_molarity.tsv", sep="\t")

# Dictionary of E0 data paths
data_paths_E0 = {
    "StringTie2 E0": main_dir + "nanoCount_E0_output.tsv",
    "NanoCount E0": main_dir + "nanoCount_E0_output.tsv",
    "OurModel E0": main_dir + "SIRV_output_sample1_aln_E0.tsv",
    # Additional data sources can be added here
}

# Load and transform E0 data
{label: getCV(path) for label, path in data_paths_E0.items()}

log_tpm_data_E0 = {label: getLogTPM(path) for label, path in data_paths_E0.items()}
combined_E0 = pd.DataFrame(log_tpm_data_E0)

plt.figure(figsize=(12, 6))
sns.violinplot(data=combined_E0)
plt.title("Violin Plot of E0 TPM Data")
plt.ylim(0, 15)
plt.ylabel("TPM (loge)")
plt.xticks(fontsize=25)  # Increase x-axis tick label font size
plt.yticks(fontsize=25)
plt.grid()
#plt.savefig('figures/'+gen_img.create_image_name('SIRV_E0'))
plt.show()

# Dictionary of E2 data paths
data_paths_E2 = {
    "StringTie2 E2": main_dir + "nanoCount_E2_output.tsv",
    "NanoCount E2": main_dir + "nanoCount_E2_output.tsv",
    "OurModel E2": main_dir + "SIRV_output_sample2_aln_E2.tsv",
    # Additional data sources can be added here
}

# {label: get_spearsman_corr(path) for label, path in data_paths_E2.items()}

# Load and transform E2 data for all conditions
combined_E2_data = {}
conditions = ["1/32", "1/4", "1", "4"]

# Initialize a dictionary to hold the data for plotting
plot_data = {f"{condition} {dataset}": [] for condition in conditions for dataset in ["StringTie2 E2", "NanoCount E2", "OurModel E2"]}

for dataset_label, path in data_paths_E2.items():
    e2_data = getE2data(path, xlsx_df)
    for condition in conditions:
        plot_data[f"{condition} {dataset_label}"] = e2_data[f"E2 {condition}"].values


# Initialize an empty list for the rows of the new DataFrame
data_rows = []

# Iterate over each key-value pair in the plot_data dictionary
for group_label, values in plot_data.items():
    # Extract the condition and dataset from the group_label
    dataset = group_label.split(' ')[1]
    condition = group_label.split(' ')[0]
    # Create a dictionary for each value and append to the list
    for value in values:
        data_rows.append({'Condition': condition, 'Dataset': dataset, 'Value': value})

# Now, create a DataFrame from the list of dictionaries
df = pd.DataFrame(data_rows)

# Plotting
plt.figure(figsize=(16, 6))
violin = sns.violinplot(x='Condition', y='Value', hue='Dataset', data=df, split=False)
plt.title("Side-by-side Violin Plot of E2 TPM Data by Condition and Dataset")
plt.ylim(0, 15)
plt.ylabel("TPM (loge)")
#plt.legend(title='Dataset')
violin.legend_.remove()
plt.xticks(fontsize=25)  # Increase x-axis tick label font size
plt.yticks(fontsize=25)
plt.grid()
plt.savefig('figures/'+gen_img.create_image_name('SIRV_E2'))
plt.show()


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Example DataFrame
data = {
    'Category': ['StringTie2', 'NanoCount', 'Our Model'],
    'Values': [0.46, 0.41, 0.31]
}
df = pd.DataFrame(data)

# Create a bar plot
plt.figure(figsize=(5, 3))  # Set the figure size
sns.barplot(x='Category', y='Values', data=df, width=0.5)

# Adding title and labels
plt.title('SIRV E0 Coefficient of Variance')
# plt.xlabel('Category')
plt.ylabel('CV')
plt.grid()
plt.savefig('figures/'+gen_img.create_image_name('SIRV_E0_CV'))
# Show the plot
plt.show()


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Example DataFrame
data = {
    'Category': ['StringTie2', 'NanoCount', 'Our Model'],
    'Values': [0.79, 0.87, 0.87]
}
df = pd.DataFrame(data)

# Create a bar plot
plt.figure(figsize=(5, 3))  # Set the figure size
sns.barplot(x='Category', y='Values', data=df, width=0.5)

# Adding title and labels
plt.title('SIRV E2 Spearman Correlation')
# plt.xlabel('Category')
plt.ylabel('Correlation')
plt.grid()
plt.savefig('figures/'+gen_img.create_image_name('SIRV_E2_Corr'))
# Show the plot
plt.show()
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Example DataFrame
# data = {
#     'Category': ['StringTie2', 'NanoCount', 'Our Model'],
#     'Values': [0.16, 0.11, 0.25]
# }
# df = pd.DataFrame(data)
#
# # Create a bar plot
# plt.figure(figsize=(6, 3))  # Set the figure size
#
# palette = {'NanoCount': 'orange', 'Our Model': 'green'}
# #sns.barplot(x='Category', y='Values', data=df, palette=palette, width=0.25)
# sns.barplot(x='Category', y='Values', data=df, width=0.5)
#
#
# # Adding title and labels
# plt.title('Replicate 2 Spearman Correlation')
# # plt.xlabel('Category')
# plt.ylabel('Correlation')
# plt.grid()
# plt.savefig('figures/'+gen_img.create_image_name('rep2_Corr'))
# # Show the plot
# plt.show()

data = {
    'Category': ['StringTie2', 'NanoCount', 'Our Model'],
    'Values': [0.14, 0.075, 0.22]
}
df = pd.DataFrame(data)

# Create a bar plot
plt.figure(figsize=(6, 3))  # Set the figure size

palette = {'NanoCount': 'orange', 'Our Model': 'green'}
#sns.barplot(x='Category', y='Values', data=df, palette=palette, width=0.25)
sns.barplot(x='Category', y='Values', data=df, width=0.5)


# Adding title and labels
plt.title('Replicate 1 Spearman Correlation')
# plt.xlabel('Category')
plt.ylabel('Correlation')
plt.grid()
plt.savefig('figures/'+gen_img.create_image_name('rep1_Corr'))
# Show the plot
plt.show()
# # Load E2 data for all conditions and datasets
# e2_data_list = []
# for label, path in data_paths_E2.items():
#     e2_data = getE2data(path, xlsx_df)
#     for condition, data in e2_data.items():
#         temp_df = pd.DataFrame({
#             'Value': data,
#             'Condition': condition,
#             'Dataset': label
#         })
#         e2_data_list.append(temp_df)
#
# # Combine all E2 data into a single DataFrame
# combined_e2_df = pd.concat(e2_data_list)
#
# # Create a new column in the DataFrame for plotting purposes
# # This column combines the condition and dataset name for a grouped violin plot
# combined_e2_df['Plot Group'] = combined_e2_df['Dataset'] + ' ' + combined_e2_df['Condition']
#
# # Now create the violin plot with seaborn
# plt.figure(figsize=(12, 6))
# sns.violinplot(x='Plot Group', y='Value', data=combined_e2_df, palette="muted")
# plt.xticks(rotation=90)  # Rotate x labels for better readability if needed
# plt.title("Side-by-side Violin Plot of E2 TPM Data by Condition and Dataset")
# plt.ylim(0, 15)
# plt.ylabel("TPM (loge)")
# plt.tight_layout()  # Adjust layout to prevent clipping of tick-labels
# plt.savefig('figures/SIRV_E2_grouped.png')
# plt.show()











#
# # Load and transform E2 data for all conditions
# combined_E2_data = {}
# for label, path in data_paths_E2.items():
#     combined_E2_data.update({
#         f"{label} E2 {cond}": data for cond, data in getE2data(path, xlsx_df).items()
#     })
#
# # Create DataFrames for E0 and E2 data
# combined_E0 = pd.DataFrame(log_tpm_data_E0)
# combined_E2 = pd.DataFrame(combined_E2_data)
#
# # Plotting
# plt.figure(figsize=(12, 6))
# sns.violinplot(data=combined_E0)
# plt.title("Violin Plot of E0 TPM Data")
# plt.ylim(0, 15)
# plt.ylabel("TPM (loge)")
# plt.savefig('figures/SIRV_E0.png')
# plt.show()
#
# plt.figure(figsize=(12, 6))
# sns.violinplot(data=combined_E2)
# plt.title("Violin Plot of E2 TPM Data")
# plt.ylim(0, 15)
# plt.ylabel("TPM (loge)")
# plt.savefig('figures/SIRV_E2.png')
# plt.show()






#
# import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
# import numpy as np
#
#
# # FOR E0
# def getE0data(file_path):
#     df = pd.read_csv(file_path, sep="\t")
#     tpm_data = df["tpm"]
#     log_tpm = np.log(tpm_data)
#     return log_tpm
#
# main_dir = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/CU_courses/Spring_24/CBMF4761/Project/RNA_Splicing/data/SIRV/Anvit_NanoCount_Out_E0_E2/'
# nc_log_tpm_E0 = getE0data(main_dir+"E0_tx_counts.tsv")
# om_log_tpm_E0 = getE0data(main_dir+"SIRV_output_sample1_aln_E0.tsv")
#
#
# # FOR E2
# def getLogTpmE2(merged_df, val):
#     filtered_df = merged_df[merged_df["E2"].astype(str) == val]
#     print(filtered_df)
#     tpm = filtered_df["tpm"]
#     log_tpm = np.log(tpm)
#     return log_tpm
#
#
# def getE2data(file_path_data, xlsx_df):
#     tsv_df = pd.read_csv(file_path_data, sep="\t")
#     tsv_df["transcript_name"] = tsv_df["transcript_name"].str[:-3]
#
#     # Merge the DataFrames on the 'transcript_name' and 'Name' columns
#     merged_df = pd.merge(
#         tsv_df, xlsx_df, left_on="transcript_name", right_on="Name", how="inner"
#     )
#
#     log_tpm_2_132 = getLogTpmE2(merged_df, "  1/32")
#     log_tpm_2_14 = getLogTpmE2(merged_df, "  1/4 ")
#     log_tpm_2_1 = getLogTpmE2(merged_df, "1")
#     log_tpm_2_4 = getLogTpmE2(merged_df, "4")
#
#     return log_tpm_2_132, log_tpm_2_14, log_tpm_2_1, log_tpm_2_4
#
#
# # Read the .xlsx file into a DataFrame
# xlsx_df = pd.read_csv(main_dir+"E2_molarity.tsv", sep="\t")
#
# nc_log_tpm_2_132, nc_log_tpm_2_14, nc_log_tpm_2_1, nc_log_tpm_2_4 = getE2data(main_dir+"E2_tx_counts.tsv", xlsx_df)
# om_log_tpm_2_132,om_log_tpm_2_14,om_log_tpm_2_1,om_log_tpm_2_4 = getE2data(main_dir+"SIRV_output_sample2_aln_E2.tsv",xlsx_df)
#
# # Remove below comment for our model
# # DO THE SAME FOR OTHER SOTA
#
#
# result_E0 = pd.DataFrame(
#     {
#         "NC E0": nc_log_tpm_E0,
#     }
# )
#
# result_E2 = pd.DataFrame(
#     {
#         "NC E2 1/32": nc_log_tpm_2_132,
#         "NC E2 1/4": nc_log_tpm_2_14,
#         "NC E2 1": nc_log_tpm_2_1,
#         "NC E2 4": nc_log_tpm_2_4,
#     }
# )
#
#
# sns.violinplot(data=result_E0)
# plt.title("Violin Plot of TPM Data")
#
# # Can change ylim for better plots
# plt.ylim(0, 15)
# plt.ylabel("TPM (loge)")
# plt.savefig('figures/SIRV_E0.png')
# plt.show()
#
# sns.violinplot(data=result_E2)
# plt.title("Violin Plot of TPM Data")
#
# # Can change ylim for better plots
# plt.ylim(0, 15)
# plt.ylabel("TPM (loge)")
# plt.savefig('figures/SIRV_E2.png')
# plt.show()