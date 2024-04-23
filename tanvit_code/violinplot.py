import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


# FOR E0
def getE0data(file_path):
    df = pd.read_csv(file_path, sep="\t")
    tpm_data = df["tpm"]
    log_tpm = np.log(tpm_data)
    return log_tpm


nc_log_tpm_E0 = getE0data("E0_tx_counts.tsv")
# Remove below comment for our model
# om_log_tpm_E0 = getE0data("E0_tx_counts.tsv from our model")
# DO THE SAME FOR OTHER SOTA


# FOR E2
def getLogTpmE2(merged_df, val):
    filtered_df = merged_df[merged_df["E2"].astype(str) == val]
    print(filtered_df)
    tpm = filtered_df["tpm"]
    log_tpm = np.log(tpm)
    return log_tpm


def getE2data(file_path_data, xlsx_df):
    tsv_df = pd.read_csv(file_path_data, sep="\t")
    tsv_df["transcript_name"] = tsv_df["transcript_name"].str[:-3]

    # Merge the DataFrames on the 'transcript_name' and 'Name' columns
    merged_df = pd.merge(
        tsv_df, xlsx_df, left_on="transcript_name", right_on="Name", how="inner"
    )

    log_tpm_2_132 = getLogTpmE2(merged_df, "  1/32")
    log_tpm_2_14 = getLogTpmE2(merged_df, "  1/4 ")
    log_tpm_2_1 = getLogTpmE2(merged_df, "1")
    log_tpm_2_4 = getLogTpmE2(merged_df, "4")

    return log_tpm_2_132, log_tpm_2_14, log_tpm_2_1, log_tpm_2_4


# Read the .xlsx file into a DataFrame
xlsx_df = pd.read_csv("E2_molarity.tsv", sep="\t")
nc_log_tpm_2_132, nc_log_tpm_2_14, nc_log_tpm_2_1, nc_log_tpm_2_4 = getE2data(
    "E2_tx_counts.tsv", xlsx_df
)
# Remove below comment for our model
# om_log_tpm_2_132,om_log_tpm_2_14,om_log_tpm_2_1,om_log_tpm_2_4 = getE2data("E2_tx_counts.tsv OUR MODEL",xlsx_df)
# DO THE SAME FOR OTHER SOTA


result_E0 = pd.DataFrame(
    {
        "NC E0": nc_log_tpm_E0,
    }
)
"""
#This is the real result after populating the data
result_E0 = pd.DataFrame(
    {
        "NC E0": nc_log_tpm_E0,
        "OM E0": om_log_tpm_E0
    }
)
"""

result_E2 = pd.DataFrame(
    {
        "NC E2 1/32": nc_log_tpm_2_132,
        "NC E2 1/4": nc_log_tpm_2_14,
        "NC E2 1": nc_log_tpm_2_1,
        "NC E2 4": nc_log_tpm_2_4,
    }
)

"""
#This is the real result after populating the data
result_E2 = pd.DataFrame(
    {
        "NC E2 1/32": nc_log_tpm_2_132,
        "OM E2 1/32": om_log_tpm_2_132,
        "NC E2 1/4": nc_log_tpm_2_14,
        "OM E2 1/4": nc_log_tpm_2_14,
        "NC E2 1": nc_log_tpm_2_1,
        "OM E2 1": nc_log_tpm_2_1,
        "NC E2 4": nc_log_tpm_2_4,
        "OM E2 4": nc_log_tpm_2_4,
    }
)
"""


sns.violinplot(data=result_E0)
plt.title("Violin Plot of TPM Data")

# Can change ylim for better plots
plt.ylim(0, 15)
plt.ylabel("TPM (loge)")
plt.show()

sns.violinplot(data=result_E2)
plt.title("Violin Plot of TPM Data")

# Can change ylim for better plots
plt.ylim(0, 15)
plt.ylabel("TPM (loge)")
plt.show()
