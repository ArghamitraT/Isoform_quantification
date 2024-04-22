import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# FOR E0
# enter file name under file_path
file_path1 = "E0_tx_counts.tsv"
df1 = pd.read_csv(file_path1, sep="\t")
tpm_data1 = df1["tpm"]
# This is the E0 data
log_tpm1 = np.log(tpm_data1)

# FOR E2
# enter file name under file_path
file_path2 = "E2_tx_counts.tsv"
df2 = pd.read_csv(file_path2, sep="\t")
tpm_data2 = df2["tpm"]
log_tpm2 = np.log(tpm_data2)

tsv_df = pd.read_csv("E2_tx_counts.tsv", sep="\t")
tsv_df["transcript_name"] = tsv_df["transcript_name"].str[:-3]

# Read the .xlsx file into a DataFrame
xlsx_df = pd.read_csv("E2_molarity.tsv", sep="\t")

# Merge the DataFrames on the 'transcript_name' and 'Name' columns
merged_df = pd.merge(
    tsv_df, xlsx_df, left_on="transcript_name", right_on="Name", how="inner"
)


def getLogTpmE2(merged_df, val):
    filtered_df = merged_df[merged_df["E2"].astype(str) == val]
    print(filtered_df)
    tpm = filtered_df["tpm"]
    log_tpm = np.log(tpm)
    return log_tpm


# Filter out rows where the 'E2' column is not empty
log_tpm_2_132 = getLogTpmE2(merged_df, "  1/32")
log_tpm_2_14 = getLogTpmE2(merged_df, "  1/4 ")
log_tpm_2_1 = getLogTpmE2(merged_df, "1")
log_tpm_2_4 = getLogTpmE2(merged_df, "4")


# Replave E0 and E2 with appropraite labels. Eg: NanoCount, Liqa...
# Add more here in the format label: data
result = pd.DataFrame(
    {
        "E0": log_tpm1,
        "E2 1/32": log_tpm_2_132,
        "E2 1/4": log_tpm_2_14,
        "E2 1": log_tpm_2_1,
        "E2 4": log_tpm_2_4,
    }
)


sns.violinplot(data=result)
plt.title("Violin Plot of TPM Data")

# Can change ylim for better plots
plt.ylim(0, 15)
plt.ylabel("TPM (loge)")
plt.show()
