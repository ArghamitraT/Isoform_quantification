import pandas as pd
import sys

# Load the TSV file
file_path = sys.argv[1]
output_path = sys.argv[2]

df = pd.read_csv(file_path, sep="\t", header=None, names=['Isoform', 'Count'])

#Count the total number of reads
read_count = df['Count'].sum()

# Count the number of reads for each isoform
isoform_counts = df['Count']
tpms = isoform_counts/read_count * 10**6

# Convert the result to a DataFrame for better readability
isoform_counts_df = tpms.reset_index()
isoform_counts_df.columns = ['transcript_id', 'TPM']
isoform_counts_df['transcript_id'] = df['Isoform']


isoform_counts_df.to_csv(output_path)