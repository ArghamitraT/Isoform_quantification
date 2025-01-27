import pandas as pd
import sys

# Load the TSV file
file_path = sys.argv[1]
output_path = sys.argv[2]

df = pd.read_csv(file_path, sep="\t", header=None, names=['Reads', 'Isoform'])

# Count the number of reads for each isoform
isoform_counts = df['Isoform'].value_counts()

# Convert the result to a DataFrame for better readability
isoform_counts_df = isoform_counts.reset_index()
isoform_counts_df.columns = ['Isoform', 'Read_Count']

isoform_counts_df.to_csv(output_path)