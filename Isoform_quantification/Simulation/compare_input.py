# To compare output abundance file with input abundance file

import pandas as pd
from scipy.stats import spearmanr
import sys

df1 = pd.read_csv(sys.argv[1], sep='\t') #output file
df2 = pd.read_csv(sys.argv[2], sep='\t') #input file

# Merge the two DataFrames on the Isoform column
merged_df = pd.merge(df1, df2, left_on='transcript_id', right_on='#transcript_id', suffixes=('_df1', '_df2'))

print(merged_df.head())

# Calculate Spearman correlation between the Read_Count columns
correlation, p_value = spearmanr(merged_df['TPM'], merged_df['tpm'])

# Output the result
print(f"Spearman Correlation: {correlation}")
print(f"P-value: {p_value}")