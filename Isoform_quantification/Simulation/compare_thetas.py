import pandas as pd
from scipy.stats import spearmanr
import sys

# Example DataFrames (replace these with your actual data)
df1 = pd.read_csv(sys.argv[1], sep='\t')
df2 = pd.read_csv(sys.argv[2], sep='\t')

# Merge the two DataFrames on the Isoform column
merged_df = pd.merge(df1, df2, on='#transcript_id', suffixes=('_df1', '_df2'))

# Calculate Spearman correlation between the Read_Count columns
correlation, p_value = spearmanr(merged_df['tpm_df1'], merged_df['tpm_df2'])

# Output the result
print(f"Spearman Correlation: {correlation}")
print(f"P-value: {p_value}")