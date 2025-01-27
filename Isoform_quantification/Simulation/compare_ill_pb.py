import pandas as pd
from scipy.stats import spearmanr
import sys

# First Ill TSV file and then PB CSV
df1 = pd.read_csv(sys.argv[1], sep='\t') #Ill file
df2 = pd.read_csv(sys.argv[2]) #PB file

# Merge the two DataFrames on the Isoform column
merged_df = pd.merge(df1, df2, on='transcript_id', suffixes=('_df1', '_df2'))

# Calculate Spearman correlation between the Read_Count columns
correlation, p_value = spearmanr(merged_df['TPM_df1'], merged_df['TPM_df2'])

# Output the result
print(f"Spearman Correlation: {correlation}")
print(f"P-value: {p_value}")