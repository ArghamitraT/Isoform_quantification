import pandas as pd
from scipy.stats import spearmanr
import sys

# First Ill TSV file and then PB CSV
df1 = pd.read_csv(sys.argv[1], sep='\t') #Ill file
df2 = pd.read_csv(sys.argv[2]) #PB file

# Merge the two DataFrames on the Isoform column
ill_isoforms = set(df1['transcript_id'])
pb_isoforms = set(df2['transcript_id'])

# Calculate Spearman correlation between the Read_Count columns
common_isoforms = ill_isoforms.intersection(pb_isoforms)
total_isoforms = ill_isoforms.union(pb_isoforms)

# Output the result
print('Number of Common Isoforms: ', len(common_isoforms))
print('Total Number of Isoforms: ', len(total_isoforms))