import pandas as pd
from scipy.stats import spearmanr
import sys

# First Ill TSV file and then PB CSV
df = pd.read_csv(sys.argv[1], sep='\t') #Ill file
#df = pd.read_csv(sys.argv[1], sep = '\t', header=None, names=['transcript_id', 'count']) #PB file


# Merge the two DataFrames on the Isoform column
ill_isoforms = set(df['transcript_id'])
num_reads = df['count'].sum() #counts for input file



# Output the result
print('Number of Isoforms: ', len(ill_isoforms))
print('Number of Reads: ', num_reads)