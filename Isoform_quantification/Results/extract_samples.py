import numpy as np
import pandas as pd
import sys
import os

isoseq = pd.read_csv(sys.argv[1], sep='\t', header=None, names=['transcript_id', 'count']) #Isoseq ground truth
rsem = pd.read_csv(sys.argv[2], sep='\t') #RSEM ground truth

expt = sys.argv[3]

df = pd.read_csv(sys.argv[4], sep='\t')

short = df[df['transcript_name'].isin(rsem['transcript_id'])]

long = df[df['transcript_name'].isin(isoseq['transcript_id'])]

# print(set(rsem['transcript_id']) >= set(df['transcript_name']))

# print(set(isoseq['transcript_id']) >= set(df['transcript_name']))

#print(short.head())
# print(df.shape)
# print(long.shape)

short.to_csv(os.path.join(expt, 'short.tsv'), sep='\t')
long.to_csv(os.path.join(expt, 'long.tsv'), sep='\t')