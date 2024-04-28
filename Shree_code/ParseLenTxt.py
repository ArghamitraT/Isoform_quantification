import sys
import pickle
import pandas as pd
import json
import os

from collections import defaultdict, OrderedDict

file = sys.argv[1]
pik_file = sys.argv[2]

# file =  './PacBio/PacBio/day0_rep1/collapsed_transcripts.read_stat.txt'
# pik_file = './PacBio/PacBio/day0_rep1/day0_rep_1_ref_len.pkl'


count = 0
df = set()
with open(file) as f:
        header = f.readline()
        for new_line in f:
            line = new_line.strip('\n')
            q, r = line.split('\t')

            df.add((r, 0))
            count+=1

df = OrderedDict(list(df))

with open(pik_file, 'wb') as handle:
     pickle.dump(df, handle)

print('Parsing Complete')