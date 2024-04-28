import sys
import pickle
import pandas as pd
import json
import os

from collections import defaultdict, OrderedDict

from json import JSONEncoder

#This allows python to serialize a custom class object
def _default(self, obj):
    return getattr(obj.__class__, "to_json", _default.default)(obj)

_default.default = JSONEncoder().default
JSONEncoder.default = _default

#If temp.txt exists, remove it
try:
    os.remove('temp.txt')
except OSError:
    pass

file = sys.argv[1]
# json_file = sys.argv[2]
pik_file = sys.argv[2]
len_file = sys.argv[3]

# file = 'PacBio/PacBio/day0_rep1/downsampled.txt'
# pik_file = 'PacBio/PacBio/day0_rep1/day0_rep1.pkl' 
# len_file = 'PacBio/PacBio/day0_rep1/day0_rep1_ref_len.pkl'

# pik_file = './PacBio/PacBio/day0_rep1/day0_rep1.pkl'

#Class to represent a read
class Read:
    def __init__(self):
        self.alignment_list = []
    def __repr__(self):
        m = ""
        for r in self.alignment_list:
            m += "\t\t{}\n".format(r)
        return m
    def add_alignment(self, alignment, **kwargs):
        self.alignment_list.append(alignment)
    def to_json(self):
        return self.__repr__()

#Class to represent an aligned segment
class Alignment:
    def __init__(self, q, r, qlen=0, align_len=0, align_score=0, secondary=False):
        """"""
        self.qname = q
        self.rname = r
        self.qlen = qlen
        self.align_len = align_len
        self.align_score = align_score
        self.secondary = secondary
        #self.read_idx = read_idx
        
    def __repr__(self):
        return "Query:{} | Reference:{} | Query len:{} | Alignment len:{} | Align Score:{} | Secondary:{}".format(
            self.qname,
            self.rname,
            self.qlen,
            self.align_len,
            self.align_score,
            self.secondary,
        )

count = 0
df = defaultdict(Read)
df_len = set()

with open(file) as f:
        header = f.readline()
        for new_line in f:
            line = new_line.strip('\n')
            q, r = line.split('\t')

            aln = Alignment(q, r)

            df[q].add_alignment(aln)
            df_len.add((r, 0))

            count+=1

            # if count%100000==0:
            #     print(count, ' lines have been processed')
            #     with open('temp.txt', 'a') as handle: 
            #         handle.write(json.dumps(df))
            #         handle.write('\n')
            #         df.clear()

with open(pik_file, 'wb') as handle:
    pickle.dump(df, handle)

df_len = OrderedDict(list(df_len))
with open(len_file, 'wb') as handle:
     pickle.dump(df_len, handle)

with open(pik_file, 'rb') as pkl:
    dp = pickle.load(pkl)

with open(len_file, 'rb') as pkl:
    dp_len = pickle.load(pkl)

# #Writing all the left-over contents in df
# with open('temp.txt', 'a') as handle: 
#     handle.write(json.dumps(df))
#     handle.write('\n')
#     df.clear()

# #Parse temp.txt to remove extra brackets and add appropriate commas
# with open('temp.txt', 'r') as temp_file:
#     with open(json_file, 'w') as out:
#         first_line = temp_file.readline()
#         out.write(first_line.strip('\n').strip('}')+',')

#         prev_line = None
#         for line in temp_file:
#             if prev_line is not None:
#                 out.write(prev_line.strip('\n').strip('}{')+',')
#             prev_line = line
#         out.write(prev_line.strip('\n').strip('{'))            

print()
print()
print('Total lines in TXT file: ',count)
print()
print('Parsing Complete')