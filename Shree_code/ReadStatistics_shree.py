import sys
import math
import numpy as np

from Bio import SeqIO

file = sys.argv[1]
file_type = sys.argv[2]


def fasta_len(fasta_file):

    #record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fastq"))

    lengths = {}
    records = SeqIO.parse(fasta_file, 'fastq')
    for record in records:   #loop through each fasta entry
        length = len(record.seq)    #get sequence length
        lengths[record.id] = length

    #access dictionary outside of loop
    return np.mean(np.array([x[1] for x in list(lengths.items())]))

def txt_len(txt_file):
    lengths = []
    with open(txt_file) as f:
        header = f.readline()
        for new_line in f:
            line = new_line.strip('\n')
            q, r = line.split('\t')
            lengths.append(len(q))
    return np.mean(lengths)


if file_type == 'fasta':
    avg = fasta_len(file)
    print(avg)

elif file_type == 'txt':
    avg = txt_len(file)
    print(avg)