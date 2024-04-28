import sys
import pickle
import pandas as pd
import json
import os
import numpy as np
import random

import matplotlib.pyplot as plt

from collections import defaultdict, Counter

from json import JSONEncoder

MAX_ISOFORM_COUNT = 200000
MAX_READ_COUNT = 500000
MAX_READS_PER_ISOFORM = 1000

aln_path = 'PacBio/PacBio/day0_rep2/collapsed_transcripts.read_stat.txt'
file_path = 'PacBio/PacBio/day0_rep2/collapsed_transcripts.flnc_count.txt'
text_path = 'PacBio/PacBio/day0_rep2/downsampled_new.txt'


def plot_normalized_read_distribution(sample, title, nameFig):
    # Normalize the values
    values = np.array(list(sample.values()), dtype=float)
    min_value = np.min(values)
    max_value = np.max(values)
    normalized_values = (values - min_value) / (max_value - min_value)

    # Define 10 bins for normalized data which ranges from 0 to 1
    bins = np.linspace(0, 1, 11)

    # Calculate the count of isoforms in each bin
    hist, bin_edges = np.histogram(normalized_values, bins=bins)

    # Normalize the histogram counts for the y-axis
    normalized_hist = hist / float(hist.sum())

    # Plot the histogram with the bin counts
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(range(len(normalized_hist)), normalized_hist, width=0.8, align='center')

    # Annotate each bar with the actual range of data values
    for bar, bin_count, bin_start, bin_end in zip(bars, hist, bin_edges[:-1], bin_edges[1:]):
        actual_start = int(bin_start * (max_value - min_value) + min_value)
        actual_end = int(bin_end * (max_value - min_value) + min_value)
        annotation_text = f'Reads #: {actual_start}-{actual_end}\nIsoforms #: {bin_count}'

        # Dynamically set the vertical position based on the height of the bar
        vertical_position = bar.get_height() if bar.get_height() > 0 else 0.1
        ax.annotate(annotation_text,
                    xy=(bar.get_x() + bar.get_width() / 2, vertical_position),
                    xytext=(0, 5),  # Adjust text position for visibility
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8, rotation=45, color='blue')

    # Set the title and labels
    plt.yscale('log')  # Log scale for y-axis
    plt.title(title)
    plt.xlabel('Normalized Value Bins')
    plt.ylabel('Isoforms # (normalized)')

    # Adjust the x-tick labels to match the bin edges and prevent overlap
    ax.set_xticks(range(len(bin_edges) - 1))
    ax.set_xticklabels([f"{int(edge):d}" for edge in bin_edges[:-1]], rotation=45)

    # Save the figure to a file and show it
    plt.tight_layout()
    plt.savefig(nameFig)
    plt.show()

def convert_text_to_counter(file_path):
    # Create an empty Counter object
    isoform_counter = Counter()

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            # Skip the header line
            if i == 0:
                continue

            # Extract the relevant data
            parts = line.strip().split()
            # Assume that the first column is the isoform identifier and the second column is the count
            isoform_id = parts[0].strip().split(',')[0]
            count_fl = int(parts[0].strip().split(',')[1])  # Convert count to integer

            # isoform_id = parts[0]
            # count_fl = int(parts[1])  # Convert count to integer

            # Update the Counter object
            isoform_counter[isoform_id] = count_fl

    return isoform_counter

def get_isoforms(sample):

    kept_isoforms = []

    count_dict = defaultdict(int)

    values = np.array(list(sample.values()), dtype=float)
    min_value = np.min(values)
    max_value = np.max(values)
    normalized_values = (values - min_value) / (max_value - min_value)


    bins = np.linspace(0, 1, 11)

    # Calculate the count of isoforms in each bin
    hist, bin_edges = np.histogram(normalized_values, bins=bins)

    count=0

    bin_num = -1

    isoforms_per_bin = {}

    for bin_count, bin_start, bin_end in zip(hist, bin_edges[:-1], bin_edges[1:]):
        bin_num+=1
        high_num_reads = 0
        actual_start = int(bin_start * (max_value - min_value) + min_value)
        actual_end = int(bin_end * (max_value - min_value) + min_value)

        count_dict[count] = (hist[count], actual_start, actual_end)
        count+=1

        isoforms_in_bin = []
        for key in sample.keys():
            if sample[key] <= actual_end and sample[key]>=actual_start:
                isoforms_in_bin.append(key)
        
        
        #Keep isoforms that are in bins with less than MAX_ISOFORM_COUNT isoforms
        if bin_count<MAX_ISOFORM_COUNT:
            kept_isoforms.append(isoforms_in_bin)
            high_num_reads = 1
            isoforms_per_bin[bin_num] = (set(isoforms_in_bin), high_num_reads)
            continue

        
        to_keep = []

        less_than_5_per = []

        #Idea: Of total isoforms kept - 65% should be the ones with the least 5% of reads, and the other 35% are sampled randomly

        #Isoforms with low number of reads (least 5% of reads)
        threshold_num_reads = int(0.01 * actual_end)
        for i in isoforms_in_bin:
            if sample[i]<=threshold_num_reads:
                less_than_5_per.append(i)
        
        if len(less_than_5_per) < 0.7*MAX_ISOFORM_COUNT:
            diff = MAX_ISOFORM_COUNT - len(less_than_5_per)
            random_sample = np.random.choice(isoforms_in_bin, diff)
            to_keep.extend(random_sample)
        
        else:
            t = int(0.7 * MAX_ISOFORM_COUNT)
            o = int(0.3 * MAX_ISOFORM_COUNT)
            r1 = np.random.choice(less_than_5_per, t)
            r2 = np.random.choice(isoforms_in_bin, o)
            to_keep.extend(r1)
            to_keep.extend(r2)
        
        isoforms_per_bin[bin_num] = (set(to_keep), high_num_reads)

        
        kept_isoforms.append(to_keep)
        

    return kept_isoforms, isoforms_per_bin, threshold_num_reads



counter_from_text_file = convert_text_to_counter(file_path)
res, isoforms_per_bin, threshold_num_reads = get_isoforms(counter_from_text_file)

binned = {}
for i in range(len(res)):
    for iso in res[i]:
        binned[iso] = i+1

all_kept = []
for i in res:
    all_kept.extend(i)

all_kept = set(all_kept)

filtered_counts = {}

for i in all_kept:
    filtered_counts[i] = counter_from_text_file[i]


# num_reads = sum(filtered_counts.values())
# per_reads = MAX_READ_COUNT/num_reads

random.seed(42)

# selection_indices = random.sample(range(num_reads), MAX_READ_COUNT)
# selection_indices.sort(reverse=True)

df = pd.DataFrame(columns = ['isoforms', 'counts', 'bins'])
df['isoforms'] = filtered_counts.keys()
df['counts'] = filtered_counts.values()

bins = []
for i in filtered_counts.keys():
    bins.append(binned[i])

df['bins'] = bins

df_downsample = df.sample(frac=0.55, weights='counts', axis=0)
d = {k:v for (k,v) in zip(df_downsample['isoforms'].astype('str'), df_downsample['counts'].astype('int'))}
# plot_normalized_read_distribution(d, 'Downsampled Read Distribution', 'long_ds_1.png')

# print(sum(d.values()))

df1 = df_downsample[df_downsample['bins']!=1]
df2 = df_downsample[df_downsample['bins']==1]
df_downsample_new = df1.sample(frac=0.19, weights='counts', axis=0)
d1 = {k:v for (k,v) in zip(df_downsample_new['isoforms'].astype('str'), df_downsample_new['counts'].astype('int'))}
d2 = {k:v for (k,v) in zip(df2['isoforms'].astype('str'), df2['counts'].astype('int'))}

d1.update(d2)
# print(sum(d1.values()))

# plot_normalized_read_distribution(d1, 'Downsampled Read Distribution', 'long_ds_2.png')

df1 = df_downsample_new[df_downsample_new['bins']>6]
df2 = df_downsample_new[df_downsample_new['bins']<=6]
df_downsample_new = df_downsample_new.sample(frac=0.009, weights='counts', axis=0)
d3 = {k:v for (k,v) in zip(df_downsample_new['isoforms'].astype('str'), df_downsample_new['counts'].astype('int'))}
d4 = {k:v for (k,v) in zip(df2['isoforms'].astype('str'), df2['counts'].astype('int'))}

# for k in d.keys():
#     if k not in d1.keys():
#         d1[k] = d[k]

d3.update(d4)
d3.update(d2)
print(sum(d3.values()))

#plot_normalized_read_distribution(d3, 'Downsampled Read Distribution', 'long_ds_3.png')


final_counts = defaultdict(int)

final_isoforms = set(d3.keys())

count = 0
low_read_num_threshold = int(0.1*MAX_READ_COUNT)
count_low_reads = 0
with open(text_path, 'w') as text_file:
    with open(aln_path) as f:
            header = f.readline()
            for new_line in f:
                line = new_line.strip('\n')
                q, r = line.split('\t')
                if r in final_isoforms: #all_kept:
                    # if r in isoforms_per_bin[0][0]:
                    #     for i in range(filtered_counts[r]):
                    #         if i >= MAX_READS_PER_ISOFORM:
                    #             break
                    #         if selection_indices==[]:
                    #             break
                    #         selection_indices.pop()
                    #         final_counts[r]+=1
                    #         text_file.write(q+'\t'+r+'\n')
                    # else:
                    #     if selection_indices==[]:
                    #         break
                    #     if count==selection_indices[-1]:
                    #         selection_indices.pop()
                    #         final_counts[r]+=1
                    #         text_file.write(q+'\t'+r+'\n')

                    #         if selection_indices==[]:
                    #             break
                    if final_counts[r]<d3[r]:
                        final_counts[r]+=1
                        text_file.write(q+'\t'+r+'\n')

                    count+=1

plot_normalized_read_distribution(final_counts, 'Downsampled Read Distribution', 'long_ds.png')
print('Downsampling Complete')

