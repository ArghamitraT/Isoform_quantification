#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import Counter
from collections import *
import numpy as np
import matplotlib.pyplot as plt
import sys

import pysam
import pandas as pd
from tqdm import tqdm

from NanoCount.Read import Read
from NanoCount.common import *

bam_path = sys.argv[1]
img_path = sys.argv[2]

def parse_bam(bam_file):
        """
        Parse Bam/Sam file, group alignments per reads, filter reads based on
        selection criteria and return a dict of valid read/alignments
        """
        # Parse bam files
        read_dict = defaultdict(Read)
        ref_len_dict = OrderedDict()
        c = Counter()

        aligned_read = bam_file

        #Default Initializations - NanoCount

        min_alignment_length = 50
        keep_suplementary = False
        min_query_fraction_aligned = 0.5
        sec_scoring_threshold = 0.95
        sec_scoring_value = "alignment_score"
        primary_score = "alignment_score"
        max_dist_3_prime = 50
        max_dist_5_prime = -1

        # (AT)
        # with pysam.AlignmentFile(alignment_file) as bam:
        with pysam.AlignmentFile(aligned_read) as bam:

            # Collect reference lengths in dict
            for name, length in zip(bam.references, bam.lengths):
                ref_len_dict[name] = length

            for idx, alignment in enumerate(bam):
                if alignment.is_unmapped:
                    c["Discarded unmapped alignments"] += 1
                #elif alignment.is_reverse:
                #    c["Discarded negative strand alignments"] += 1
                # elif not keep_suplementary and alignment.is_supplementary:
                #     c["Discarded supplementary alignments"] += 1
                # elif min_alignment_length > 0 and alignment.query_alignment_length < min_alignment_length:
                #     c["Discarded short alignments"] += 1
                # elif max_dist_3_prime >= 0 and alignment.reference_end <= ref_len_dict[alignment.reference_name] - max_dist_3_prime:
                #     c["Discarded alignment with invalid 3 prime end"] += 1
                # elif max_dist_5_prime >= 0 and alignment.reference_start >= max_dist_5_prime:
                #     c["Discarded alignment with invalid 5 prime end"] += 1
                else:
                    c["Valid alignments"] += 1
                    read_dict[alignment.query_name].add_pysam_alignment(pysam_aligned_segment=alignment, read_idx=idx)


        # Filter alignments
        filtered_read_dict = defaultdict(Read)

        #dictionary containing primary and secondary alignments
        #all_alignments_dict = defaultdict(Read)

        c = Counter()
        for query_name, read in read_dict.items():
            # Check if best alignment is valid
            best_alignment = read.get_best_alignment(primary_score=primary_score)

            # In case the primary alignment was removed by filters
            if best_alignment:
                if best_alignment.align_score == 0:
                    c["Reads with zero score"] += 1
                if best_alignment.align_len == 0:
                    c["Reads with zero len"] += 1
                elif best_alignment.query_fraction_aligned < min_query_fraction_aligned:
                    c["Reads with low query fraction aligned"] += 1
                else:
                    filtered_read_dict[query_name].add_alignment(best_alignment)
                    #all_alignments_dict[query_name].add_alignment(best_alignment)
                    c["Reads with valid best alignment"] += 1
                    for alignment in read.get_secondary_alignments_list(primary_score=primary_score):
                        
                        #Get all secondary alignments, not just "high quality" ones
                        #all_alignments_dict[query_name].add_alignment(alignment)

                        # Filter out secondary alignments based on minimap alignment score
                        if sec_scoring_value == "alignment_score" and alignment.align_score / best_alignment.align_score < sec_scoring_threshold:
                            c["Invalid secondary alignments"] += 1

                        # Filter out secondary alignments based on minimap alignment length
                        elif sec_scoring_value == "alignment_length" and alignment.align_len / best_alignment.align_len < sec_scoring_threshold:
                            c["Invalid secondary alignments"] += 1

                        # Select valid secondary alignments
                        else:
                            c["Valid secondary alignments"] += 1
                            filtered_read_dict[query_name].add_alignment(alignment)
            else:
                c["Reads without best alignment"] += 1

        return filtered_read_dict, ref_len_dict

def read_dict_to_counter(read_dict):
    counts = defaultdict(int)

    for r in list(read_dict.keys()):
        for a in read_dict[r].alignment_list:
            counts[a.rname]+=1

    #print(type(read_dict['m54284U_210712_004945/6947739/ccs'].alignment_list[0].rname))

    return Counter(counts)

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


def plot_read_distribution(sample, title, nameFig):
    # Assuming 'sample' is a dictionary and 'title' and 'nameFig' are defined
    max_key = max(sample, key=sample.get)
    max_value = sample[max_key]

    # Prepare the data
    sorted_items = sorted(sample.items(), key=lambda item: item[1])
    mid_index = len(sorted_items) // 2
    left_half = sorted_items[:mid_index]
    right_half = sorted_items[mid_index:]
    ordered_keys = [k for k, v in left_half] + [max_key] + [k for k, v in reversed(right_half)]
    ordered_values = [sample[k] for k in ordered_keys]

    # Set the figure size and bar width
    fig_width = 20  # You might need to adjust this based on your actual data
    bar_width = fig_width / len(ordered_values)  # Adjust bar width based on number of data points

    # Create the figure and plot the bars
    plt.figure(figsize=(fig_width, 5))
    plt.bar(range(len(ordered_values)), ordered_values, width=bar_width)

    # Set a logarithmic scale if the data varies widely
    plt.yscale('log')

    # Set the y-axis limit to slightly above the max value
    plt.ylim(0.1, max_value * 10)  # Adjust the lower limit if you have zero or negative values

    # Add a title and other plot settings
    plt.title(title)
    plt.tight_layout()

    # Save the figure to a file and show it
    plt.savefig('figures/' + nameFig)
    plt.show()


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

with open(bam_path, 'rb') as file:
     read_dict, len_dict = parse_bam(file) #pickle.load(file)
counts = read_dict_to_counter(read_dict)
print(counts)
plot_normalized_read_distribution(counts, 'Read Distribution for Sample', img_path)


# with open('theta_PacBio_day0_rep1.pkl', 'rb') as file:
#     theta_sample2 = pickle.load(file)

# # Find common isoforms between the two samples
# common_isoforms_dict = {isoform: theta_sample2[isoform] for isoform in theta_sample1 if isoform in theta_sample2}
# common_isoforms = Counter(common_isoforms_dict)


# Plot read distribution for theta_sample1
# plot_read_distribution(theta_sample1, 'Read Distribution for Sample 1', 'theta_sample1.png')

# Plot read distribution for theta_sample2
# plot_normalized_read_distribution(theta_sample2, 'Read Distribution for Sample 2', 'theta_sample2.png')

# Plot overlapping read distribution for common isoforms
#plot_normalized_read_distribution(common_isoforms, 'Overlapping Small Read Distribution', 'small.png')

# Display the plot

#file_path = 'PacBio/PacBio/day0_rep1/collapsed_transcripts.flnc_count.txt'
# Replace 'path_to_your_text_file.txt' with the actual path to your text file
#counter_from_text_file = convert_text_to_counter(file_path)

#plot_normalized_read_distribution(counter_from_text_file, 'Long Read Distribution', 'long.png')

