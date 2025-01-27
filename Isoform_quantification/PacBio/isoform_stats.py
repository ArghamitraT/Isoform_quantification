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

bam_path_short = sys.argv[1]
bam_path_long = sys.argv[2]

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
                elif alignment.is_reverse:
                    c["Discarded negative strand alignments"] += 1
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
                        #if sec_scoring_value == "alignment_score" and alignment.align_score / best_alignment.align_score < sec_scoring_threshold:
                        #    c["Invalid secondary alignments"] += 1

                        # Filter out secondary alignments based on minimap alignment length
                        if sec_scoring_value == "alignment_length" and alignment.align_len / best_alignment.align_len < sec_scoring_threshold: #elif
                            c["Invalid secondary alignments"] += 1

                        # Select valid secondary alignments
                        else:
                            c["Valid secondary alignments"] += 1
                            filtered_read_dict[query_name].add_alignment(alignment)
            else:
                c["Reads without best alignment"] += 1

        return filtered_read_dict, ref_len_dict


with open(bam_path_short, 'rb') as file:
     read_dict_short, len_dict_short = parse_bam(file)

short_isoforms = set()
for r in list(read_dict_short.keys()):
    for a in read_dict_short[r].alignment_list:
        short_isoforms.add(a.rname)
count_short = len(short_isoforms)

with open(bam_path_long, 'rb') as file:
     read_dict_long, len_dict_long = parse_bam(file)

long_isoforms = set()
for r in list(read_dict_long.keys()):
    for a in read_dict_long[r].alignment_list:
        long_isoforms.add(a.rname)
count_long = len(long_isoforms)

common_isoforms = set(long_isoforms).intersection(set(short_isoforms))
count_common = len(common_isoforms)

print("Isoforms in short alignment file = ", count_short)
print("Isoforms in long alignment file = ", count_long)
print("Common Isoforms = ", count_common)
