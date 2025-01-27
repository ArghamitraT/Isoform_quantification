#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~ #
# Standard library imports
from collections import *

# Third party imports
import pysam
import pandas as pd
from tqdm import tqdm

# Local imports
from NanoCount.Read import Read
from NanoCount.common import *

# ~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~ #
class ExtractStatistics:

    # ~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~ #
    def __init__(
        self,
        alignment_file: str="",
        filter_bam_out: str = "",
        min_alignment_length: int = 50,
        keep_suplementary: bool = False,
        min_query_fraction_aligned: float = 0.5,
        sec_scoring_threshold: float = 0.95,
        sec_scoring_value: str = "alignment_score",
        convergence_target: float = 0.005,
        max_em_rounds: int = 100,
        extra_tx_info: bool = False,
        primary_score: str = "alignment_score",
        max_dist_3_prime: int = 50,
        max_dist_5_prime: int = -1,
        verbose: bool = False,
        quiet: bool = False,
    ):

        # Init package
        opt_summary_dict = opt_summary(local_opt=locals())
        self.log = get_logger(name="ExtractStatistics", verbose=verbose, quiet=quiet)

        self.log.warning("Checking options and input files")
        log_dict(opt_summary_dict, self.log.debug, "Options summary")

        # Save args in self variables
        self.alignment_file = alignment_file
        self.filter_bam_out = filter_bam_out
        self.min_alignment_length = min_alignment_length
        self.min_query_fraction_aligned = min_query_fraction_aligned
        self.sec_scoring_threshold = sec_scoring_threshold
        self.sec_scoring_value = sec_scoring_value
        self.convergence_target = convergence_target
        self.max_em_rounds = max_em_rounds
        self.extra_tx_info = extra_tx_info
        self.primary_score = primary_score
        self.keep_suplementary = keep_suplementary
        self.max_dist_5_prime = max_dist_5_prime
        self.max_dist_3_prime = max_dist_3_prime


        # Collect all alignments grouped by read name
        self.log.info("Parse Bam file and filter low quality alignments")

        #Getting reads and isoforms
        read_dict, ref_len_dict = self._parse_bam(read='short')
        self.read_dict = read_dict
        #self.all_reads_dict = all_reads_dict

        #print(read_dict)
        

        compatibility_dict = self.get_compatibility_modified() #(AT)
        #print(compatibility_dict)
        
        ambiguous_mappings = 0
        total_mappings = 0
        for i in compatibility_dict.keys():
            total_mappings+=1
            if len(compatibility_dict[i])>1: #This means that there is ambiguity (read maps to more than one isoform)
                ambiguous_mappings +=1
        print('Total number of reads = ', total_mappings)
        print('Number of ambiguous reads = ', ambiguous_mappings)
        print('Fraction of ambigous reads = ', ambiguous_mappings/total_mappings)

    
    def is_iterable(self, obj):
        """
        Check if obj is iterable but not a string.
        """
        return isinstance(obj, (list, tuple, set, dict, frozenset)) and not isinstance(obj, (str, bytes))

    def get_compatibility_modified(self):
        """"""
        compatibility_dict = defaultdict(dict)

        #for read_name, read in self.read_dict.items():
        for read_name, read in self.read_dict.items():
            # Check if read is iterable and not a string; if not, make it a list
            if not self.is_iterable(read):
                read = [read]  # Wrap non-iterable read in a list

            for alignment in read:
                for string in alignment.alignment_list:
                    score = 1.0 / alignment.n_alignment
                    # compatibility_dict[read_name][string.rname] = \
                    #     ((score * string.align_len)/self.ref_len_dict[string.rname])
                    compatibility_dict[read_name][string.rname] = score

        return compatibility_dict

    # ~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~ #
    
    def _parse_bam(self, read='short'):
        """
        Parse Bam/Sam file, group alignments per reads, filter reads based on
        selection criteria and return a dict of valid read/alignments
        """
        # Parse bam files
        read_dict = defaultdict(Read)
        ref_len_dict = OrderedDict()
        c = Counter()

        aligned_read = self.alignment_file

        # (AT)
        # with pysam.AlignmentFile(self.alignment_file) as bam:
        with pysam.AlignmentFile(aligned_read) as bam:

            # Collect reference lengths in dict
            for name, length in zip(bam.references, bam.lengths):
                ref_len_dict[name] = length

            for idx, alignment in enumerate(bam):
                if alignment.is_unmapped:
                    c["Discarded unmapped alignments"] += 1
                elif alignment.is_reverse:
                    c["Discarded negative strand alignments"] += 1
                # elif not self.keep_suplementary and alignment.is_supplementary:
                #     c["Discarded supplementary alignments"] += 1
                # elif self.min_alignment_length > 0 and alignment.query_alignment_length < self.min_alignment_length:
                #     c["Discarded short alignments"] += 1
                # elif self.max_dist_3_prime >= 0 and alignment.reference_end <= ref_len_dict[alignment.reference_name] - self.max_dist_3_prime:
                #     c["Discarded alignment with invalid 3 prime end"] += 1
                # elif self.max_dist_5_prime >= 0 and alignment.reference_start >= self.max_dist_5_prime:
                #     c["Discarded alignment with invalid 5 prime end"] += 1
                else:
                    c["Valid alignments"] += 1
                    read_dict[alignment.query_name].add_pysam_alignment(pysam_aligned_segment=alignment, read_idx=idx)

        # Write filtered reads counters
        log_dict(
            d=c,
            logger=self.log.info,
            header="Summary of alignments parsed in input bam file",
        )

        # Filter alignments
        filtered_read_dict = defaultdict(Read)

        #dictionary containing primary and secondary alignments
        #all_alignments_dict = defaultdict(Read)

        c = Counter()
        for query_name, read in read_dict.items():
            # Check if best alignment is valid
            best_alignment = read.get_best_alignment(primary_score=self.primary_score)

            # In case the primary alignment was removed by filters
            if best_alignment:
                if best_alignment.align_score == 0:
                    c["Reads with zero score"] += 1
                if best_alignment.align_len == 0:
                    c["Reads with zero len"] += 1
                elif best_alignment.query_fraction_aligned < self.min_query_fraction_aligned:
                    c["Reads with low query fraction aligned"] += 1
                else:
                    filtered_read_dict[query_name].add_alignment(best_alignment)
                    #all_alignments_dict[query_name].add_alignment(best_alignment)
                    c["Reads with valid best alignment"] += 1
                    for alignment in read.get_secondary_alignments_list(primary_score=self.primary_score):
                        
                        #Get all secondary alignments, not just "high quality" ones
                        #all_alignments_dict[query_name].add_alignment(alignment)

                        # Filter out secondary alignments based on minimap alignment score
                        if self.sec_scoring_value == "alignment_score" and alignment.align_score / best_alignment.align_score < self.sec_scoring_threshold:
                            c["Invalid secondary alignments"] += 1

                        # Filter out secondary alignments based on minimap alignment length
                        elif self.sec_scoring_value == "alignment_length" and alignment.align_len / best_alignment.align_len < self.sec_scoring_threshold:
                            c["Invalid secondary alignments"] += 1

                        # Select valid secondary alignments
                        else:
                            c["Valid secondary alignments"] += 1
                            filtered_read_dict[query_name].add_alignment(alignment)
            else:
                c["Reads without best alignment"] += 1

        if not "Valid secondary alignments" in c:
            self.log.error("No valid secondary alignments found in bam file. Were the reads aligned with minimap `-p 0 -N 10` options ?")

        # Write filtered reads counters
        log_dict(d=c, logger=self.log.info, header="Summary of reads filtered")

        # (AT)
        return filtered_read_dict, ref_len_dict #, all_alignments_dict


ExtractStatistics(alignment_file=sys.argv[1])