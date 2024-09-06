#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~ #
# Standard library imports
from collections import *

# Third party imports
import pysam
import pandas as pd
from tqdm import tqdm
import hashlib
import pickle

# Local imports
from NanoCount.Read import Read
from NanoCount.common import *

# ~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~ #
class NanoCount:

    # ~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~ #
    def __init__(
        self,
        short_read_file: str,   #(AT)
        long_read_file: str,    #(AT)
        alignment_file: str="",
        count_file: str = "",
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
        """
        Estimate abundance of transcripts using an EM
        * alignment_file
            Sorted and indexed BAM or SAM file containing aligned ONT dRNA-Seq reads including secondary alignments
        * count_file
            Output file path where to write estimated counts (TSV format)
        * filter_bam_out
            Optional output file path where to write filtered reads selected by NanoCount to perform quantification estimation (BAM format)
        * min_alignment_length
            Minimal length of the alignment to be considered valid
        * min_query_fraction_aligned
            Minimal fraction of the primary alignment query aligned to consider the read valid
        * sec_scoring_threshold
            Fraction of the alignment score or the alignment length of secondary alignments compared to the primary alignment to be considered valid
            alignments
        * sec_scoring_value
            Value to use for score thresholding of secondary alignments either "alignment_score" or "alignment_length"
        * convergence_target
            Convergence target value of the cummulative difference between abundance values of successive EM round to trigger the end of the EM loop.
        * max_em_rounds
            Maximum number of EM rounds before triggering stop
        * extra_tx_info
            Add transcripts length and zero coverage transcripts to the output file (required valid bam/sam header)
        * primary_score
            Method to pick the best alignment for each read. By default ("alignment_score") uses the best alignment score (AS optional field), but it can be changed to
            use either the primary alignment defined by the aligner ("primary") or the longest alignment ("alignment_length"). choices = [primary, alignment_score, alignment_length]
        * keep_suplementary
            Retain any supplementary alignments and considered them like secondary alignments. Discarded by default.
        * max_dist_3_prime
            Maximum distance of alignment end to 3 prime of transcript. In ONT dRNA-Seq reads are assumed to start from the polyA tail (-1 to deactivate)
        * max_dist_5_prime
            Maximum distance of alignment start to 5 prime of transcript. In conjunction with max_dist_3_prime it can be used to select near full transcript reads
            only (-1 to deactivate).
        * verbose
            Increase verbosity for QC and debugging
        * quiet
            Reduce verbosity
        """

        # Init package
        opt_summary_dict = opt_summary(local_opt=locals())
        self.log = get_logger(name="Nanocount", verbose=verbose, quiet=quiet)

        self.log.warning("Checking options and input files")
        log_dict(opt_summary_dict, self.log.debug, "Options summary")

        # Save args in self variables
        self.short_read_file = short_read_file   # (AT)
        self.long_read_file = long_read_file     # (AT)
        self.alignment_file = alignment_file
        self.count_file = count_file
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

        self.log.warning("Initialise Nanocount")

        # Collect all alignments grouped by read name
        self.log.info("Parse Bam file and filter low quality alignments")

        # (AT) start
        # UNCOMMENT
        # read_dict_short, ref_len_dict_short = self._parse_bam(read='short')
        # read_dict_long, ref_len_dict_long = self._parse_bam(read='long')

        # COMMENT
        with open('pkl_files/read_dict_short.pkl', 'rb') as file:
            read_dict_short = pickle.load(file)
        with open('pkl_files/read_dict_long.pkl', 'rb') as file:
            read_dict_long = pickle.load(file)
        with open('pkl_files/ref_len_dict_short.pkl', 'rb') as file:
            ref_len_dict_short = pickle.load(file)
        with open('pkl_files/ref_len_dict_long.pkl', 'rb') as file:
            ref_len_dict_long = pickle.load(file)

        # Initialize merged_dict with an appropriate default factory
        self.read_dict = defaultdict()
        # Add all items from the first dictionary to self.read_dict
        for key, value in read_dict_long.items():
            self.read_dict[key] = value

        # Merge with items from the second dictionary, combining values if key exists
        for key, value in read_dict_short.items():
            if key in self.read_dict:
                # Combine values; this could be a list, a sum, etc., depending on your needs; Here, we're creating a list of values
                if not isinstance(self.read_dict[key], list):
                    self.read_dict[key] = [self.read_dict[key]]
                self.read_dict[key].append(value)
            else:
                # If the key is not in self.read_dict, simply add it
                self.read_dict[key] = value

        # merging the transcript length for both long and short read
        self.ref_len_dict = OrderedDict()
        self.ref_len_dict = {k: v for k, v in ref_len_dict_long.items() if k not in ref_len_dict_short}
        self.ref_len_dict.update(ref_len_dict_short)

        # Add all items from the first dictionary to self.read_dict
        for key, value in ref_len_dict_long.items():
            self.ref_len_dict[key] = value
        # (AT) end

        if self.filter_bam_out:
            self.log.info("Write selected alignments to BAM file")
            self._write_bam()

        # Generate compatibility dict grouped by reads
        self.log.info("Generate initial read/transcript compatibility index")
        compatibility_dict = self.get_compatibility_modified() #(AT)

        # (AT) start edit
        hashed_compatibility_dict = defaultdict(dict)
        # Iterate over the compatibility_dict to populate hashed_compatibility_dict
        for key, value in compatibility_dict.items():
            # Use the values (which are dictionaries) to create a hash key
            value_tuple = tuple(value.items())  # This creates a hashable tuple of the dictionary items
            hashed_key = self.hash_key(value_tuple)

            # Merge dictionaries with the same hash, or add them if not present
            if hashed_key in hashed_compatibility_dict:
                # Update the existing dictionary with the new values, if they don't already exist
                for subkey, subvalue in value.items():
                    if subkey not in hashed_compatibility_dict[hashed_key]:
                        hashed_compatibility_dict[hashed_key][subkey] = subvalue
            else:
                hashed_compatibility_dict[hashed_key] = value
        self.compatibility_dict =  hashed_compatibility_dict
        # (AT) end edit

        # EM loop to calculate abundance and update read-transcript compatibility
        self.log.warning("Start EM abundance estimate")

        self.em_round = 0
        self.convergence = 1

        with tqdm(
            unit=" rounds",
            unit_scale=True,
            desc="\tProgress",
            disable=(quiet or verbose),
        ) as pbar:
            # Iterate until convergence threshold or max EM round are reached
            while self.convergence > self.convergence_target and self.em_round < self.max_em_rounds:
                self.em_round += 1
                # Calculate abundance from compatibility assignments
                self.abundance_dict = self._calculate_abundance()
                # Update compatibility assignments
                self.compatibility_dict = self._update_compatibility()
                # Update counter
                pbar.update(1)
                self.log.debug("EM Round: {} / Convergence value: {}".format(self.em_round, self.convergence))

        self.log.info("Exit EM loop after {} rounds".format(self.em_round))
        self.log.info("Convergence value: {}".format(self.convergence))
        if not self.convergence <= self.convergence_target:
            self.log.error("Convergence target ({}) could not be reached after {} rounds".format(self.convergence_target, self.max_em_rounds))

        # Write out results
        self.log.warning("Summarize data")

        self.log.info("Convert results to dataframe")
        self.count_df = pd.DataFrame(self.abundance_dict.most_common(), columns=["transcript_name", "raw"])
        self.count_df.set_index("transcript_name", inplace=True, drop=True)

        self.log.info("Compute estimated counts and TPM")
        self.count_df["est_count"] = self.count_df["raw"] * len(self.read_dict)
        self.count_df["tpm"] = self.count_df["raw"] * 1000000

        # Add extra transcript info is required
        if self.extra_tx_info:
            tx_df = self._get_tx_df()
            self.count_df = pd.merge(self.count_df, tx_df, left_index=True, right_index=True, how="outer")

        # Cleanup and sort
        self.count_df.sort_values(by="raw", ascending=False, inplace=True)
        self.count_df.fillna(value=0, inplace=True)
        self.count_df.index.name = "transcript_name"

        if self.count_file:
            self.log.info("Write file")
            self.count_df.to_csv(self.count_file, sep="\t")

    # (AT)
    def hash_key(self, elements):
        # Create a stable hash of the sorted elements using a cryptographic hash function
        # Convert the elements tuple to a string and encode to bytes before hashing
        return hashlib.sha256(str(sorted(elements)).encode('utf-8')).hexdigest()

    # (AT)
    def is_iterable(self, obj):
        """
        Check if obj is iterable but not a string.
        """
        return isinstance(obj, (list, tuple, set, dict, frozenset)) and not isinstance(obj, (str, bytes))

    def get_compatibility_modified(self):
        """"""
        compatibility_dict = defaultdict(dict)

        for read_name, read in self.read_dict.items():
            # Check if read is iterable and not a string; if not, make it a list
            if not self.is_iterable(read):
                read = [read]  # Wrap non-iterable read in a list

            for alignment in read:
                for string in alignment.alignment_list:
                    score = 1.0 / alignment.n_alignment
                    compatibility_dict[read_name][string.rname] = \
                        ((score * string.align_len)/self.ref_len_dict[string.rname])

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

        # (AT) which file to align (can be optimized)
        if read == 'short':
            aligned_read = self.short_read_file
        else:
            aligned_read = self.long_read_file

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
                elif not self.keep_suplementary and alignment.is_supplementary:
                    c["Discarded supplementary alignments"] += 1
                elif self.min_alignment_length > 0 and alignment.query_alignment_length < self.min_alignment_length:
                    c["Discarded short alignments"] += 1
                elif self.max_dist_3_prime >= 0 and alignment.reference_end <= ref_len_dict[alignment.reference_name] - self.max_dist_3_prime:
                    c["Discarded alignment with invalid 3 prime end"] += 1
                elif self.max_dist_5_prime >= 0 and alignment.reference_start >= self.max_dist_5_prime:
                    c["Discarded alignment with invalid 5 prime end"] += 1
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
        c = Counter()

        for query_name, read in read_dict.items():
            # Check if best alignment is valid
            best_alignment = read.get_best_alignment(primary_score=self.primary_score)

            # In case the primary alignment was removed by filters
            if best_alignment:
                if best_alignment.align_score == 0:
                    c["Reads with zero score"] += 1
                elif best_alignment.align_len == 0:
                    c["Reads with zero len"] += 1
                elif best_alignment.query_fraction_aligned < self.min_query_fraction_aligned:
                    c["Reads with low query fraction aligned"] += 1
                else:
                    filtered_read_dict[query_name].add_alignment(best_alignment)
                    c["Reads with valid best alignment"] += 1
                    for alignment in read.get_secondary_alignments_list(primary_score=self.primary_score):

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
        return filtered_read_dict, ref_len_dict

    def _write_bam(self):
        """"""
        c = Counter()

        # Make list of alignments idx to select
        selected_read_idx = set()
        for read in self.read_dict.values():
            for alignment in read.alignment_list:
                selected_read_idx.add(alignment.read_idx)
                c["Alignments to select"] += 1

        # Select from original bam file and write to output bam file
        with pysam.AlignmentFile(self.alignment_file) as bam_in:
            with pysam.AlignmentFile(self.filter_bam_out, "wb", template=bam_in) as bam_out:
                for read_idx, alignment in enumerate(bam_in):
                    if read_idx in selected_read_idx:
                        bam_out.write(alignment)
                        c["Alignments written"] += 1
                    else:
                        c["Alignments skipped"] += 1

        log_dict(d=c, logger=self.log.info, header="Summary of alignments written to bam")

    
    def _get_compatibility(self):
        """"""
        compatibility_dict = defaultdict(dict)

        for read_name, read in self.read_dict.items():
            # Check if read is iterable and not a string; if not, make it a list
            if not self.is_iterable(read):
                read = [read]  # Wrap non-iterable read in a list

            for alignment in read:
                for string in alignment.alignment_list:
                    compatibility_dict[read_name][string.rname] = \
                        score = 1.0 / alignment.n_alignment

        return compatibility_dict

    def _calculate_abundance(self):
        """
        Calculate the abundance of the transcript set based on read-transcript compatibilities
        """
        abundance_dict = Counter()
        total = 0
        convergence = 0

        for read_name, comp in self.compatibility_dict.items():
            for ref_name, score in comp.items():
                abundance_dict[ref_name] += score
                total += score

        for ref_name in abundance_dict.keys():
            abundance_dict[ref_name] = abundance_dict[ref_name] / total

            if self.em_round > 1:
                convergence += abs(self.abundance_dict[ref_name] - abundance_dict[ref_name])

        if self.em_round == 1:
            self.convergence = 1
        else:
            self.convergence = convergence

        return abundance_dict

    def _update_compatibility(self):
        """
        Update read-transcript compatibility based on transcript abundances
        """
        compatibility_dict = defaultdict(dict)

        for read_name, comp in self.compatibility_dict.items():
            total = 0
            for ref_name in comp.keys():
                total += self.abundance_dict[ref_name]

            for ref_name in comp.keys():
                compatibility_dict[read_name][ref_name] = self.abundance_dict[ref_name] / total

        return compatibility_dict

    def _get_tx_df(self):
        """
        Extract transcript info from bam file header
        """
        try:
            with pysam.AlignmentFile(self.alignment_file) as bam:
                references = bam.references
                lengths = bam.lengths

            return pd.DataFrame(index=references, data=lengths, columns=["transcript_length"])
        # If any error return empty DataFrame silently
        except Exception:
            return pd.DataFrame()




















#####################################################
# """ ORIGINAL NanoCount code """

# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# # ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~ #
# # Standard library imports
# from collections import *

# # Third party imports
# import pysam
# import pandas as pd
# from tqdm import tqdm
# import time

# # Local imports
# from NanoCount.Read import Read
# from NanoCount.common import *

# # ~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~ #

# class Read_new:
#     def __init__(self):
#         self.alignment_list = []

#     def __repr__(self):
#         m = ""
#         for r in self.alignment_list:
#             m += "\t\t{}\n".format(r)
#         return m

#     def add_alignment(self, alignment):
#         self.alignment_list.append(alignment)

#     def to_json(self):
#         return self.__repr__()

# # Define the Alignment class
# class Alignment:
#     def __init__(self, q, r, qlen=0, align_len=0, align_score=0, secondary=False):
#         self.qname = q
#         self.rname = r
#         self.qlen = qlen
#         self.align_len = align_len
#         self.align_score = align_score
#         self.secondary = secondary

#     def __repr__(self):
#         return f"Query:{self.qname} | Reference:{self.rname} | Query len:{self.qlen} | Alignment len:{self.align_len} | Align Score:{self.align_score} | Secondary:{self.secondary}"

# class EnhancedOrderedDict(OrderedDict):
#     def add(self, key, value=None):
#         """Adds a key with a default value if not specified."""
#         if key not in self:
#             self[key] = value


# class NanoCount:

#     # ~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~ #
#     def __init__(
#         self,
#         alignment_file: str,
#         count_file: str = "",
#         filter_bam_out: str = "",
#         min_alignment_length: int = 50,
#         keep_suplementary: bool = False,
#         min_query_fraction_aligned: float = 0.5,
#         sec_scoring_threshold: float = 0.95,
#         sec_scoring_value: str = "alignment_score",
#         convergence_target: float = 0.005,
#         max_em_rounds: int = 100,
#         extra_tx_info: bool = False,
#         primary_score: str = "alignment_score",
#         max_dist_3_prime: int = 50,
#         max_dist_5_prime: int = -1,
#         verbose: bool = False,
#         quiet: bool = False,
#         downsampled_files="",

#     ):
#         """
#         Estimate abundance of transcripts using an EM
#         * alignment_file
#             Sorted and indexed BAM or SAM file containing aligned ONT dRNA-Seq reads including secondary alignments
#         * count_file
#             Output file path where to write estimated counts (TSV format)
#         * filter_bam_out
#             Optional output file path where to write filtered reads selected by NanoCount to perform quantification estimation (BAM format)
#         * min_alignment_length
#             Minimal length of the alignment to be considered valid
#         * min_query_fraction_aligned
#             Minimal fraction of the primary alignment query aligned to consider the read valid
#         * sec_scoring_threshold
#             Fraction of the alignment score or the alignment length of secondary alignments compared to the primary alignment to be considered valid
#             alignments
#         * sec_scoring_value
#             Value to use for score thresholding of secondary alignments either "alignment_score" or "alignment_length"
#         * convergence_target
#             Convergence target value of the cummulative difference between abundance values of successive EM round to trigger the end of the EM loop.
#         * max_em_rounds
#             Maximum number of EM rounds before triggering stop
#         * extra_tx_info
#             Add transcripts length and zero coverage transcripts to the output file (required valid bam/sam header)
#         * primary_score
#             Method to pick the best alignment for each read. By default ("alignment_score") uses the best alignment score (AS optional field), but it can be changed to
#             use either the primary alignment defined by the aligner ("primary") or the longest alignment ("alignment_length"). choices = [primary, alignment_score, alignment_length]
#         * keep_suplementary
#             Retain any supplementary alignments and considered them like secondary alignments. Discarded by default.
#         * max_dist_3_prime
#             Maximum distance of alignment end to 3 prime of transcript. In ONT dRNA-Seq reads are assumed to start from the polyA tail (-1 to deactivate)
#         * max_dist_5_prime
#             Maximum distance of alignment start to 5 prime of transcript. In conjunction with max_dist_3_prime it can be used to select near full transcript reads
#             only (-1 to deactivate).
#         * verbose
#             Increase verbosity for QC and debugging
#         * quiet
#             Reduce verbosity
#         """

#         # Init package
#         opt_summary_dict = opt_summary(local_opt=locals())
#         self.log = get_logger(name="Nanocount", verbose=verbose, quiet=quiet)

#         self.log.warning("Checking options and input files")
#         log_dict(opt_summary_dict, self.log.debug, "Options summary")

#         # Save args in self variables
#         self.alignment_file = alignment_file
#         self.count_file = count_file
#         self.filter_bam_out = filter_bam_out
#         self.min_alignment_length = min_alignment_length
#         self.min_query_fraction_aligned = min_query_fraction_aligned
#         self.sec_scoring_threshold = sec_scoring_threshold
#         self.sec_scoring_value = sec_scoring_value
#         self.convergence_target = convergence_target
#         self.max_em_rounds = max_em_rounds
#         self.extra_tx_info = extra_tx_info
#         self.primary_score = primary_score
#         self.keep_suplementary = keep_suplementary
#         self.max_dist_5_prime = max_dist_5_prime
#         self.max_dist_3_prime = max_dist_3_prime

#         self.downsampled_files = downsampled_files
#         # self.temp_read_dicts = defaultdict(Read)
#         # self.read_dict = defaultdict(Read) # (AT)

#         self.log.warning("Initialise Nanocount")

#         # Collect all alignments grouped by read name
#         self.log.info("Parse Bam file and filter low quality alignments")
#         start = time.time()
#         self.read_dict = self._parse_bam()  # (AT)
#         # self.load_data()    # (AT)
#         end = time.time()
#         print(f"Time taken to run the parsing was {end - start} seconds")



#         if self.filter_bam_out:
#             self.log.info("Write selected alignments to BAM file")
#             self._write_bam()

#         # Generate compatibility dict grouped by reads
#         self.log.info("Generate initial read/transcript compatibility index")
#         self.compatibility_dict = self._get_compatibility()

#         # EM loop to calculate abundance and update read-transcript compatibility
#         self.log.warning("Start EM abundance estimate")

#         self.em_round = 0
#         self.convergence = 1

#         with tqdm(
#             unit=" rounds",
#             unit_scale=True,
#             desc="\tProgress",
#             disable=(quiet or verbose),
#         ) as pbar:
#             # Iterate until convergence threshold or max EM round are reached
#             while self.convergence > self.convergence_target and self.em_round < self.max_em_rounds:
#                 self.em_round += 1
#                 # Calculate abundance from compatibility assignments
#                 self.abundance_dict = self._calculate_abundance()
#                 # Update compatibility assignments
#                 self.compatibility_dict = self._update_compatibility()
#                 # Update counter
#                 pbar.update(1)
#                 self.log.debug("EM Round: {} / Convergence value: {}".format(self.em_round, self.convergence))

#         self.log.info("Exit EM loop after {} rounds".format(self.em_round))
#         self.log.info("Convergence value: {}".format(self.convergence))
#         if not self.convergence <= self.convergence_target:
#             self.log.error("Convergence target ({}) could not be reached after {} rounds".format(self.convergence_target, self.max_em_rounds))

#         # Write out results
#         self.log.warning("Summarize data")

#         self.log.info("Convert results to dataframe")
#         self.count_df = pd.DataFrame(self.abundance_dict.most_common(), columns=["transcript_name", "raw"])
#         self.count_df.set_index("transcript_name", inplace=True, drop=True)

#         self.log.info("Compute estimated counts and TPM")
#         self.count_df["est_count"] = self.count_df["raw"] * len(self.read_dict)
#         self.count_df["tpm"] = self.count_df["raw"] * 1000000

#         # Add extra transcript info is required
#         if self.extra_tx_info:
#             tx_df = self._get_tx_df()
#             self.count_df = pd.merge(self.count_df, tx_df, left_index=True, right_index=True, how="outer")

#         # Cleanup and sort
#         self.count_df.sort_values(by="raw", ascending=False, inplace=True)
#         self.count_df.fillna(value=0, inplace=True)
#         self.count_df.index.name = "transcript_name"

#         if self.count_file:
#             self.log.info("Write file")
#             self.count_df.to_csv(self.count_file, sep="\t")

#     def load_data(self):
#         with open(self.downsampled_files) as f:
#             header = f.readline()
#             for line in f:
#                 line = line.strip('\n')
#                 q, r = line.split('\t')
#                 aln = Alignment(q, r)
#                 #self.temp_read_dicts[q].add_alignment(aln)
#                 self.read_dict[q].add_alignment(aln)
#                 #self.temp_read_lengths.add((r, 0))
#         # self.temp_read_lengths = OrderedDict(list(self.temp_read_lengths))
#         # self.all_read_dicts['sample2'] = self.temp_read_dicts
#         # self.all_ref_len_dicts['sample2'] = self.temp_read_lengths
#         print('Data Loaded and Processed')


#     # ~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~ #
#     def _parse_bam(self):
#         """
#         Parse Bam/Sam file, group alignments per reads, filter reads based on
#         selection criteria and return a dict of valid read/alignments
#         """
#         # Parse bam files
#         read_dict = defaultdict(Read)
#         ref_len_dict = OrderedDict()
#         c = Counter()
#         with pysam.AlignmentFile(self.alignment_file) as bam:

#             # Collect reference lengths in dict
#             for name, length in zip(bam.references, bam.lengths):
#                 ref_len_dict[name] = length

#             for idx, alignment in enumerate(bam):
#                 if alignment.is_unmapped:
#                     c["Discarded unmapped alignments"] += 1
#                 elif alignment.is_reverse:
#                     c["Discarded negative strand alignments"] += 1
#                 elif not self.keep_suplementary and alignment.is_supplementary:
#                     c["Discarded supplementary alignments"] += 1
#                 elif self.min_alignment_length > 0 and alignment.query_alignment_length < self.min_alignment_length:
#                     c["Discarded short alignments"] += 1
#                 elif self.max_dist_3_prime >= 0 and alignment.reference_end <= ref_len_dict[alignment.reference_name] - self.max_dist_3_prime:
#                     c["Discarded alignment with invalid 3 prime end"] += 1
#                 elif self.max_dist_5_prime >= 0 and alignment.reference_start >= self.max_dist_5_prime:
#                     c["Discarded alignment with invalid 5 prime end"] += 1
#                 else:
#                     c["Valid alignments"] += 1
#                     read_dict[alignment.query_name].add_pysam_alignment(pysam_aligned_segment=alignment, read_idx=idx)

#         # Write filtered reads counters
#         log_dict(
#             d=c,
#             logger=self.log.info,
#             header="Summary of alignments parsed in input bam file",
#         )

#         # Filter alignments
#         filtered_read_dict = defaultdict(Read)
#         c = Counter()

#         for query_name, read in read_dict.items():
#             # Check if best alignment is valid
#             best_alignment = read.get_best_alignment(primary_score=self.primary_score)

#             # In case the primary alignment was removed by filters
#             if best_alignment:
#                 if best_alignment.align_score == 0:
#                     c["Reads with zero score"] += 1
#                 elif best_alignment.align_len == 0:
#                     c["Reads with zero len"] += 1
#                 elif best_alignment.query_fraction_aligned < self.min_query_fraction_aligned:
#                     c["Reads with low query fraction aligned"] += 1
#                 else:
#                     filtered_read_dict[query_name].add_alignment(best_alignment)
#                     c["Reads with valid best alignment"] += 1
#                     for alignment in read.get_secondary_alignments_list(primary_score=self.primary_score):

#                         # Filter out secondary alignments based on minimap alignment score
#                         if self.sec_scoring_value == "alignment_score" and alignment.align_score / best_alignment.align_score < self.sec_scoring_threshold:
#                             c["Invalid secondary alignments"] += 1

#                         # Filter out secondary alignments based on minimap alignment length
#                         elif self.sec_scoring_value == "alignment_length" and alignment.align_len / best_alignment.align_len < self.sec_scoring_threshold:
#                             c["Invalid secondary alignments"] += 1

#                         # Select valid secondary alignments
#                         else:
#                             c["Valid secondary alignments"] += 1
#                             filtered_read_dict[query_name].add_alignment(alignment)
#             else:
#                 c["Reads without best alignment"] += 1

#         if not "Valid secondary alignments" in c:
#             self.log.error("No valid secondary alignments found in bam file. Were the reads aligned with minimap `-p 0 -N 10` options ?")

#         # Write filtered reads counters
#         log_dict(d=c, logger=self.log.info, header="Summary of reads filtered")
#         return filtered_read_dict

#     def _write_bam(self):
#         """"""
#         c = Counter()

#         # Make list of alignments idx to select
#         selected_read_idx = set()
#         for read in self.read_dict.values():
#             for alignment in read.alignment_list:
#                 selected_read_idx.add(alignment.read_idx)
#                 c["Alignments to select"] += 1

#         # Select from original bam file and write to output bam file
#         with pysam.AlignmentFile(self.alignment_file) as bam_in:
#             with pysam.AlignmentFile(self.filter_bam_out, "wb", template=bam_in) as bam_out:
#                 for read_idx, alignment in enumerate(bam_in):
#                     if read_idx in selected_read_idx:
#                         bam_out.write(alignment)
#                         c["Alignments written"] += 1
#                     else:
#                         c["Alignments skipped"] += 1

#         log_dict(d=c, logger=self.log.info, header="Summary of alignments written to bam")

#     def _get_compatibility(self):
#         """"""
#         compatibility_dict = defaultdict(dict)
#         for read_name, read in self.read_dict.items():
#             for alignment in read.alignment_list:
#                 compatibility_dict[read_name][alignment.rname] = score = 1.0 / read.n_alignment

#         return compatibility_dict

#     def _calculate_abundance(self):
#         """
#         Calculate the abundance of the transcript set based on read-transcript compatibilities
#         """
#         abundance_dict = Counter()
#         total = 0
#         convergence = 0

#         for read_name, comp in self.compatibility_dict.items():
#             for ref_name, score in comp.items():
#                 abundance_dict[ref_name] += score
#                 total += score

#         for ref_name in abundance_dict.keys():
#             abundance_dict[ref_name] = abundance_dict[ref_name] / total

#             if self.em_round > 1:
#                 convergence += abs(self.abundance_dict[ref_name] - abundance_dict[ref_name])

#         if self.em_round == 1:
#             self.convergence = 1
#         else:
#             self.convergence = convergence

#         return abundance_dict

#     def _update_compatibility(self):
#         """
#         Update read-transcript compatibility based on transcript abundances
#         """
#         compatibility_dict = defaultdict(dict)

#         for read_name, comp in self.compatibility_dict.items():
#             total = 0
#             for ref_name in comp.keys():
#                 total += self.abundance_dict[ref_name]

#             for ref_name in comp.keys():
#                 compatibility_dict[read_name][ref_name] = self.abundance_dict[ref_name] / total

#         return compatibility_dict

#     def _get_tx_df(self):
#         """
#         Extract transcript info from bam file header
#         """
#         try:
#             with pysam.AlignmentFile(self.alignment_file) as bam:
#                 references = bam.references
#                 lengths = bam.lengths

#             return pd.DataFrame(index=references, data=lengths, columns=["transcript_length"])
#         # If any error return empty DataFrame silently
#         except Exception:
#             return pd.DataFrame()
