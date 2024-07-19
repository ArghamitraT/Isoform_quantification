"""
NOTE 2: this file has some project specific parsing changes. might not need in the long run
NOTE 1:This file implements EM, MAP and calls GD.
The math and detail explanation is on this file: https://drive.google.com/file/d/1LGLhGvn3KRAYunf995ZVAYA4w2lgRYfr/view?usp=sharing
"""


# ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~ #

# Third party imports
import pysam
import pandas as pd
from tqdm import tqdm
import hashlib
import pickle
from DirichletOptimizer import DirichletModel
import numpy as np
import generate_corr as gen_img
import time


# Local imports
from NanoCount.Read import Read
from NanoCount.common import *

# import types
#
# # Assuming Read is a class or function you have defined somewhere
# def Read():
#     pass
#
# # Manually adding Read to __main__
# import sys
# sys.modules['__main__'].Read = Read
#
# with open('pkl_files/PacBio_Day0_rep1.pkl', 'rb') as file:
#     all_read_dicts = pickle.load(file)
#
# print()


# ~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~ #
class Read:
    def __init__(self):
        self.alignment_list = []

    def __repr__(self):
        m = ""
        for r in self.alignment_list:
            m += "\t\t{}\n".format(r)
        return m

    def add_alignment(self, alignment):
        self.alignment_list.append(alignment)

    def to_json(self):
        return self.__repr__()

# Define the Alignment class
class Alignment:
    def __init__(self, q, r, qlen=0, align_len=0, align_score=0, secondary=False):
        self.qname = q
        self.rname = r
        self.qlen = qlen
        self.align_len = align_len
        self.align_score = align_score
        self.secondary = secondary

    def __repr__(self):
        return f"Query:{self.qname} | Reference:{self.rname} | Query len:{self.qlen} | Alignment len:{self.align_len} | Align Score:{self.align_score} | Secondary:{self.secondary}"

class EnhancedOrderedDict(OrderedDict):
    def add(self, key, value=None):
        """Adds a key with a default value if not specified."""
        if key not in self:
            self[key] = value
class Expec_Max:

    # ~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~ #
    def __init__(
            self,
            file_names = [],
            alignment_file: str = "",
            count_file: str = "",
            filter_bam_out: str = "",
            min_alignment_length: int = 50,
            keep_suplementary: bool = False,
            min_query_fraction_aligned: float = 0.5,
            sec_scoring_threshold: float = 0.95,
            sec_scoring_value: str = "alignment_score",
            convergence_target: float = 0.01855,
            #convergence_target: float = 0.5,
            max_em_rounds: int = 10,
            extra_tx_info: bool = False,
            primary_score: str = "alignment_score",
            max_dist_3_prime: int = 50,
            max_dist_5_prime: int = -1,
            verbose: bool = False,
            quiet: bool = False,
            downsampled_files =[],

    ):
        """
        NOTE: File file does not use any hashing and runs EM on both ambiguous and unambiguous reads. This can be improved in future
        EXPLANATION:
        Estimate abundance of transcripts using an EM and MAP
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
        self.file_names_list = file_names
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
        self.all_read_dicts = {}
        self.all_ref_len_dicts = {}
        self.all_Yri = {}
        self.all_theta = {}
        self.all_Zri = {}
        self.all_n = {}


        self.downsampled_files = downsampled_files
        self.temp_read_dicts = defaultdict(Read)
        self.temp_read_lengths = EnhancedOrderedDict()

        self.log.warning("Initialise Nanocount")

        # Collect all alignments grouped by read name
        self.log.info("Parse Bam file and filter low quality alignments")

        #  (AT) UNCOMMENT
        # Loop over all file names provided and parses the reads with
        # for index, file_name in enumerate(self.file_names_list, start=1):
        #     start = time.time()
        #     # Parse the BAM file
        #     read_dict, ref_len_dict = self._parse_bam(file_name=file_name)
        #     # Store the dictionaries in the universal dictionary with a sample key
        #     sample_key = f'sample{index}'
        #     self.all_read_dicts[sample_key] = read_dict
        #     self.all_ref_len_dicts[sample_key] = ref_len_dict
        #     end = time.time()
        #     print(f"Time taken to run the code was {end - start} seconds")

        # (AT) COMMENT
        # because parsing takes a bit of time, saved couple of file to develop working code
        """ SET 1 """
        # with open('pkl_files/Illumina_Day0_rep1_read.pkl', 'rb') as file:
        #     self.all_read_dicts['sample1'] = pickle.load(file)
        # with open('pkl_files/Illumina_Day0_rep1_ref_len.pkl', 'rb') as file:
        #     self.all_ref_len_dicts['sample1'] = pickle.load(file)
        # self.load_data()


        """ SET 2 """
        with open('pkl_files/Illumina_Day0_rep2_read.pkl', 'rb') as file:
            self.all_read_dicts['sample1'] = pickle.load(file)
        with open('pkl_files/Illumina_Day0_rep2_ref_len.pkl', 'rb') as file:
            self.all_ref_len_dicts['sample1'] = pickle.load(file)
        self.load_data()



        # with open('read_dict_long.pkl', 'rb') as file:
        #     self.all_read_dicts['sample1'] = pickle.load(file)
        # with open('ref_len_dict_long.pkl', 'rb') as file:
        #     self.all_ref_len_dicts['sample1'] = pickle.load(file)
        # with open('read_nano_ambi.pickle', 'rb') as file:
        #     self.all_read_dicts['sample1'] = pickle.load(file)
        # with open('ref_len_nano_ambi.pickle', 'rb') as file:
        #     self.all_ref_len_dicts['sample1'] = pickle.load(file)
        # with open('read_dict_short.pkl', 'rb') as file:
        #     self.all_read_dicts['sample2'] = pickle.load(file)
        # with open('ref_len_dict_short.pkl', 'rb') as file:
        #     self.all_ref_len_dicts['sample2'] = pickle.load(file)

        if self.filter_bam_out:
            self.log.info("Write selected alignments to BAM file")
            self._write_bam()

        # Generate compatibility dict grouped by reads
        self.log.info("Generate initial read/transcript compatibility index")

        """ EXPLANATION
            * Yri --> Binary Compatibility Matrix
            * theta --> isoform percentage/abundance/quantity
            * Zri --> Expectation of A (true asignment) matrix
            * n --> # of reads for each isoform
        """
        # All the initial calculation
        for sample_key in self.all_read_dicts:
            self.all_Yri[sample_key] = self.get_compatibility_modified(sample_key)
            self.all_theta[sample_key] = self.calculate_theta_0(self.all_Yri[sample_key])
            self.all_Zri[sample_key] = self.calculate_Z(self.all_Yri[sample_key], self.all_theta[sample_key])
            self.all_n[sample_key] = self.calculate_n(self.all_Zri[sample_key])

        # merging all theta
        # all_theta_keys = self.get_all_keys(self.all_theta)
        # self.update_dicts_with_all_keys(self.all_theta, all_theta_keys)
        # self.update_dicts_with_all_keys(self.all_n, all_theta_keys)

        # find out the reads mathced to more than 1 isoform
        # more_than_one = {key: val for key, val in compatibility_dict_long.items() if len(val) > 1}

        # EM loop to calculate abundance and update read-transcript compatibility
        self.log.warning("Start EM abundance estimate")

        self.em_round = 0
        self.convergence = 1

        # storing the values for each iteration for downstream analysis and plotting
        theta_history = {}  # Dictionary to store theta values for each sample at each iteration
        alpha_history = []  # List to store alpha values at each iteration
        convergence_history = []  # List to store convergence values at each iteration

        start = time.time()
        with tqdm(
                unit=" rounds",
                unit_scale=True,
                desc="\tProgress",
                disable=(quiet or verbose),
        ) as pbar:

            # Initialize the Dirichlet optimizer with the theta data
            dirichlet_optimizer = DirichletModel(all_theta=self.all_theta)
            # Call the method to return alpha and isoform_index to access appropriate alpha to each isoform
            self.alpha, self.isoform_to_index = dirichlet_optimizer.return_alpha()

            # Initialize storage for thetas, alphas, and convergence values
            for sample_key in self.all_read_dicts:
                theta_history[sample_key] = []
                theta_history[sample_key].append(self.all_theta[sample_key])
            alpha_history.append(np.mean(self.alpha))
            convergence_history.append(self.convergence)


            """ EM and Gradient descent """
            # Iterate until convergence threshold or max EM round are reached
            while self.convergence > self.convergence_target and self.em_round < self.max_em_rounds:
            #while self.convergence > 0.005 and self.em_round < 100:

                self.convergence = 0
                self.em_round += 1

                # Gradient Descent
                dirichlet_optimizer.update_alpha()
                # EM
                sample_num = 0
                for sample_key in self.all_read_dicts:
                    self.all_theta[sample_key], convergence = self.update_theta(sample_key)
                    self.all_Zri[sample_key] = self.calculate_Z(self.all_Yri[sample_key], self.all_theta[sample_key])
                    self.all_n[sample_key] = self.calculate_n(self.all_Zri[sample_key])
                    self.convergence +=convergence
                    theta_history[sample_key].append(self.all_theta[sample_key]) # store the values for downstream
                    sample_num +=1

                # self.update_dicts_with_all_keys(self.all_n, all_theta_keys)
                # store the values for downstream
                alpha_history.append(np.mean(self.alpha))
                self.convergence = self.convergence/sample_num
                convergence_history.append(self.convergence)
                print('\nConvergence: ', self.convergence)
                pbar.update(1)
                self.log.debug("EM Round: {} / Convergence value: {}".format(self.em_round, self.convergence))


        self.log.info("Exit EM loop after {} rounds".format(self.em_round))
        self.log.info("Convergence value: {}".format(self.convergence))
        if not self.convergence <= self.convergence_target:
            self.log.error(
                "Convergence target ({}) could not be reached after {} rounds".format(self.convergence_target,
                                                                                      self.max_em_rounds))
        # Plot some figures
        # gen_img.plot_EM_results(alpha_history, convergence_history, theta_history)
        # Write out results
        self.log.warning("Summarize data")

        # Assuming 'all_theta' is a dictionary of Counters, each Counter corresponding to a sample
        # and 'read_dict' is a dictionary of reads for each sample.
        indx = 0
        for sample, theta in self.all_theta.items():
            self.log.info(f"Processing {sample}")
            count_df = pd.DataFrame(theta.most_common(), columns=["transcript_name", "raw"])
            count_df.set_index("transcript_name", inplace=True, drop=True)

            self.log.info("Compute estimated counts and TPM")
            # Adjusted to use the length of read_dict for the current sample
            count_df["est_count"] = count_df["raw"] * len(self.all_read_dicts[sample])
            count_df["tpm"] = count_df["raw"] * 1000000 / sum(theta.values())

            # Add extra transcript info if required
            if self.extra_tx_info:
                tx_df = self._get_tx_df()
                count_df = pd.merge(count_df, tx_df, left_index=True, right_index=True, how="outer")

            # Cleanup and sort
            count_df.sort_values(by="raw", ascending=False, inplace=True)
            count_df.fillna(value=0, inplace=True)
            count_df.index.name = "transcript_name"

            # get Spearman's correlation
            # gen_img.spearman_corr(count_df, sample)

            # The output file could be named according to the sample
            file_name = (self.count_file+'_'+sample+'_'+
                         self.file_names_list[indx].split('/')[-1].split('.')[0])
            file_name_timestamp = gen_img.create_image_name(file_name, format="")
            count_file = f"{file_name_timestamp}.tsv" if self.count_file else None
            if count_file:
                self.log.info(f"Write file for {sample}")
                count_df.to_csv(count_file, sep="\t")
            indx+=1
        end = time.time()
        print(f"Time taken to run the code was {end - start} seconds")

     # ~~~~~~~~~~~~~~NEW FUNCTIONS (AT)~~~~~~~~~~~~~~ #

    def load_data(self):
        with open(self.downsampled_files[0]) as f:
            header = f.readline()
            for line in f:
                line = line.strip('\n')
                q, r = line.split('\t')
                aln = Alignment(q, r)
                self.temp_read_dicts[q].add_alignment(aln)
                self.temp_read_lengths.add((r, 0))
        self.temp_read_lengths = OrderedDict(list(self.temp_read_lengths))
        self.all_read_dicts['sample2'] = self.temp_read_dicts
        self.all_ref_len_dicts['sample2'] = self.temp_read_lengths
        print('Data Loaded and Processed')
        #self.all_read_dicts['sample2']


    # Find the union of all keys across all samples for a given attribute
    def get_all_keys(self,attribute_dict):
        all_keys = set()
        for nested_dict in attribute_dict.values():
            all_keys.update(nested_dict.keys())
        return all_keys

    # Update each sample's dictionary with all keys, setting default value to 0 for new keys
    def update_dicts_with_all_keys(self, attribute_dict, all_keys):
        for sample_key, nested_dict in attribute_dict.items():
            for key in all_keys:
                nested_dict.setdefault(key, 0)

    def update_theta(self, sample_key):
        """ Calculates the mode of the dirichlet """

        """ EXPLANATION:
            * Eqn 5 --> {\hat{\theta}_i  = \frac{n_i + \alpha_i - 1}
            {\sum_{i=1}^{I} (\alpha_i + n_i)-I}}"""

        sample_counts = self.all_n[sample_key]
        convergence = 0
        theta_old = self.all_theta[sample_key]

        # Prepare the arrays for n_i and alpha_i
        n_i_sample = []
        alpha_i_sample = []
        isoform_order = []  # To preserve the order of isoforms

        for isoform in sample_counts:
            index = self.isoform_to_index[isoform]
            n_i_sample.append(sample_counts[isoform])
            alpha_i_sample.append(self.alpha[index])
            isoform_order.append(isoform)

        n_i_sample = np.array(n_i_sample)
        alpha_i_sample = np.array(alpha_i_sample)

        # Calculate theta_hat
        sum_ni_alphai = n_i_sample + alpha_i_sample
        # mode
        numerator = n_i_sample + alpha_i_sample-1
        denominator = np.sum(sum_ni_alphai) - len(n_i_sample)
        # mean
        # numerator = n_i_sample + alpha_i_sample
        # denominator = np.sum(sum_ni_alphai)

        theta_hat_sample = numerator / denominator
        theta_hat_sample[theta_hat_sample < 0] = 0
        print('\ntheta sum: ',sum(theta_hat_sample))

        # Store the results in a Counter with the same isoform order
        theta_hat = Counter({isoform: theta_hat for isoform, theta_hat in zip(isoform_order, theta_hat_sample)})

        if self.em_round > 0:
            for ref_name in theta_old.keys():
                convergence += abs(theta_old[ref_name] - theta_hat[ref_name])

        return theta_hat, convergence

    def calculate_n(self, compatibility_dict):
        """
        Sums up the total assignment for isoform I
        """
        """ EXPLANATION:
            * Eqn 4 -->  n_i = \sum_{r=1}^{R} z_{ri}"""

        abundance_dict = Counter()
        total = 0

        # sums over all the reads and calculates the total abundance of  each isoform
        for read_name, comp in compatibility_dict.items():
            for ref_name, score in comp.items():
                abundance_dict[ref_name] += score
                total += score
        return abundance_dict


    def calculate_Z(self, old_compatibility_dict, theta):
        """
        The EM assignment: Update read-transcript compatibility based on transcript abundances (expectation of A or the ture assignmen)
        """
        """ EXPLANATION:
            * Eqn 3 --> z_{ri}^{t=0} = {y_{ri} \theta_i^{t=0}} /
            {\sum_{i=1}^{I} y_{ri} \theta_i^{t=0}} """

        Z_ri = defaultdict(dict)

        # Loop through each read in Yri
        for read, isoform_values in old_compatibility_dict.items():
            # Calculate the denominator for the Z_ri formula
            denominator = sum(Y_ri * theta[isoform] for isoform, Y_ri in isoform_values.items())
            # Compute Z_ri for each isoform associated with the read
            Z_ri[read] = {}
            for isoform, Y_ri in isoform_values.items():
                theta_i = theta[isoform]  # Get theta for isoform
                Z_ri[read][isoform] = Y_ri * theta_i / denominator  # Calculate Z_ri according to the formula
        return Z_ri


    def calculate_theta_0(self, compatibility_dict):
        """
        Calculates the initial model parameter, or isoform percentage
        """
        """ EXPLANATION:
            * Eqn 1 --> p_{i} = \frac{1}{R} \sum_{r=1}^{R} y_{ri}
            * Eqn 2 --> \theta_{i}^{t=0} = \frac{p_{i}}{\sum_{i=1}^{I} p_{i}} """

        abundance_dict = Counter()
        total = 0

        for read_name, comp in compatibility_dict.items():
            for ref_name, score in comp.items():
                abundance_dict[ref_name] += score
                total += score

        for ref_name in abundance_dict.keys():
            abundance_dict[ref_name] = abundance_dict[ref_name] / total

        return abundance_dict

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


    def get_compatibility_modified(self, sample_key):
        """
         For every read gives the compatible isoforms and normalize them by N-K+1
        """

        """ EXPLANATION:
        * read_dict -> "Binary Compatibility Matrix: Y_{ri}
        * ref_len_dict -> length of the isoforms """
        compatibility_dict = defaultdict(dict)
        read_dict = self.all_read_dicts[sample_key]
        ref_len_dict = self.all_ref_len_dicts[sample_key]

        for read_name, read in read_dict.items():
            # Check if read is iterable and not a string; if not, make it a list
            if not self.is_iterable(read):
                read = [read]  # Wrap non-iterable read in a list

            """ EXPLANATION:
            * if the read length is k and isoform length is n, there are n-k+1 positions the read can come from.
            so we need to multiply the score by (1/n-k+1). This is the way to incorporate isoform length """
            for alignment in read:
                for string in alignment.alignment_list:
                    # (AT) just modified temporarily for PacBio data
                    # try:
                    #     compatibility_dict[read_name][string.rname] = 1 / (ref_len_dict[string.rname]-string.align_len+1)
                    # except:
                    #     compatibility_dict[read_name][string.rname] = 1
                    #     print('problem')
                    try:
                        compatibility_dict[read_name][string.rname.split('|')[0]] = 1 / (ref_len_dict[string.rname]-string.align_len+1)
                    except:
                        compatibility_dict[read_name][string.rname.split('|')[0]] = 1
                        print('problem')
        return compatibility_dict



    # ~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~ #

    def _parse_bam(self, file_name):
        """
        Parse Bam/Sam file, group alignments per reads, filter reads based on
        selection criteria and return a dict of valid read/alignments
        """
        # Parse bam files
        read_dict = defaultdict(Read)
        ref_len_dict = OrderedDict()
        c = Counter()

        aligned_read = file_name

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
                elif self.max_dist_3_prime >= 0 and alignment.reference_end <= ref_len_dict[
                    alignment.reference_name] - self.max_dist_3_prime:
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
            self.log.error(
                "No valid secondary alignments found in bam file. Were the reads aligned with minimap `-p 0 -N 10` options ?")

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
                    compatibility_dict[read_name][string.rname.split('|')[0]] = score = 1.0 / alignment.n_alignment
                    # compatibility_dict[read_name][string.rname] = \
                    #     score = 1.0 / alignment.n_alignment

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
