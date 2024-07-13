"""
This file implements EM, VI and calls GD.
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
import generate_images as gen_img
from scipy.special import psi, gammaln
import time
import pickle

# Local imports
from NanoCount.Read import Read
from NanoCount.common import *


# ~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~ #
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
            convergence_target: float = 0.001,
            max_em_rounds: int = 100,
            extra_tx_info: bool = False,
            primary_score: str = "alignment_score",
            max_dist_3_prime: int = 50,
            max_dist_5_prime: int = -1,
            verbose: bool = False,
            quiet: bool = False,
            alpha_initial: float = 3,
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
        self.all_Phi_ri = {}
        self.all_n = {}
        self.all_alpha = {}
        self.all_alpha_prime = {}
        self.all_isoform_indices = {}
        self.all_read_iso_prob = {}
        self.exp_log_theta = {}
        self.elbo = {}
        self.alpha_initial = alpha_initial

        self.log.warning("Initialise Nanocount")

        # Collect all alignments grouped by read name
        self.log.info("Parse Bam file and filter low quality alignments")
        
        start = time.time()

        #  (AT) UNCOMMENT
        # Loop over all file names provided and parses the reads with
        for index, file_name in enumerate(self.file_names_list, start=1):
            # Parse the BAM file
            read_dict, ref_len_dict = self._parse_bam(file_name=file_name)
            # Store the dictionaries in the universal dictionary with a sample key
            sample_key = f'sample{index}'
            self.all_read_dicts[sample_key] = read_dict
            self.all_ref_len_dicts[sample_key] = ref_len_dict
        
        end = time.time()
        interval = (end-start)/60
        print(f"time: time to parse the reads {interval} min")

        # Set the environment variables
        # os.environ['PYDEVD_WARN_EVALUATION_TIMEOUT'] = '10'  # Increase timeout to 10 seconds
        # os.environ['PYDEVD_THREAD_DUMP_ON_WARN_EVALUATION_TIMEOUT'] = 'true'
        # os.environ['PYDEVD_UNBLOCK_THREADS_TIMEOUT'] = '10'  # Unblock threads after 10 seconds

        # with open('pkl_files/Pac_illu_set1_read_dicts', 'wb') as file:
        #     pickle.dump(self.all_read_dicts, file)
        # with open('pkl_files/Pac_illu_set1_ref_len_dicts', 'wb') as file:
        #     pickle.dump(self.all_ref_len_dicts, file)

        # print("Data has been saved.")



        # (AT) COMMENT
        # because parsing takes a bit of time, saved couple of file to develop working code
        # with open('read_dict_long.pkl', 'rb') as file:
        #     self.all_read_dicts['sample1'] = pickle.load(file)
        # with open('ref_len_dict_long.pkl', 'rb') as file:
        #     self.all_ref_len_dicts['sample1'] = pickle.load(file)
        # with open('pkl_files/read_nano_ambi.pickle', 'rb') as file:
        #     self.all_read_dicts['sample1'] = pickle.load(file)
        # with open('pkl_files/ref_len_nano_ambi.pickle', 'rb') as file:
        #     self.all_ref_len_dicts['sample1'] = pickle.load(file)
        # with open('pkl_files/read_dict_short.pkl', 'rb') as file:
        #     self.all_read_dicts['sample2'] = pickle.load(file)
        # with open('pkl_files/ref_len_dict_short.pkl', 'rb') as file:
        #     self.all_ref_len_dicts['sample2'] = pickle.load(file)

        if self.filter_bam_out:
            self.log.info("Write selected alignments to BAM file")
            self._write_bam()

        # Generate compatibility dict grouped by reads
        # self.log.info("Generate initial read/transcript compatibility index")

        """ EXPLANATION
            * Yri --> Binary Compatibility Matrix
            * theta --> isoform percentage/abundance/quantity
            * Phi_ri --> Expectation of A (true asignment) matrix
            * n --> # of reads for each isoform
        """
        # All the initial calculation
        for sample_key in self.all_read_dicts:
            self.all_Phi_ri[sample_key], self.all_read_iso_prob[sample_key] = self.get_compatibility_modified(sample_key)
            self.all_theta[sample_key], self.all_alpha[sample_key], self.all_isoform_indices[sample_key] \
                = self.calculate_theta_and_alpha_prime_0(self.all_Phi_ri[sample_key])

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

        # with tqdm(
        #         unit=" rounds",
        #         unit_scale=True,
        #         desc="\tProgress",
        #         disable=(quiet or verbose),
        # ) as pbar:

        # Initialize the Dirichlet optimizer with the theta data
        dirichlet_optimizer = DirichletModel(self.all_theta, self.all_isoform_indices, self.all_alpha, self.exp_log_theta)

        """ EM and Gradient descent """
        # Iterate until convergence threshold or max EM round are reached
        while self.convergence > self.convergence_target and self.em_round < self.max_em_rounds:
            self.convergence = 0
            self.em_round += 1

            # EM
            sample_num = 0
            for sample_key in self.all_read_dicts:
                self.all_alpha_prime[sample_key] = self.update_alpha_prime(sample_key)
                self.exp_log_theta[sample_key] = self.calculate_exp_log_theta_m(sample_key)
                self.all_Phi_ri[sample_key] = self.update_Phi_ri(sample_key)
                self.all_theta[sample_key], convergence = self.update_exp_theta(sample_key)
                self.elbo[sample_key] = self.calculate_elbo(sample_key)
                sample_num +=1
                self.convergence += convergence
                print(f"ELBO_sample_{sample_num} {self.elbo[sample_key]}")
                print(f"Convergence_sample_{sample_num} {convergence}")


            # UPDATE ALPHA
            dirichlet_optimizer.update_alpha()
            self.convergence = self.convergence/len(self.all_alpha)
            self.log.info("EM_loop {}".format(self.em_round))
            self.log.info("EM_convergence {} ".format(self.convergence))

            # store the values for downstream
            # alpha_history.append(np.mean(self.alpha))
            convergence_history.append(self.convergence)
            # pbar.update(1)
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
            file_name = (self.count_file + '_' + sample + '_' +
                         self.file_names_list[indx].split('/')[-1].split('.')[0])
            file_name_timestamp = gen_img.create_image_name(file_name, format="")
            count_file = f"{file_name_timestamp}.tsv" if self.count_file else None
            if count_file:
                self.log.info(f"Write file for {sample}")
                count_df.to_csv(count_file, sep="\t")
            indx += 1
        # end = time.time()
        # print(f"Time taken to run the code was {end - start} seconds")


     # ~~~~~~~~~~~~~~NEW FUNCTIONS (AT)~~~~~~~~~~~~~~ #
    def calculate_elbo(self, sample_key):

        # Extract necessary variables
        Phi_nm = self.all_Phi_ri[sample_key]
        Pnm = self.all_read_iso_prob[sample_key]
        alpha = self.all_alpha[sample_key]
        alpha_prime = self.all_alpha_prime[sample_key]
        exp_log_theta = self.exp_log_theta[sample_key]
        isoform_indices = self.all_isoform_indices[sample_key]

        # Initialize ELBO components
        elbo = 0.0

        # Calculate the first component: \sum_{n=1}^{N} \sum_{m=1}^{M} \phi_{nm} \left( \log \frac{p_{nm}}{\phi_{nm}} + \psi(\alpha_m) - \psi\left(\sum_{m=1}^{M} \alpha_m'\right) \right)
        sum_alpha_prime = sum(alpha_prime)
        # psi_sum_alpha_prime = psi(sum_alpha_prime)
        for n, phi_n in Phi_nm.items():
            for m, phi_nm in phi_n.items():
                p_nm = Pnm[n][m]
                elbo += phi_nm * (np.log(p_nm / phi_nm) + exp_log_theta[m])
                if np.isnan(elbo):
                    raise ValueError(f"NaN encountered in first component for isoform {m} and read {n}")


        # Calculate the second component: \log \frac{\Gamma(\sum_{m=1}^{M} \alpha_m)}{\Gamma(\sum_{m=1}^{M} \alpha_m')}
        sum_alpha = sum(alpha)
        log_gamma_sum_alpha = gammaln(sum_alpha)
        log_gamma_sum_alpha_prime = gammaln(sum_alpha_prime)
        elbo += log_gamma_sum_alpha - log_gamma_sum_alpha_prime

        # Calculate the third component: \sum_{m=1}^{M} \log \frac{\Gamma(\alpha_m')}{\Gamma(\alpha_m)}
        for isoform, index in isoform_indices.items():
            log_gamma_alpha_prime_m = gammaln(alpha_prime[index])
            log_gamma_alpha_m = gammaln(alpha[index])
            elbo += log_gamma_alpha_prime_m - log_gamma_alpha_m

        # Calculate the fourth component: \sum_{m=1}^{M} (\alpha_m - \alpha_m') \left( \exp_log_theta[isoform] \right)
        for isoform, index in isoform_indices.items():
            elbo += (alpha[index] - alpha_prime[index]) * exp_log_theta[isoform]

        return elbo

    def update_alpha_prime(self, sample_key):
        """
        Calculates the initial model parameter (isoform percentage) and alpha prime.
        """
        all_alpha = self.all_alpha[sample_key]
        all_Zri = self.all_Phi_ri[sample_key]
        isoform_indices = self.all_isoform_indices[sample_key]

        # Initialize all_alpha_prime
        all_alpha_prime = all_alpha

        # Iterate over each read in all_Zri
        for read_key, isoform_dict in all_Zri.items():
            # Iterate over each isoform and its proportion (phi)
            for isoform, phi in isoform_dict.items():
                index = isoform_indices[isoform]
                # Update alpha_prime for the isoform
                all_alpha_prime[index] += phi

        return all_alpha_prime

    def update_exp_theta(self, sample_key):
        """
        Calculate the expected value E_Q(theta)[theta_m] for a Dirichlet distribution.

        Parameters:
        alpha_prime (list or array): The alpha prime values for the isoforms.
        isoform_indices (dict): Mapping of isoform names to their indices.

        Returns:
        Counter: The expected theta values as a Counter.
        """
        convergence = 0
        theta_old = self.all_theta[sample_key]
        alpha_prime = self.all_alpha_prime[sample_key]
        isoform_indices = self.all_isoform_indices[sample_key]
        sum_alpha_prime = sum(alpha_prime)
        expected_theta = Counter()

        for isoform, index in isoform_indices.items():
            expected_theta[isoform] = alpha_prime[index] / sum_alpha_prime

        if self.em_round > 0:
            for ref_name in theta_old.keys():
                convergence += abs(theta_old[ref_name] - expected_theta[ref_name])


        return expected_theta, convergence

    def calculate_exp_log_theta_m(self, sample_key):

        """
        Calculate the expected log value E_Q(theta)[log(theta_m)] for a Dirichlet distribution.

        Parameters:
        alpha_prime (Counter): The alpha prime values for the isoforms.
        isoform_indices (dict): Mapping of isoform names to their indices.

        Returns:
        Counter: The expected log theta values as a Counter.
        """
        alpha_prime = self.all_alpha_prime[sample_key]
        isoform_indices = self.all_isoform_indices[sample_key]
        sum_alpha_prime = np.sum(alpha_prime)
        eq_log_theta_m = Counter()

        for isoform, index in isoform_indices.items():
            eq_log_theta_m[isoform] = psi(alpha_prime[index]) - psi(sum_alpha_prime)

        return eq_log_theta_m

    def update_Phi_ri(self, sample_key):

        """
           Update phi_nm using the equation:
           phi_nm ∝ p_nm * exp(E_Q(theta)[log(theta_m)])
           """
        updated_phi = {}
        exp_log_theta_m = self.exp_log_theta[sample_key]

        for read, isoform_probs in self.all_read_iso_prob[sample_key].items():
            updated_phi[read] = {}
            for isoform, p_nm in isoform_probs.items():
                updated_phi[read][isoform] = p_nm * np.exp(exp_log_theta_m[isoform])

            # Normalize phi values for the current read
            total_phi = np.sum(list(updated_phi[read].values()))
            for isoform in updated_phi[read]:
                updated_phi[read][isoform] /= total_phi

        return updated_phi
    def calculate_theta_and_alpha_prime_0(self, all_Zri):
        """
        Calculates the initial model parameter (isoform percentage) and alpha prime.
        """
        # Initialize the abundance dictionary and total
        abundance_dict = Counter()
        total = 0
        isoform_indices = {}
        isoform_counter = 0

        # Calculate abundance and create isoform indices
        for read_name, comp in all_Zri.items():
            for ref_name, score in comp.items():
                abundance_dict[ref_name] += score
                total += score
                if ref_name not in isoform_indices:
                    isoform_indices[ref_name] = isoform_counter
                    isoform_counter += 1

        # Normalize the abundance dictionary
        for ref_name in abundance_dict.keys():
            abundance_dict[ref_name] = abundance_dict[ref_name] / total

        all_alpha = [self.alpha_initial] * len(abundance_dict)

        return abundance_dict, all_alpha, isoform_indices

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
        numerator = n_i_sample + alpha_i_sample - 1
        denominator = np.sum(numerator) - len(n_i_sample)
        theta_hat_sample = numerator / denominator

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
        isoform_indices = {}
        isoform_counter = 0

        for read_name, comp in compatibility_dict.items():
            for ref_name, score in comp.items():
                abundance_dict[ref_name] += score
                #compatibility_dict[ref_name] += score
                total += score
                isoform_indices[ref_name] = isoform_counter
                isoform_counter += 1

        for ref_name in abundance_dict.keys():
            abundance_dict[ref_name] = abundance_dict[ref_name] / total

        return abundance_dict, isoform_indices

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
        read_isoform_prob = defaultdict(dict)
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
                    compatibility_dict[read_name][string.rname] = 1 / len(alignment.alignment_list)

                    # Check for potential division by zero
                    denominator = ref_len_dict[string.rname] - string.align_len + 1
                    if denominator == 0:
                        p_nm = np.nan  # This will trigger the NaN check below
                    else:
                        p_nm = 1 / denominator

                    # Check if p_nm is NaN or less than 1
                    if np.isnan(p_nm) or p_nm < 1:
                        read_isoform_prob[read_name][string.rname] = 1
                        #print(f"read len {ref_len_dict[string.rname]}, isoform len {string.align_len}")
                    else:
                        read_isoform_prob[read_name][string.rname] = p_nm
        return compatibility_dict, read_isoform_prob



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
                # elif not self.keep_suplementary and alignment.is_supplementary:
                #     c["Discarded supplementary alignments"] += 1
                # elif self.min_alignment_length > 0 and alignment.query_alignment_length < self.min_alignment_length:
                #     c["Discarded short alignments"] += 1
                # elif self.max_dist_3_prime >= 0 and alignment.reference_end <= ref_len_dict[
                #     alignment.reference_name] - self.max_dist_3_prime:
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