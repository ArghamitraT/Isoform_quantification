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
import numpy as np
from scipy.special import psi, gammaln
import time
import pickle
from scipy.stats import spearmanr, pearsonr
from contextlib import redirect_stdout

from DirichletOptimizer import DirichletModel 
#from dirichlet_pyro import DirichletModel


# Local imports
from NanoCount.Read import Read
from NanoCount.common import *

crnt_tm = datetime.datetime.now()


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
            max_em_rounds: int = 10, 
            extra_tx_info: bool = False,
            primary_score: str = "alignment_score",
            max_dist_3_prime: int = 50,
            max_dist_5_prime: int = -1,
            verbose: bool = False,
            quiet: bool = False,
            alpha_initial: float = 1, 
            GD_lr: float = 0.01,
            process: str ='expectation_log_theta', #(AT) 'expectation_log_theta' or 'log_expectation_theta'
            load: int = 0,
            load_filename: str = ""
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

        print("Checking options and input files")
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
        self.expectation_log_theta = {}
        self.elbo = {}
        self.alpha_initial = alpha_initial
        self.GD_lr = GD_lr
        self.process=process
        self.load=load
        self.load_filename = load_filename
        new_max_em_rounds = max_em_rounds

        print("Initialise Nanocount")
        print("Parse Bam file and filter low quality alignments")
        
        start = time.time()

        start1 = time.time() # comment

        if load:
            self.load_state(max_em_rounds, count_file)
        
        else:
            self.initialize_model()
        
        # comment
        end1 = time.time()
        interval = (end1-start1)/60
        print(f"time_initialize {interval} min")

        start1 = time.time() # comment

         # # Initialize the Dirichlet optimizer with the theta data
        dirichlet_optimizer = DirichletModel(self.all_alpha, self.GD_lr, self.process)
        
        # comment
        end1 = time.time()
        interval = (end1-start1)/60
        print(f"time_initialize_diri {interval} min")
        
        
        start1 = time.time() # comment

        # filename for saving variables
        token = (self.count_file.split('/')[-1]).split('_')[-1]
        model_save_path = '/'.join(self.count_file.split('/')[:-3])+'/weights/'
        saved_state_filename = self.create_saved_state_filename(result="allWeights")
        final_save_path = self.create_image_name(model_save_path+saved_state_filename+'_token_'+token, format=".pkl")
        saved_EMstat_filename = self.create_saved_state_filename(result="allEMstats")
        final_EMstat_path = self.create_image_name(model_save_path+saved_EMstat_filename+'_token_'+token, format=".csv")
        saved_GDloss_filename = self.create_saved_state_filename(result="allGDloss")
        final_GDloss_path = self.create_image_name(model_save_path+saved_GDloss_filename+'_token_'+token, format=".pkl")

        # Initialize an empty DataFrame to store the stats
        stats_df = pd.DataFrame(columns=[
            'EM_loop', 'Alpha_summation', 'EM_convergence', 
            'Spearman_corr_theta1_theta2', 
            'Spearman_corr_theta1_alpha', 
            'Spearman_corr_theta2_alpha',
            'pearson_corr_theta1_theta2',
            'pearson_corr_theta1_alpha', 
            'pearson_corr_theta2_alpha'
        ] + [f'ELBO_sample_{i+1}' for i in range(len(self.all_read_dicts))] 
        + [f'Convergence_sample_{i+1}' for i in range(len(self.all_read_dicts))])
        # Initialize an empty DataFrame to store the GD loss history
        GDloss_history_df = pd.DataFrame(columns=['GD_Loss'])

        # comment
        end1 = time.time()
        interval = (end1-start1)/60
        print(f"time_household_before_EM {interval} min")
        
        """ EM and Gradient descent """
        # Iterate until convergence threshold or max EM round are reached
        #while self.convergence > self.convergence_target and self.em_round < self.max_em_rounds: ## (AT)
        while self.em_round < self.max_em_rounds:
            self.convergence = 0
            self.em_round += 1
            elbo_values = {}
            convergence_values = {}

            start1 = time.time() # comment

            # EM
            sample_num = 0
            for sample_key in self.all_Yri:
                self.all_alpha_prime[sample_key] = self.update_alpha_prime(sample_key)
                self.expectation_log_theta[sample_key] = self.calculate_expectation_log_theta_m(sample_key)
                self.all_Phi_ri[sample_key] = self.update_Phi_ri(sample_key)
                self.all_theta[sample_key], convergence = self.update_expectation_theta(sample_key)
                self.elbo[sample_key] = self.calculate_elbo(sample_key)
                sample_num +=1
                self.convergence += convergence

                # Store ELBO and Convergence for this sample
                elbo_values[f'ELBO_sample_{sample_num}'] = self.elbo[sample_key]
                convergence_values[f'Convergence_sample_{sample_num}'] = convergence

                print(f"ELBO_sample_{sample_num} {self.elbo[sample_key]}")
                print(f"Convergence_sample_{sample_num} {convergence}")
            
            # comment
            end1 = time.time()
            interval = (end1-start1)/60
            print(f"time_EM_calculation {interval} min")

            start1 = time.time() # comment

            # UPDATE ALPHA (AT)
            # self.all_alpha, GDloss_history = dirichlet_optimizer.update_alpha(expectation_log_theta=self.expectation_log_theta, 
            #                                                   all_theta=self.all_theta)
            
            # comment
            end1 = time.time()
            interval = (end1-start1)/60
            print(f"time_GDloss {interval} min")


            print("alpha_summation ", np.sum(list(self.all_alpha.values())))
            self.convergence = self.convergence/sample_num
            print(f"EM_loop {self.em_round}" )
            print(f"EM_convergence {self.convergence}")

            start1 = time.time() # comment

            # generate model fit correlation (AT)
            # spearman_corr_theta1_theta2, pearson_corr_theta1_theta2 = self.spearman_pearson_corr('theta1_theta2')
            # spearman_corr_theta1_alpha, pearson_corr_theta1_alpha = self.spearman_pearson_corr('theta1_alpha')
            # spearman_corr_theta2_alpha, pearson_corr_theta2_alpha = self.spearman_pearson_corr('theta2_alpha')

            spearman_corr_theta1_theta2, pearson_corr_theta1_theta2 = 0, 0
            spearman_corr_theta1_alpha, pearson_corr_theta1_alpha = self.spearman_pearson_corr('theta1_alpha')
            spearman_corr_theta2_alpha, pearson_corr_theta2_alpha = 0, 0


            # comment
            end1 = time.time()
            interval = (end1-start1)/60
            print(f"time_corr_calculation {interval} min")

            start1 = time.time() # comment

            # Create a DataFrame for the current iteration
            current_stats = pd.DataFrame([{
                'EM_loop': self.em_round,
                'Alpha_summation': np.sum(list(self.all_alpha.values())),
                'EM_convergence': self.convergence,
                'Spearman_corr_theta1_theta2': spearman_corr_theta1_theta2,
                'Spearman_corr_theta1_alpha': spearman_corr_theta1_alpha,
                'Spearman_corr_theta2_alpha': spearman_corr_theta2_alpha,
                'pearson_corr_theta1_theta2': pearson_corr_theta1_theta2,
                'pearson_corr_theta1_alpha': pearson_corr_theta1_alpha, 
                'pearson_corr_theta2_alpha': pearson_corr_theta2_alpha,
                **elbo_values,
                **convergence_values
            }])

            # Concatenate the current iteration's DataFrame with the main DataFrame
            stats_df = pd.concat([stats_df, current_stats], ignore_index=True)

            # (AT)
            # Create a DataFrame for the loss history 
            # GDloss_history_current = pd.DataFrame({
            #     'GD_Loss': GDloss_history
            # })
            # # Append the loss history of the current iteration to the loss history DataFrame
            # GDloss_history_df = pd.concat([GDloss_history_df, GDloss_history_current], ignore_index=True)

            # Save the state after each iteration
            with open(final_GDloss_path, 'wb') as f:
                pickle.dump(GDloss_history_df, f)
            stats_df.to_csv(final_EMstat_path, index=False)
            self.save_state(final_save_path)

            # comment
            end1 = time.time()
            interval = (end1-start1)/60
            print(f"time_rest {interval} min")

        
        end = time.time()
        interval = (end-start)/60
        print(f"time_EM {interval} min")

        print(f"Exit EM loop after {self.em_round} rounds")
        print(f"Convergence value: {self.convergence}")
        if not self.convergence <= self.convergence_target:
            print(f"Convergence target {self.convergence_target} could not be reached after {self.max_em_rounds} rounds")
    
        # Write out results
        print("Summarize data")
        indx = 0
        for sample, theta in self.all_theta.items():
            print(f"Processing {sample}")
            count_df = pd.DataFrame(theta.most_common(), columns=["transcript_name", "raw"])
            count_df.set_index("transcript_name", inplace=True, drop=True)

            print("Compute estimated counts and TPM")
            # Adjusted to use the length of read_dict for the current sample
            count_df["est_count"] = count_df["raw"] * len(self.all_Yri[sample])
            count_df["tpm"] = count_df["raw"] * 1000000 / sum(theta.values())

            # Add extra transcript info if required
            if self.extra_tx_info:
                tx_df = self._get_tx_df()
                count_df = pd.merge(count_df, tx_df, left_index=True, right_index=True, how="outer")

            # Cleanup and sort
            count_df.sort_values(by="raw", ascending=False, inplace=True)
            count_df.fillna(value=0, inplace=True)
            count_df.index.name = "transcript_name"

            file_name = self.create_saved_state_filename(result=self.count_file + '_' + sample)
            file_name_timestamp = self.create_image_name(file_name, format="")
            count_file = f"{file_name_timestamp}.tsv" if self.count_file else None
            if count_file:
                print(f"Write file for {sample}")
                count_df.to_csv(count_file, sep="\t")
            indx += 1
        # end = time.time()
        # print(f"Time taken to run the code was {end - start} seconds")


     # ~~~~~~~~~~~~~~NEW FUNCTIONS~~~~~~~~~~~~~~ #
    def create_saved_state_filename(self, result):
        
        # Iterate through the list of file names
        for index, file_path in enumerate(self.file_names_list, start=1):
            # Split the file path by '/' and take the last part (the file name)
            file_name = file_path.split('/')[-1]
            # Extract a specific part of the file name if necessary (e.g., removing extension)
            file_identifier = ''.join(file_name.split('_')).split('.')[0]
            # Construct the string
            result += f"_file{index}_{file_identifier}"
        result = f"{result}_GDlr_{self.GD_lr}_AlphaInitial_{self.alpha_initial}_EMround_{self.max_em_rounds}"
        return result


     # Add these methods to your class
    def save_state(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    # @staticmethod
    def load_state(self,  new_max_em_rounds, count_file):
        # with open(self.load_filename, 'rb') as f:
        #     return pickle.load(f)
        print(f"loading_weights_from {self.load_filename}")
        with open(self.load_filename, 'rb') as f:
            loaded_obj = pickle.load(f)
            # Update the current instance's __dict__ with the loaded object's __dict__
            self.__dict__.update(loaded_obj.__dict__)
        self.max_em_rounds = new_max_em_rounds
        self.count_file = count_file
        

    def create_image_name(self, name, format=".png"):
        image_name = (name+"_" + str(crnt_tm.year) + "_" + str(crnt_tm.month) + "_" + str(crnt_tm.day) + "_"
                    + time.strftime("%H_%M_%S") + format)
        return image_name

    def fraction_to_float(self, fraction):
        # Assuming fraction is in string format '1/2', '3/4', etc.
        try:
            numerator, denominator = fraction.split('/')
            return float(numerator) / float(denominator)
        except ValueError:
            # If fraction is not in expected format, return it as a float
            return float(fraction)

    def spearman_pearson_corr(self, comparison_criteria):
        
        if comparison_criteria == 'theta1_theta2':
            model_weight1 = self.all_theta['sample1']
            model_weight2 = self.all_theta['sample2']
        elif comparison_criteria == 'theta1_alpha':
            model_weight1 = self.all_theta['sample1']
            model_weight2 = self.all_alpha
        elif comparison_criteria == 'theta2_alpha':
            model_weight1 = self.all_theta['sample2']
            model_weight2 = self.all_alpha
        
        # if we are correlating theta with alpha, alpha needs to be normalized
        if comparison_criteria == 'theta1_alpha' or comparison_criteria == 'theta2_alpha':
            s = sum(model_weight2.values())
            items = model_weight2.items()
            model_weight2 = {k:(v/s) for (k,v) in items}

        # Initialize lists to store matched TPM and molarity values
        matched_tpm = []
        matched_molarity = []

        # Iterate over first variable to find matching molarity values in second variable
        for tname, tpm in model_weight1.items():
            # Find the corresponding molarity value in second variable
            try:
                # If a matching isoform is found in second variable
                molarity_value = model_weight2[tname]
                # Append the tpm and molarity to the respective lists
                matched_tpm.append(np.log(tpm * 1000000  + 1))
                matched_molarity.append(np.log(molarity_value* 1000000+1))
            except KeyError:
                continue

        # Calculate Spearman's correlation using the matched lists
        correlation, p_value = spearmanr(matched_tpm, matched_molarity)
        spearman_correlation = correlation
        # Output the results
        if comparison_criteria == 'theta1_theta2':
            print(f'Spearman_corr_theta1_theta2 {correlation}')
        elif comparison_criteria == 'theta1_alpha':
            print(f'Spearman_corr_theta1_alpha {correlation}')
        elif comparison_criteria == 'theta2_alpha':
            print(f'Spearman_corr_theta2_alpha {correlation}')
        
        # Calculate Pearsons's correlation using the matched lists
        correlation, p_value = pearsonr(matched_tpm, matched_molarity)
        pearson_correlation = correlation
        # Output the results
        if comparison_criteria == 'theta1_theta2':
            print(f'Pearson_corr_theta1_theta2 {correlation}')
        elif comparison_criteria == 'theta1_alpha':
            print(f'Pearson_corr_theta1_alpha {correlation}')
        elif comparison_criteria == 'theta2_alpha':
            print(f'Pearson_corr_theta2_alpha {correlation}')
        
        return spearman_correlation, pearson_correlation

    
    def calculate_elbo(self, sample_key):
        """
        Calculate the Evidence Lower BOund (ELBO) for a given sample.

        Parameters:
        sample_key (str): The key for the sample in the dataset.

        Returns:
        float: The ELBO value.
        """
        # Extract necessary variables
        Phi_nm = self.all_Phi_ri[sample_key]
        Pnm = self.all_read_iso_prob[sample_key]
        alpha_prime = self.all_alpha_prime[sample_key]
        expectation_log_theta = self.expectation_log_theta[sample_key]

        # Initialize ELBO components
        elbo = 0.0

        # Calculate the first component: \sum_{n=1}^{N} \sum_{m=1}^{M} \phi_{nm} \left( \log \frac{p_{nm}}{\phi_{nm}} + \psi(\alpha_m) - \psi\left(\sum_{m=1}^{M} \alpha_m'\right) \right)
        sum_alpha_prime = sum(alpha_prime.values())
        # Define a small value
        epsilon = 1e-10

        for n, phi_n in Phi_nm.items():
            for m, phi_nm in phi_n.items():
                p_nm = Pnm[n][m]
                phi_nm_adjusted = phi_nm if phi_nm != 0 else epsilon # Adjust phi_nm to avoid zero values
                try:
                    elbo += phi_nm * (np.log(p_nm / phi_nm_adjusted) + expectation_log_theta[m])
                except Exception as e:
                    print(f"Exception encountered: {e} for isoform {m} and read {n}")
                if np.isnan(elbo):
                    raise ValueError(f"NaN encountered in first component for isoform {m} and read {n}")

        # Calculate the second component: \log \frac{\Gamma(\sum_{m=1}^{M} \alpha_m)}{\Gamma(\sum_{m=1}^{M} \alpha_m')}
        alpha = {k: self.all_alpha[k] for k in alpha_prime.keys()}
        sum_alpha = sum(alpha.values())
        log_gamma_sum_alpha = gammaln(sum_alpha)
        log_gamma_sum_alpha_prime = gammaln(sum_alpha_prime)
        elbo += log_gamma_sum_alpha - log_gamma_sum_alpha_prime

        # Calculate the third component: \sum_{m=1}^{M} \log \frac{\Gamma(\alpha_m')}{\Gamma(\alpha_m)}
        for isoform in alpha_prime.keys():
            log_gamma_alpha_prime_m = gammaln(alpha_prime[isoform])
            log_gamma_alpha_m = gammaln(alpha[isoform])
            elbo += log_gamma_alpha_prime_m - log_gamma_alpha_m

        # Calculate the fourth component: \sum_{m=1}^{M} (\alpha_m - \alpha_m') \left( \expectation_log_theta[isoform] \right)
        for isoform in alpha_prime.keys():
            elbo += (alpha[isoform] - alpha_prime[isoform]) * expectation_log_theta[isoform]

        return elbo

    def update_alpha_prime(self, sample_key):
        """
        Calculates the initial model parameter (isoform percentage) and alpha prime.
        """
        all_alpha = self.all_alpha
        all_Zri = self.all_Phi_ri[sample_key]
        all_theta = self.all_theta[sample_key]

        # Initialize all_alpha_prime as an empty dictionary
        all_alpha_prime = {isoform: all_alpha[isoform] for isoform in all_theta}

        # Iterate over each read in all_Zri
        for read_key, isoform_dict in all_Zri.items():
            # Iterate over each isoform and its proportion (phi)
            for isoform, phi in isoform_dict.items():
                if isoform in all_alpha_prime:
                    # Update alpha_prime for the isoform
                    all_alpha_prime[isoform] += phi

        return all_alpha_prime

    def update_expectation_theta(self, sample_key):
        """
        Calculate the expected value E_Q(theta)[theta_m] for a Dirichlet distribution.

        Parameters:
        sample_key (str): The key for the sample in the dataset.

        Returns:
        Counter: The expected theta values as a Counter.
        """
        convergence = 0
        theta_old = self.all_theta[sample_key]
        alpha_prime = self.all_alpha_prime[sample_key]
        sum_alpha_prime = sum(alpha_prime.values())
        expected_theta = Counter()

        for isoform, alpha_value in alpha_prime.items():
            expected_theta[isoform] = alpha_value / sum_alpha_prime

        if self.em_round > 0:
            for ref_name in theta_old.keys():
                convergence += abs(theta_old[ref_name] - expected_theta[ref_name])

        return expected_theta, convergence

    def calculate_expectation_log_theta_m(self, sample_key):
        """
        Calculate the expected log value E_Q(theta)[log(theta_m)] for a Dirichlet distribution.

        Parameters:
        sample_key (str): The key for the sample in the dataset.

        Returns:
        Counter: The expected log theta values as a Counter.
        """
        alpha_prime = self.all_alpha_prime[sample_key]
        sum_alpha_prime = np.sum(list(alpha_prime.values()))
        eq_log_theta_m = Counter()

        for isoform, alpha_value in alpha_prime.items():
            eq_log_theta_m[isoform] = psi(alpha_value) - psi(sum_alpha_prime)

        return eq_log_theta_m
    

    def update_Phi_ri(self, sample_key):

        """
           Update phi_nm using the equation:
           phi_nm ∝ p_nm * exp(E_Q(theta)[log(theta_m)])
           """
        updated_phi = {}
        expectation_log_theta_m = self.expectation_log_theta[sample_key]

        for read, isoform_probs in self.all_read_iso_prob[sample_key].items():
            updated_phi[read] = {}
            for isoform, p_nm in isoform_probs.items():
                updated_phi[read][isoform] = p_nm * np.exp(expectation_log_theta_m[isoform])

            # Normalize phi values for the current read
            total_phi = np.sum(list(updated_phi[read].values()))
            for isoform in updated_phi[read]:
                updated_phi[read][isoform] /= total_phi

        return updated_phi


    def assign_alpha(self):
        # Initialize a set to store all unique isoforms
        isoforms_set = set()
        
        # Iterate over all samples and add isoforms to the set
        for sample in self.all_theta.values():
            isoforms_set.update(sample.keys())
        
        # Number of unique isoforms
        n = len(isoforms_set)

        # Set the fixed seed for reproducibility
        np.random.seed(5)
        
        # Step 1: Initialize alpha values in an unconstrained space
        alpha_raw = np.random.randn(n)
        
        # Step 2: Apply the softmax function to ensure they sum to 1
        alpha_softmax = np.exp(alpha_raw) / np.sum(np.exp(alpha_raw))
        
        # Step 3: Scale the resulting values to have the desired fixed sum
        alpha_scaled = alpha_softmax * self.alpha_initial
        
        # Initialize a dictionary to store alpha values for each isoform
        alpha = {}
        
        # Assign the scaled alpha values to each isoform
        for isoform, value in zip(isoforms_set, alpha_scaled):
            alpha[isoform] = value  # Convert tensor to float
        
        return alpha
    

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

        # return abundance_dict, all_alpha, isoform_indices
        return abundance_dict

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


    def calculate_Z(self, old_compatibility_dict, theta, Pnm):
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
            total = 0
            Z_ri[read] = {}
            for isoform, Y_ri in isoform_values.items():
                p_nm = Pnm[read][isoform]
                theta_i = theta[isoform]  # Get theta for isoform
                temp_val = p_nm * theta_i / denominator
                Z_ri[read][isoform] = temp_val  # Calculate Z_ri according to the formula
                total +=temp_val
            # Normalize phi values for the current read
            for isoform in Z_ri[read]:
                Z_ri[read][isoform] /= total
        
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

    
    def hash_key(self, elements):
        # Create a stable hash of the sorted elements using a cryptographic hash function
        # Convert the elements tuple to a string and encode to bytes before hashing
        return hashlib.sha256(str(sorted(elements)).encode('utf-8')).hexdigest()

    
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

        read_len_dispro = 0
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
                    
                    if np.isnan(p_nm) or p_nm < 0:
                        read_isoform_prob[read_name][string.rname] = 1
                        # print(f"read len {ref_len_dict[string.rname]}, isoform len {string.align_len}") #(AT)
                        read_len_dispro+=1

                    else:
                        read_isoform_prob[read_name][string.rname] = p_nm
        print(f"total_read_len_longer_than_reflen {read_len_dispro}")
        return compatibility_dict, read_isoform_prob


    def initialize_model(self):

        """ EXPLANATION
            * Yri --> Binary Compatibility Matrix
            * theta --> isoform percentage/abundance/quantity
            * Phi_ri --> Expectation of A (true asignment) matrix
            * n --> # of reads for each isoform
        """

        #(AT)
        # simulation_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/round4_small/' 
        simulation_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/round8/'
        for index, file_name in enumerate(self.file_names_list, start=1):
            sample_key = f'sample{index}'
            file_name = file_name.split('/')[-1]
            
            dir = os.path.join(simulation_dir, f"{file_name}_read_dict.pkl")
            with open(dir, 'rb') as file:
                self.all_Yri[sample_key] = pickle.load(file)
                print(f"loading {dir}")
            dir = os.path.join(simulation_dir, f"{file_name}_read_iso_prob.pkl")
            with open(dir, 'rb') as file:
                self.all_read_iso_prob[sample_key] = pickle.load(file)
                print(f"loading {dir}")

        # All the initial calculation
        for sample_key in self.all_Yri:
            #self.all_Yri[sample_key], self.all_read_iso_prob[sample_key] = self.get_compatibility_modified(sample_key)
            self.all_theta[sample_key] \
                = self.calculate_theta_and_alpha_prime_0(self.all_Yri[sample_key])
            demo_phi = self.calculate_Z(self.all_Yri[sample_key], self.all_theta[sample_key], self.all_read_iso_prob[sample_key])
            self.all_Phi_ri[sample_key] = demo_phi

        self.all_alpha = self.assign_alpha()

        # COMMENT
        print("Initiation")

        # find out the reads mathced to more than 1 isoform
        # more_than_one = {key: val for key, val in compatibility_dict_long.items() if len(val) > 1}

        # EM loop to calculate abundance and update read-transcript compatibility
        print("Start EM abundance estimate")

        self.em_round = 0
        self.convergence = 1

        # Initialize the Dirichlet optimizer with the theta data
        #self.dirichlet_optimizer = DirichletModel(self.all_alpha, self.GD_lr, self.process)

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