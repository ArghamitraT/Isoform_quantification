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

from DirichletOptimizer_vector import DirichletModel 
#from dirichlet_pyro import DirichletModel


# Local imports
from NanoCount.Read import Read
from NanoCount.common import *

crnt_tm = datetime.datetime.now()


# ~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~ #
class process_bam:

    # ~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~ #
    def __init__(
            self,
            file_names = [],
            alignment_file: str = "",
            count_file: str = "",
            filter_bam_out: str = "",
            min_alignment_length: int = 50,
            keep_supplementary: bool = False,
            keep_neg_strand: bool = False,
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
            process: str ='theta', #(AT) 'expectation_log_theta' or 'theta'
            dirichlet_builtin: int = 1, 
            load: int = 0,
            load_filename: str = "",
            experiment_num: int = 4,
            parse_original: int = 0,
            process_dir: str = "",
            bam_dir: str = ""
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
        self.keep_supplementary = keep_supplementary
        self.max_dist_5_prime = max_dist_5_prime
        self.max_dist_3_prime = max_dist_3_prime
        self.keep_neg_strand = keep_neg_strand
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
        self.all_rnames = {}
        self.all_readName = {}
        self.theta_names = {}
        self.experiment_num = experiment_num
        self.dirichlet_builtin=dirichlet_builtin
        self.parse_original = parse_original
        self.process_dir = process_dir
        self.bam_dir = bam_dir


        print("Initialise Nanocount")
        print("Parse Bam file and filter low quality alignments")
        
        start = time.time()
        start1 = time.time() # comment

        if load:
            self.load_state(max_em_rounds, count_file)
        
        else:
            self.initialize_model()
        
        
        # start1 = time.time() # comment
        
        # ## Initialize the Dirichlet optimizer with the theta data
        # dirichlet_optimizer = DirichletModel(self.all_alpha, self.GD_lr, self.dirichlet_builtin, self.process)
        
        # # comment
        # end1 = time.time()
        # interval = (end1-start1)/60
        # print(f"time_initialize_diri {interval} min")
        # start1 = time.time() # comment
        
        # # filename for saving variables
        # token = (self.count_file.split('/')[-1]).split('_')[-1]
        # model_save_path = '/'.join(self.count_file.split('/')[:-3])+'/weights/'
        
        # if self.experiment_num == 5:
        #     saved_state_filename = self.create_saved_state_filename_exp5(result="allWeights")
        #     final_save_path = self.create_image_name(model_save_path+saved_state_filename+'_token_'+token, format=".pkl")
        #     saved_EMstat_filename = self.create_saved_state_filename_exp5(result="allEMstats")
        #     final_EMstat_path = self.create_image_name(model_save_path+saved_EMstat_filename+'_token_'+token, format=".csv")
        #     saved_GDloss_filename = self.create_saved_state_filename_exp5(result="allGDloss")
        #     final_GDloss_path = self.create_image_name(model_save_path+saved_GDloss_filename+'_token_'+token, format=".pkl")
        # else:
        #     saved_state_filename = self.create_saved_state_filename(result="allWeights")
        #     final_save_path = self.create_image_name(model_save_path+saved_state_filename+'_token_'+token, format=".pkl")
        #     saved_EMstat_filename = self.create_saved_state_filename(result="allEMstats")
        #     final_EMstat_path = self.create_image_name(model_save_path+saved_EMstat_filename+'_token_'+token, format=".csv")
        #     saved_GDloss_filename = self.create_saved_state_filename(result="allGDloss")
        #     final_GDloss_path = self.create_image_name(model_save_path+saved_GDloss_filename+'_token_'+token, format=".pkl")

        # # Initialize an empty DataFrame to store the stats
        # stats_df = pd.DataFrame(columns=[
        #     'EM_loop', 'Alpha_summation', 'EM_convergence', 
        #     'Spearman_corr_theta1_theta2', 
        #     'Spearman_corr_theta1_alpha', 
        #     'Spearman_corr_theta2_alpha',
        #     'pearson_corr_theta1_theta2',
        #     'pearson_corr_theta1_alpha', 
        #     'pearson_corr_theta2_alpha'
        # ] + [f'ELBO_sample_{i+1}' for i in range(len(self.all_read_dicts))] 
        # + [f'Convergence_sample_{i+1}' for i in range(len(self.all_read_dicts))])
        
        # # Initialize an empty DataFrame to store the GD loss history
        # GDloss_history_df = pd.DataFrame(columns=['GD_Loss'])
        
        # # comment
        # end1 = time.time()
        # interval = (end1-start1)/60
        # print(f"time_household_before_EM {interval} min")
        
        # """ EM and Gradient descent """
        # # Iterate until convergence threshold or max EM round are reached
        # #while self.convergence > self.convergence_target and self.em_round < self.max_em_rounds: ## (AT)
        # while self.em_round < self.max_em_rounds:
        #     self.convergence = 0
        #     self.em_round += 1
        #     elbo_values = {}
        #     convergence_values = {}

        #     start1 = time.time() # comment

        #     # EM
        #     sample_num = 0
        #     for sample_key in self.all_read_dicts:
        #         self.all_alpha_prime[sample_key] = self.update_alpha_prime(sample_key)
        #         self.expectation_log_theta[sample_key] = self.calculate_expectation_log_theta_m(sample_key)
        #         self.all_Phi_ri[sample_key] = self.update_Phi_ri(sample_key)
        #         self.all_theta[sample_key], convergence = self.update_expectation_theta(sample_key)
        #         self.elbo[sample_key] = self.calculate_elbo(sample_key)
        #         sample_num +=1
        #         self.convergence += convergence

        #         # Store ELBO and Convergence for this sample
        #         elbo_values[f'ELBO_sample_{sample_num}'] = self.elbo[sample_key]
        #         convergence_values[f'Convergence_sample_{sample_num}'] = convergence

        #         print(f"ELBO_sample_{sample_num} {self.elbo[sample_key]}")
        #         print(f"Convergence_sample_{sample_num} {convergence}")

        #     # comment
        #     end1 = time.time()
        #     interval = (end1-start1)/60
        #     print(f"time_EM_calculation {interval} min")

        #     start1 = time.time() # comment

        #     # UPDATE ALPHA
        #     if self.experiment_num == 4 or experiment_num == 5:
        #         self.all_alpha, GDloss_history = dirichlet_optimizer.update_alpha(self.all_unique_isoforms, self.theta_names, expectation_log_theta=self.expectation_log_theta, 
        #                                                       all_theta=self.all_theta)
            
        #     # comment
        #     end1 = time.time()
        #     interval = (end1-start1)/60
        #     print(f"time_GDloss {interval} min")

        #     print("alpha_summation ", np.sum(self.all_alpha))
        #     self.convergence = self.convergence/sample_num
        #     print(f"EM_loop {self.em_round}" )
        #     print(f"EM_convergence {self.convergence}")

        #     start1 = time.time() # comment

        #     # generate model fit correlation
        #     if self.experiment_num == 4 or experiment_num == 5:
        #         spearman_corr_theta1_theta2, pearson_corr_theta1_theta2 = self.spearman_pearson_corr('theta1_theta2')
        #         spearman_corr_theta1_alpha, pearson_corr_theta1_alpha = self.spearman_pearson_corr('theta1_alpha')
        #         spearman_corr_theta2_alpha, pearson_corr_theta2_alpha = self.spearman_pearson_corr('theta2_alpha')
        #     else:
        #         spearman_corr_theta1_theta2, pearson_corr_theta1_theta2 = 0, 0
        #         spearman_corr_theta1_alpha, pearson_corr_theta1_alpha = self.spearman_pearson_corr('theta1_alpha')
        #         spearman_corr_theta2_alpha, pearson_corr_theta2_alpha = 0, 0


        #     # comment
        #     end1 = time.time()
        #     interval = (end1-start1)/60
        #     print(f"time_corr_calculation {interval} min")

        #     start1 = time.time() # comment

        #     # Create a DataFrame for the current iteration
        #     current_stats = pd.DataFrame([{
        #         'EM_loop': self.em_round,
        #         'Alpha_summation': np.sum(self.all_alpha),
        #         'EM_convergence': self.convergence,
        #         'Spearman_corr_theta1_theta2': spearman_corr_theta1_theta2,
        #         'Spearman_corr_theta1_alpha': spearman_corr_theta1_alpha,
        #         'Spearman_corr_theta2_alpha': spearman_corr_theta2_alpha,
        #         'pearson_corr_theta1_theta2': pearson_corr_theta1_theta2,
        #         'pearson_corr_theta1_alpha': pearson_corr_theta1_alpha, 
        #         'pearson_corr_theta2_alpha': pearson_corr_theta2_alpha,
        #         **elbo_values,
        #         **convergence_values
        #     }])

        #     # Concatenate the current iteration's DataFrame with the main DataFrame
        #     stats_df = pd.concat([stats_df, current_stats], ignore_index=True)

        #     if self.experiment_num == 4 or experiment_num == 5:
        #         # Create a DataFrame for the loss history
        #         GDloss_history_current = pd.DataFrame({
        #             'GD_Loss': GDloss_history
        #         })
        #         # Append the loss history of the current iteration to the loss history DataFrame
        #         GDloss_history_df = pd.concat([GDloss_history_df, GDloss_history_current], ignore_index=True)
                
        #         # Check if the directory exists
        #         dir_path = os.path.dirname(final_GDloss_path)
        #         if not os.path.exists(dir_path):
        #             os.makedirs(dir_path)
                    

        #         with open(final_GDloss_path, 'wb') as f:
        #             pickle.dump(GDloss_history_df, f)

        #     # Save the state after each iteration
        #     stats_df.to_csv(final_EMstat_path, index=False)
        #     self.save_state(final_save_path)

        #     # comment
        #     end1 = time.time()
        #     interval = (end1-start1)/60
        #     print(f"time_rest {interval} min")


        # print(f"Exit EM loop after {self.em_round} rounds")
        # print(f"Convergence value: {self.convergence}")
        # if not self.convergence <= self.convergence_target:
        #     print(f"Convergence target {self.convergence_target} could not be reached after {self.max_em_rounds} rounds")
    
        # # Write out results
        # print("Summarize data")
        # for sample, theta_array in self.all_theta.items():
        #     print(f"Processing {sample}")
            
        #     # Retrieve the isoform names for the current sample
        #     theta_names = self.theta_names[sample]

        #     # Convert `theta_array` and `theta_names` into a dictionary
        #     theta_dict = {name: value for name, value in zip(theta_names, theta_array)}

        #     # Create DataFrame using the theta dictionary
        #     count_df = pd.DataFrame.from_dict(theta_dict, orient='index', columns=["raw"])
        #     count_df.index.name = "transcript_name"

        #     print("Compute estimated counts and TPM")
        #     # Adjusted to use the length of read_dict for the current sample
        #     count_df["est_count"] = count_df["raw"] * len(self.all_read_dicts[sample])
        #     count_df["tpm"] = count_df["raw"] * 1000000 / sum(theta_dict.values())

        #     # Add extra transcript info if required
        #     if self.extra_tx_info:
        #         tx_df = self._get_tx_df()
        #         count_df = pd.merge(count_df, tx_df, left_index=True, right_index=True, how="outer")

        #     # Cleanup and sort
        #     count_df.sort_values(by="raw", ascending=False, inplace=True)
        #     count_df.fillna(value=0, inplace=True)

        #     # Create file names and save the DataFrame to a TSV file
        #     if self.experiment_num == 5:
        #         file_name = self.create_saved_state_filename_exp5(result=self.count_file + '_' + sample)
        #     else:
        #         file_name = self.create_saved_state_filename(result=self.count_file + '_' + sample)

            
        #     file_name_timestamp = self.create_image_name(file_name, format="")
        #     count_file = f"{file_name_timestamp}.tsv" if self.count_file else None

        #     if count_file:
        #         print(f"Write file for {sample}")
                
        #         dir_path = os.path.dirname(count_file)
        #         if not os.path.exists(dir_path):
        #             os.makedirs(dir_path)
        #         count_df.to_csv(count_file, sep="\t")
            
        # end = time.time()
        # interval = (end-start)/60
        # print(f"time_total {interval} min")



     # ~~~~~~~~~~~~~~NEW FUNCTIONS~~~~~~~~~~~~~~ #
    def create_saved_state_filename_exp5(self, result):
        
        # Iterate through the list of file names
        for index, sample_file_path in enumerate(self.file_names_list, start=1):
            result += f"_file{index}_"
            for file_path in enumerate(sample_file_path, start=1):
                # Split the file path by '/' and take the last part (the file name)
                file_name = file_path[1].split('/')[-1]
                # Extract a specific part of the file name if necessary (e.g., removing extension)
                file_identifier = ''.join(file_name.split('_')).split('.')[0]
                # Construct the string
                result += f"{file_identifier}"
        result = f"{result}_GDlr_{self.GD_lr}_AlphaInitial_{self.alpha_initial}_EMround_{self.max_em_rounds}"
        return result
    
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
            model_weight1 = np.array(self.all_theta['sample1'])  # Array of weights for sample1
            model_weight2 = np.array(self.all_theta['sample2'])  # Array of weights for sample2
            theta_names_1 = np.array(self.theta_names['sample1'])
            theta_names_2 = np.array(self.theta_names['sample2'])
        elif comparison_criteria == 'theta1_alpha':
            model_weight1 = np.array(self.all_theta['sample1'])
            model_weight2 = np.array(self.all_alpha)
            theta_names_1 = np.array(self.theta_names['sample1'])
            theta_names_2 = np.array(self.all_unique_isoforms)  # Assuming `all_alpha` aligns with all unique isoforms
        elif comparison_criteria == 'theta2_alpha':
            model_weight1 = np.array(self.all_theta['sample2'])
            model_weight2 = np.array(self.all_alpha)
            theta_names_1 = np.array(self.theta_names['sample2'])
            theta_names_2 = np.array(self.all_unique_isoforms)
            
        # Normalize `model_weight2` if necessary
        if comparison_criteria in ['theta1_alpha', 'theta2_alpha']:
            model_weight2 = model_weight2/np.sum(model_weight2)
        
        # Find the common isoforms and their indices in both arrays
        common_isoforms, idx1, idx2 = np.intersect1d(theta_names_1, theta_names_2, return_indices=True)

        # Extract the common values from the arrays using the indices
        weight1_values = model_weight1[idx1]
        weight2_values = model_weight2[idx2]

        # Apply logarithmic transformation
        matched_tpm = np.log(weight1_values * 1000000  + 1)
        matched_molarity = np.log(weight2_values * 1000000  + 1)

        # Calculate Spearman's and Pearson's correlation using the matched arrays
        spearman_correlation, _ = spearmanr(matched_tpm, matched_molarity)
        pearson_correlation, _ = pearsonr(matched_tpm, matched_molarity)

        # Output the results
        if comparison_criteria == 'theta1_theta2':
            print(f'Spearman_corr_theta1_theta2 {spearman_correlation}')
            print(f'Pearson_corr_theta1_theta2 {pearson_correlation}')
        elif comparison_criteria == 'theta1_alpha':
            print(f'Spearman_corr_theta1_alpha {spearman_correlation}')
            print(f'Pearson_corr_theta1_alpha {pearson_correlation}')
        elif comparison_criteria == 'theta2_alpha':
            print(f'Spearman_corr_theta2_alpha {spearman_correlation}')
            print(f'Pearson_corr_theta2_alpha {pearson_correlation}')

        return spearman_correlation, pearson_correlation

    def spearman_corr(self, comparison_criteria):
        
        if comparison_criteria == 'theta1_theta2':
            model_weight1 = self.all_theta['sample1']
            model_weight2 = self.all_theta['sample2']
        elif comparison_criteria == 'theta1_alpha':
            model_weight1 = self.all_theta['sample1']
            model_weight2 = self.all_alpha
        elif comparison_criteria == 'theta2_alpha':
            model_weight1 = self.all_theta['sample2']
            model_weight2 = self.all_alpha

        # Convert dictionary to DataFrame
        if isinstance(model_weight1, dict):
            model_weight1 = pd.DataFrame(list(model_weight1.items()), columns=['transcript_name', 'tpm'])
            model_weight2 = pd.DataFrame(list(model_weight2.items()), columns=['transcript_name', 'tpm'])
        
        # if we are correlating theta with alpha, alpha needs to be normalized
        if comparison_criteria == 'theta1_alpha' or comparison_criteria == 'theta2_alpha':
            model_weight2['tpm'] = model_weight2['tpm'] / model_weight2['tpm'].sum()

        # Initialize lists to store matched TPM and molarity values
        matched_tpm = []
        matched_molarity = []

        # Iterate over first variable to find matching molarity values in second variable
        for index, row in model_weight1.iterrows():
            # Extract the cleaned isoform name and tpm value
            cleaned_name = row['transcript_name']
            tpm_value = row['tpm']

            # Find the corresponding molarity value in second variable
            molarity_value = model_weight2.loc[model_weight2['transcript_name'] == cleaned_name, 'tpm'].values

            # If a matching isoform is found in second variable
            if len(molarity_value) > 0:
                # Append the tpm and molarity to the respective lists
                matched_tpm.append(np.log(tpm_value * 1000000  + 1))
                matched_molarity.append(np.log(molarity_value[0]* 1000000+1))  # Take the first match in case of multiple

        # Calculate Spearman's correlation using the matched lists
        correlation, p_value = spearmanr(matched_tpm, matched_molarity)

        # Output the results
        if comparison_criteria == 'theta1_theta2':
            print(f'Spearman_corr_theta1_theta2 {correlation}')
        elif comparison_criteria == 'theta1_alpha':
            print(f'Spearman_corr_theta1_alpha {correlation}')
        elif comparison_criteria == 'theta2_alpha':
            print(f'Spearman_corr_theta2_alpha {correlation}')
        
        return correlation
        

    def calculate_elbo(self, sample_key):
        """
        Calculate the Evidence Lower BOund (ELBO) for a given sample.

        Parameters:
        sample_key (str): The key for the sample in the dataset.

        Returns:
        float: The ELBO value.
        """
       
        # Extract necessary data
        Phi_nm = self.all_Phi_ri[sample_key]  # 2D array of size (num_reads, num_isoforms)
        Pnm = self.all_read_iso_prob[sample_key]  # 2D array of size (num_reads, num_isoforms)
        expectation_log_theta_m = self.expectation_log_theta[sample_key]  # 1D array of size (num_isoforms)
        theta_names = self.theta_names[sample_key]  # Array of isoform names for the specific sample
        rnames = self.all_rnames[sample_key]  # Array of isoform names for each read
        
        # Small epsilon to prevent log(0) and division by zero
        epsilon = 1e-10

        #### Calculate the first component: ####
        # \sum_{n=1}^{N} \sum_{m=1}^{M} \phi_{nm} \left( \log \frac{p_{nm}}{\phi_{nm}} + \psi(\alpha_m) - \psi\left(\sum_{m=1}^{M} \alpha_m'\right) \right)

        # Step 1: Ensure `theta_names` is sorted and create an index mapping using `np.searchsorted`
        sorted_theta_names = np.array(theta_names)
        sort_indices = np.argsort(sorted_theta_names)
        sorted_theta_names = sorted_theta_names[sort_indices]
        # Map `rnames` to indices in `sorted_theta_names`
        isoform_indices = np.searchsorted(sorted_theta_names, rnames)
        valid_mask = (isoform_indices < len(sorted_theta_names)) & (sorted_theta_names[isoform_indices] == rnames)
        # Step 2: Filter valid indices for the isoforms
        valid_indices = sort_indices[isoform_indices[valid_mask]]
        valid_Pnm = Pnm[valid_mask]
        valid_Phi_nm = Phi_nm[valid_mask]
        log_theta_values = expectation_log_theta_m[valid_indices]
        # Adjust `Phi_nm` values to avoid log(0) and division by zero
        Phi_nm_adjusted = np.where(valid_Phi_nm == 0, epsilon, valid_Phi_nm)
        # Step 3: Calculate the ELBO components: log(p_nm / phi_nm) + expectation_log_theta_m
        log_term = np.log(valid_Pnm / Phi_nm_adjusted) + log_theta_values
        elbo_contributions = valid_Phi_nm * log_term

        elbo = np.sum(elbo_contributions)
        
        
        #### Calculate the second component:   ####
        # \log \frac{\Gamma(\sum_{m=1}^{M} \alpha_m)}{\Gamma(\sum_{m=1}^{M} \alpha_m')}

        alpha_prime = self.all_alpha_prime[sample_key]  # Dictionary of alpha prime values
        unique = self.all_unique_isoforms  # Array of isoform names for `alpha`
        alpha_values = self.all_alpha  # Array of alpha values corresponding to `unique`
        # Step 1: Find indices of `theta_names` in `unique`
        isoform_indices = np.searchsorted(unique, theta_names)
        # Step 2: Extract alpha values that correspond to `theta_names`
        selected_alpha_values = alpha_values[isoform_indices]
        # Step 3: Compute sums
        sum_alpha_prime = np.sum(alpha_prime)
        sum_alpha = np.sum(selected_alpha_values)
        # Step 4: Compute log-gamma sums
        log_gamma_sum_alpha = gammaln(sum_alpha)
        log_gamma_sum_alpha_prime = gammaln(sum_alpha_prime)

        elbo += log_gamma_sum_alpha - log_gamma_sum_alpha_prime
    

        ### Calculate the third component: ###
        # \sum_{m=1}^{M} \log \frac{\Gamma(\alpha_m')}{\Gamma(\alpha_m)}
        # Step 3: Calculate gammaln for both alpha_prime and alpha
        log_gamma_alpha_prime = gammaln(alpha_prime)
        log_gamma_alpha = gammaln(selected_alpha_values)
        # Step 4: Compute the ELBO contribution
        elbo += np.sum(log_gamma_alpha_prime - log_gamma_alpha)


        ### Calculate the fourth component: ###
        # \sum_{m=1}^{M} (\alpha_m - \alpha_m') \left( \expectation_log_theta[isoform] \right)
        differences = selected_alpha_values - alpha_prime
        elbo += np.sum(differences * expectation_log_theta_m)

        return elbo


    def update_alpha_prime(self, sample_idx):

        """
        Implements following equation: \alpha_m' = \alpha_m + \sum_n \phi_{nm}
        """
        
        # Extract current alpha, theta, and other sample-specific values
        all_alpha = self.all_alpha
        all_Phi_ri = self.all_Phi_ri[sample_idx]
        all_rnames = self.all_rnames[sample_idx]
        all_theta = self.theta_names[sample_idx]

        # Ensure `self.all_unique_isoforms` is sorted for `np.searchsorted`
        sorted_unique_isoforms = np.array(self.all_unique_isoforms)

        # Map `all_theta` isoforms to indices using `np.searchsorted`
        theta_indices = np.searchsorted(sorted_unique_isoforms, all_theta)
        valid_theta_mask = np.isin(all_theta, sorted_unique_isoforms)
        theta_indices = theta_indices[valid_theta_mask]

        # Initialize `all_alpha_prime` directly using `all_alpha`
        all_alpha_prime = np.zeros_like(all_alpha)
        all_alpha_prime[theta_indices] = all_alpha[theta_indices]

        # Map `all_rnames` to indices using `np.searchsorted`
        rname_indices = np.searchsorted(sorted_unique_isoforms, all_rnames)
        valid_rname_mask = np.isin(all_rnames, sorted_unique_isoforms)
        valid_rname_indices = rname_indices[valid_rname_mask]
        valid_Phi_ri = all_Phi_ri[valid_rname_mask]

        all_alpha_prime = all_alpha_prime.astype(float)

        # Accumulate `phi` values into `all_alpha_prime`
        np.add.at(all_alpha_prime, valid_rname_indices, valid_Phi_ri)

        # Return the subset of `all_alpha_prime` corresponding to isoforms in `all_theta`
        return all_alpha_prime[theta_indices]

    def update_expectation_theta(self, sample_key):
        """
        Calculate the expected value E_Q(theta)[theta_m] for a Dirichlet distribution.

        Parameters:
        sample_key (str): The key for the sample in the dataset.

        Returns:
        Counter: The expected theta values as a Counter.
        """
        # Extract data as NumPy arrays
        theta_old = self.all_theta[sample_key]  # Assuming this is a 1D NumPy array
        alpha_prime = self.all_alpha_prime[sample_key]  # Assuming this is a 1D NumPy array

        # Calculate sum of alpha_prime
        sum_alpha_prime = np.sum(alpha_prime)

        # Step 1: Calculate expected_theta using vectorized operations
        expected_theta = alpha_prime / sum_alpha_prime

        # Step 2: Calculate convergence if em_round > 0
        convergence = 0
        if self.em_round > 0:
            convergence = np.sum(np.abs(theta_old - expected_theta))

        return expected_theta, convergence

    def calculate_expectation_log_theta_m(self, sample_key):
        """
        Calculate the expected log value E_Q(theta)[log(theta_m)] for a Dirichlet distribution.

        Parameters:
        sample_key (str): The key for the sample in the dataset.

        Returns:
        Counter: The expected log theta values as a Counter.
        """
        # Extract alpha_prime as a NumPy array for the given sample
        alpha_prime = self.all_alpha_prime[sample_key]  # Assuming this is a NumPy array
        sum_alpha_prime = np.sum(alpha_prime)

        # Calculate the expected log theta values using vectorized psi operations
        eq_log_theta_m = psi(alpha_prime) - psi(sum_alpha_prime)

        return eq_log_theta_m

    def update_Phi_ri(self, sample_key):

        """
        Update phi_nm using the equation:
        phi_nm ∝ p_nm * exp(E_Q(theta)[log(theta_m)])
        """
        # Extract the necessary data
        p_nm_matrix = self.all_read_iso_prob[sample_key]  # 1D array of size (172619,)
        rnames = self.all_rnames[sample_key]  # Array of isoform names for each read
        expectation_log_theta_m = self.expectation_log_theta[sample_key]  # 1D array matching `theta_names`
        theta_names = self.theta_names[sample_key]  # Isoform names for the specific sample
        read_ids = self.all_readName[sample_key]  # Contains read IDs

        # Step 1: Ensure `theta_names` is sorted and create an index mapping using `np.searchsorted`
        sorted_theta_names = np.array(theta_names)
        sort_indices = np.argsort(sorted_theta_names)
        sorted_theta_names = sorted_theta_names[sort_indices]

        # Find the indices in `sorted_theta_names` corresponding to `rnames`
        isoform_indices = np.searchsorted(sorted_theta_names, rnames)
        valid_mask = (isoform_indices < len(sorted_theta_names)) & (sorted_theta_names[isoform_indices] == rnames)

        # Step 2: Calculate unnormalized phi values for valid indices
        valid_indices = sort_indices[isoform_indices[valid_mask]]
        log_theta_values = expectation_log_theta_m[valid_indices]
        unnormalized_phi = p_nm_matrix[valid_mask] * np.exp(log_theta_values)

        # Step 3: Normalize phi values for each read using `read_ids`
        valid_read_ids = read_ids[valid_mask]

        # Sum unnormalized phi values per unique read using `np.bincount`
        unique_read_ids, inverse_indices = np.unique(valid_read_ids, return_inverse=True)
        sum_phi_per_read = np.bincount(inverse_indices, weights=unnormalized_phi)

        # Normalize unnormalized phi values by their respective read sums
        normalized_phi_values = unnormalized_phi / sum_phi_per_read[inverse_indices]

        # Step 4: Insert the normalized values back into the full `normalized_phi` array
        normalized_phi = np.zeros_like(p_nm_matrix)
        normalized_phi = normalized_phi.astype(float)
        normalized_phi[valid_mask] = normalized_phi_values

        return normalized_phi


    def assign_alpha(self):
        # Step 1: Use `self.all_unique_isoforms` directly
        unique_isoforms = self.all_unique_isoforms

        # Step 2: Number of unique isoforms
        n = len(unique_isoforms)

        # Step 3: Initialize alpha values in an unconstrained space
        np.random.seed(5)  # Set the fixed seed for reproducibility
        alpha_raw = np.random.randn(n)

        # Step 4: Apply the softmax function to ensure they sum to 1
        alpha_softmax = np.exp(alpha_raw) / np.sum(np.exp(alpha_raw))

        # Step 5: Scale the resulting values to have the desired fixed sum
        alpha_scaled = alpha_softmax * self.alpha_initial

        # Step 6: Return the alpha values
        return alpha_scaled
    
    def assign_alpha_constant(self):
        # Step 1: Use `self.all_unique_isoforms` directly
        unique_isoforms = self.all_unique_isoforms

        # Step 2: Number of unique isoforms
        n = len(unique_isoforms)

        # Step 5: Scale the resulting values to have the desired fixed sum
        alpha_scaled = np.full(n, 3)

        # Step 6: Return the alpha values
        return alpha_scaled
    

    def calculate_theta_and_alpha_prime_0(self, all_rnames, compatibility_list):
        """
        Calculates the initial model parameter (isoform percentage) using vectors.
        """
        # Identify unique reference names and assign indices
        unique_refs, inverse_indices = np.unique(all_rnames, return_inverse=True)
        
        # Calculate total abundance scores for each unique reference name
        abundance_scores = np.zeros(len(unique_refs))
        abundance_scores = abundance_scores.astype(float)

        np.add.at(abundance_scores, inverse_indices, compatibility_list)

        # Normalize the abundance scores
        total_score = np.sum(abundance_scores)
        normalized_abundances = abundance_scores / total_score

        return unique_refs, normalized_abundances

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


    def calculate_Z(self, all_read_names, all_rnames, compatibility_list, p_nm, all_theta):
        """
        Vectorized EM assignment to update read-transcript compatibility based on transcript abundances.
        """
        # Find unique reads and map them to indices
        unique_reads, read_indices = np.unique(all_read_names, return_inverse=True)

        # Map isoforms (rnames) to indices in unique_refs (isoforms used in theta)
        unique_rnames, isoform_indices = np.unique(all_rnames, return_inverse=True)

        # Get theta values for each isoform
        theta_for_isoforms = all_theta[isoform_indices]

        # Step 1: Calculate denominator for each read
        compat_theta_product = compatibility_list * theta_for_isoforms
        denominators = np.zeros(len(unique_reads))
        denominators = denominators.astype(float)

        np.add.at(denominators, read_indices, compat_theta_product)

        # Step 2: Calculate Z_ri values
        Z_values = (p_nm * theta_for_isoforms) / denominators[read_indices]

        # Step 3: Normalize Z_values for each read
        normalization_factors = np.zeros(len(unique_reads))
        normalization_factors = normalization_factors.astype(float)

        np.add.at(normalization_factors, read_indices, Z_values)
        Z_values /= normalization_factors[read_indices]

        # Return Z values as a structured array
        return Z_values


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
        For every read, this function provides the compatible isoforms and normalizes them by N-K+1, returning vectors.
        """
        # Fetch reference lengths and names
        read_dict = self.all_read_dicts[sample_key]
        ref_len_dict = self.all_ref_len_dicts[sample_key]

    
        # Create a mapping of reference names to indices for quick lookup
        ref_len_dict_keys = np.array(list(ref_len_dict.keys()))
        
        ref_lengths = np.fromiter((ref_len_dict[key] for key in ref_len_dict_keys), dtype=float)
        ref_name_to_index = dict(zip(ref_len_dict_keys, np.arange(len(ref_len_dict_keys))))

        all_read_names = []
        all_alignments = []
        all_alignment_lens = []
        all_rnames = []

        # Collect all data into lists for vectorized processing
        for read_name, read in read_dict.items():
            if not self.is_iterable(read):
                read = [read]  # Wrap non-iterable read in a list

            for alignment in read:
                alignments = alignment.alignment_list
                alignment_lens = np.array([string.align_len for string in alignments])
                rnames = [string.rname for string in alignments]

                # Collect the read name and alignment details for each alignment
                all_read_names.extend([read_name] * len(alignments))
                all_alignments.extend(alignments)
                all_alignment_lens.extend(alignment_lens)
                all_rnames.extend(rnames)

        # Convert all collected data into numpy arrays for vectorized processing
        all_alignment_lens = np.array(all_alignment_lens)
        all_rnames = np.array(all_rnames)
        all_read_names = np.array(all_read_names)

        # Sort all arrays based on `all_read_names` and `all_rnames` to maintain consistent ordering
        sort_indices = np.lexsort((all_rnames, all_read_names))  # Sort by `all_read_names`, then by `all_rnames`
        all_read_names = all_read_names[sort_indices]
        all_rnames = all_rnames[sort_indices]
        all_alignment_lens = all_alignment_lens[sort_indices]

        # Find the indices of reference names in ref_len_dict_keys
        ref_indices = np.vectorize(ref_name_to_index.get)(all_rnames)
        ref_lengths_for_alignments = ref_lengths[ref_indices]

        # Vectorized compatibility scores
        unique_read_names, counts = np.unique(all_read_names, return_counts=True)
        compatibility_list = np.repeat(1 / counts, counts)  # Compatibility score is the inverse of the count of alignments for each read

        # Vectorized probability calculation
        denominators = ref_lengths_for_alignments - all_alignment_lens + 1
        p_nm = np.where(denominators <= 0, np.nan, 1 / denominators)

        # Handle problematic probabilities
        problematic_indices = np.isnan(p_nm) | (p_nm < 0)
        p_nm[problematic_indices] = 1  # Replace problematic probabilities with 1
        
        return all_read_names, all_rnames, compatibility_list, p_nm
    

    def initialize_model(self):
        start = time.time()

        # Loop over all file names provided and parses the reads with
        for index, file_name in enumerate(self.file_names_list, start=1):
            main_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/new_minimap_file'
            file_name = main_dir+'/'+file_name
            print(f"sample_{index} {file_name}")
            # Parse the BAM file
            if self.parse_original:
                read_dict, ref_len_dict = self._parse_bam_original(file_name=file_name)
            else:
                read_dict, ref_len_dict = self._parse_bam(file_name=file_name)


            # Check if the read_dict file already exists
            ## (AT)
            main_dir =  self.process_dir 
            pkl_file_name_read_dict = file_name.split('/')[-1].split('.')[0] + '_read_dicts.pkl'
            pkl_file_name_ref_len_dict = file_name.split('/')[-1].split('.')[0] + '_ref_len_dicts.pkl'

            read_dict_path = os.path.join(main_dir, pkl_file_name_read_dict)
            if not os.path.exists(read_dict_path):
                with open(read_dict_path, 'wb') as file:
                    pickle.dump(read_dict, file)
                    print(f"Saved {read_dict_path}")
            else:
                print(f"File {read_dict_path} already exists. Skipping save.")

            # Check if the ref_len_dict file already exists
            ref_len_dict_path = os.path.join(main_dir, pkl_file_name_ref_len_dict)
            if not os.path.exists(ref_len_dict_path):
                with open(ref_len_dict_path, 'wb') as file:
                    pickle.dump(ref_len_dict, file)
                    print(f"Saved {ref_len_dict_path}")
            else:
                print(f"File {ref_len_dict_path} already exists. Skipping save.")


            # Store the dictionaries in the universal dictionary with a sample key
            sample_key = f'sample{index}'
            self.all_read_dicts[sample_key] = read_dict
            self.all_ref_len_dicts[sample_key] = ref_len_dict
        

    # ~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~ #

    def _parse_bam_original(self, file_name):
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
        # with pysam.AlignmentFile(self.alignment_file) as bam:
            # Collect reference lengths in dict
            for name, length in zip(bam.references, bam.lengths):
                ref_len_dict[name] = length

            for idx, alignment in enumerate(bam):
                if alignment.is_unmapped:
                    c["Discarded unmapped alignments"] += 1
                elif not self.keep_neg_strand and alignment.is_reverse:
                    c["Discarded negative strand alignments"] += 1
                elif not self.keep_supplementary and alignment.is_supplementary:
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
        return filtered_read_dict, ref_len_dict

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
        


"""

        # Loop over all file names provided and parses the reads with
        # for index, file_name in enumerate(self.file_names_list, start=1):
        #     print(f"sample_{index} {file_name}")
        #     # Parse the BAM file
        #     read_dict, ref_len_dict = self._parse_bam(file_name=file_name)


        #     # Check if the read_dict file already exists
        #     
        #     main_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln_pklfiles'
        #     pkl_file_name_read_dict = file_name.split('/')[-1].split('.')[0] + '_read_dicts.pkl'
        #     pkl_file_name_ref_len_dict = file_name.split('/')[-1].split('.')[0] + '_ref_len_dicts.pkl'

        #     read_dict_path = os.path.join(main_dir, pkl_file_name_read_dict)
        #     if not os.path.exists(read_dict_path):
        #         with open(read_dict_path, 'wb') as file:
        #             pickle.dump(read_dict, file)
        #             print(f"Saved {read_dict_path}")
        #     else:
        #         print(f"File {read_dict_path} already exists. Skipping save.")

        #     # Check if the ref_len_dict file already exists
        #     ref_len_dict_path = os.path.join(main_dir, pkl_file_name_ref_len_dict)
        #     if not os.path.exists(ref_len_dict_path):
        #         with open(ref_len_dict_path, 'wb') as file:
        #             pickle.dump(ref_len_dict, file)
        #             print(f"Saved {ref_len_dict_path}")
        #     else:
        #         print(f"File {ref_len_dict_path} already exists. Skipping save.")


        #     # Store the dictionaries in the universal dictionary with a sample key
        #     sample_key = f'sample{index}'
        #     self.all_read_dicts[sample_key] = read_dict
        #     self.all_ref_len_dicts[sample_key] = ref_len_dict
"""