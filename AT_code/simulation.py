import numpy as np
import hashlib
import math
import random
from collections import defaultdict
import os
import pickle
import time
from scipy.linalg import cholesky
from scipy.stats import spearmanr, rankdata
from NanoCount.Read import Read
from NanoCount.common import *


timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")

def calculate_theta_and_alpha_prime_0(all_Zri):
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

# Function to generate samples from a Dirichlet distribution
# def generate_dirichlet_samples(alpha, num_samples):
#     """
#     Generate samples from a Dirichlet distribution.
    
#     Args:
#         alpha (list or array): The alpha parameters for the Dirichlet distribution.
#         num_samples (int): Number of samples to generate.
    
#     Returns:
#         list: List of generated samples, where each sample is a list of 150k elements.
#     """
#     return [np.random.dirichlet(alpha) for _ in range(num_samples)]

def generate_correlated_sample(sample1, desired_corr):
    """
    Generate a sample that has a desired Spearman correlation with sample1.
    
    Args:
        sample1 (array): The first sample (used as the reference for correlation).
        desired_corr (float): The desired Spearman correlation with sample1.
    
    Returns:
        array: A second sample that has the desired Spearman correlation with sample1.
    """
    # Step 1: Generate an initial uncorrelated sample
    sample2 = np.random.rand(len(sample1))
    
    # Step 2: Rank the first and second samples
    rank1 = rankdata(sample1)
    rank2 = rankdata(sample2)
    
    # Step 3: Adjust the ranks of the second sample to achieve the desired Spearman correlation
    rank2_adjusted = rank1 * desired_corr + rank2 * (1 - desired_corr)
    
    # Step 4: Create a new sample2 by assigning values based on the new ranks
    sample2_adjusted = np.interp(rank2_adjusted, np.sort(rank2), np.sort(sample2))

    # Step 5: Normalize sample2 so it sums to 1 (to meet Dirichlet properties)
    sample2_normalized = np.abs(sample2_adjusted) / np.sum(np.abs(sample2_adjusted))
    
    return sample2_normalized


def generate_correlated_dirichlet_samples(alpha, desired_corr, num_samples):
    """
    Generate a list of Dirichlet samples with the first two samples having a fixed Spearman correlation.
    
    Args:
        alpha (list or array): The alpha parameters for the Dirichlet distribution.
        desired_corr (float): Desired Spearman correlation between the first two samples.
        num_samples (int): Number of samples to generate.
    
    Returns:
        list: List of generated samples, where the first two samples have the desired Spearman correlation.
    """
    # Step 1: Generate the first sample from Dirichlet
    sample1 = np.random.dirichlet(alpha)
    
    # Step 2: Generate a second sample that has the desired Spearman correlation with the first sample
    sample2 = generate_correlated_sample(sample1, desired_corr)
    
    # Generate remaining samples normally from Dirichlet distribution
    remaining_samples = [np.random.dirichlet(alpha) for _ in range(num_samples - 2)]
    
    # Return all samples with the first two having the desired correlation
    return [sample1, sample2] + remaining_samples

def generate_dirichlet_samples(alpha, desired_corr, num_samples):
    """
    Generate a list of Dirichlet samples with the first two samples having a fixed Spearman correlation.
    
    Args:
        alpha (list or array): The alpha parameters for the Dirichlet distribution.
        desired_corr (float): Desired Spearman correlation between the first two samples.
        num_samples (int): Number of samples to generate.
    
    Returns:
        list: List of generated samples, where the first two samples have the desired Spearman correlation.
    """
    # Step 1: Generate the first sample from Dirichlet
    sample1 = np.random.dirichlet(alpha)
    
    # Step 2: Generate a second sample that has the desired Spearman correlation with the first sample
    sample2 = sample1
    
    # Generate remaining samples normally from Dirichlet distribution
    remaining_samples = [np.random.dirichlet(alpha) for _ in range(num_samples - 2)]
    
    # Return all samples with the first two having the desired correlation
    return [sample1, sample2] + remaining_samples

# Function to generate a unique hash for each element
def generate_isoform_hash(index, element_num):
    """
    Generate a unique truncated hash for a given index using SHA-256.
    
    Args:
        index (str): The string to hash.
        length (int): Number of characters to return (default is 3).
        element_num (int): The number of isoforms (unique items).
        min_length (int): The minimum number of characters required.
    
    Returns:
        str: A unique hash string of the specified length.
    """
    # Calculate the minimum number of characters needed
    extra_length = 4
    min_length = math.ceil(math.log(element_num) / math.log(16)) + extra_length
    hash_object = hashlib.sha256(str(index).encode())
    return hash_object.hexdigest()[:min_length]

# Function to create a dictionary of elements with unique hashes for each sample
def create_sample_hashes(samples, element_num):
    """
    Create a dictionary with consistent hashed names for each isoform across all samples.
    Each isoform will have the same hash across samples.

    Args:
        samples (list): List of generated samples.
        element_num (int): Minimum number of characters required for the hash.
    
    Returns:
        dict: Dictionary where keys are sample names and values are dictionaries of hashed elements.
    """
    # Create a dictionary for all unique isoform hashes
    isoform_hashes = {}

    # Generate consistent unique hashes for each isoform index (only once)
    num_isoforms = len(samples[0])  # Assume all samples have the same number of isoforms
    for idx in range(num_isoforms):
        unique_hash = generate_isoform_hash(f"isoform_{idx}", element_num)
        isoform_hashes[idx] = unique_hash

    sample_dict = {}

    # Loop over each sample
    for i, sample in enumerate(samples):
        sample_name = f"sample_{i+1}"
        sample_dict[sample_name] = {}

        # Use pre-generated hashes to ensure consistency across samples
        for idx in range(num_isoforms):
            unique_hash = isoform_hashes[idx]  # Use the pre-generated hash for this isoform
            value = sample[idx]  # Get the isoform value for this sample
            sample_dict[sample_name][unique_hash] = value  # Ensure all isoforms are included

    return sample_dict

def assign_alpha(element_num, alpha_initial):
        
        # Set the fixed seed for reproducibility
        np.random.seed(5)
        
        # Step 1: Initialize alpha values in an unconstrained space
        alpha_raw = np.random.randn(element_num)
        
        # Step 2: Apply the softmax function to ensure they sum to 1
        alpha_softmax = np.exp(alpha_raw) / np.sum(np.exp(alpha_raw))
        
        # Step 3: Scale the resulting values to have the desired fixed sum
        alpha_scaled = alpha_softmax * alpha_initial
        
        return alpha_scaled

def find_common_isoforms(sample_dict, sample1_name, sample2_name):
    """
    Find the common isoform names (hashes) between two samples.
    
    Args:
        sample_dict (dict): Dictionary of samples with hashed isoform names.
        sample1_name (str): Name of the first sample (e.g., "sample_1").
        sample2_name (str): Name of the second sample (e.g., "sample_2").
    
    Returns:
        set: A set of common isoform names (hashes) between the two samples.
    """
    sample1_isoforms = set(sample_dict[sample1_name].keys())
    sample2_isoforms = set(sample_dict[sample2_name].keys())
    
    # Find the intersection of the two sets
    common_isoforms = sample1_isoforms.intersection(sample2_isoforms)
    
    return common_isoforms

# Function to generate a unique read hash
def generate_read_hash(index):
    """
    Generate a unique hash for a read using SHA-256.
    
    Args:
        index (str): The string to hash.
    
    Returns:
        str: A unique hash string.
    """
    hash_object = hashlib.sha256(str(index).encode())
    return hash_object.hexdigest()

# Function to generate read counts proportionate to isoform proportions for multiple samples
def generate_read_counts_multiple_samples(sample_dicts, total_reads_per_sample, isoform_lengths, simulation_dir, 
                                          min_ambiguity, max_ambiguity, min_read_length, avg_read_length, max_read_length, read_name):
    """
    Generate read counts for each isoform in multiple samples proportionate to their presence in the sample.
    
    Args:
        sample_dicts (dict): Dictionary where each key is a sample name, and the value is a dictionary of isoform proportions.
        total_reads_per_sample (int): Total number of reads to generate for each sample.
    
    Returns:
        dict: Dictionary where each key is a sample, and the value is a dictionary of reads and their assigned isoforms.
    """

    # Loop over each sample in sample_dicts
    i = 0
    for sample_name, isoform_proportions in sample_dicts.items():
        start = time.time()
        print(f"Generating reads for {sample_name}...")
        # Generate read counts for this sample
        # read_dict, read_probabilities = generate_read_counts(isoform_proportions, total_reads_per_sample, 
        #                         min_ambiguity, max_ambiguity, isoform_lengths, min_read_length, avg_read_length, max_read_length)
        read_dict, read_probabilities = generate_read_counts(isoform_proportions, total_reads_per_sample, 
                                min_ambiguity, max_ambiguity, isoform_lengths, min_read_length, avg_read_length, max_read_length)
    
        ## comment
        if i ==0:
            theta_S1 = calculate_theta_and_alpha_prime_0(read_dict)
            print(f"theta_S1_num {len(theta_S1)}")
        if i == 1:
            theta_S2 = calculate_theta_and_alpha_prime_0(read_dict)
            print(f"theta_S2_num {len(theta_S2)}")

        # Define the filename for saving
        file_path = os.path.join(simulation_dir, f"{read_name}{timestamp}_{sample_name}_read_dict.pkl")
        # Save the read_dict as a pickle file
        with open(file_path, 'wb') as f:
            pickle.dump(read_dict, f)

        # Define the filename for saving
        file_path = os.path.join(simulation_dir, f"{read_name}{timestamp}_{sample_name}_read_iso_prob.pkl")
        # Save the read_dict as a pickle file
        with open(file_path, 'wb') as f:
            pickle.dump(read_probabilities, f)
        i+=1
        
        end = time.time()
        interval = (end-start)/60
        print(f"time_to_generate {total_reads_per_sample} reads {interval} min")
    
    
    ## comment
    all_keys = set(theta_S1.keys()).union(set(theta_S2.keys()))
    values_S1 = [theta_S1.get(key, 0) for key in all_keys]
    values_S2 = [theta_S2.get(key, 0) for key in all_keys]
    corr, _ = spearmanr(values_S1, values_S2)
    print(f"SP corr {corr}")
    # Get the keys (isoforms) for each counter
    keys_S1 = set(theta_S1.keys())
    keys_S2 = set(theta_S2.keys())

    # Find the intersection of the two sets (common isoforms)
    common_isoforms = keys_S1.intersection(keys_S2)
    print(f"common_isoform_num {len(common_isoforms)}")
    print()
    return corr



def generate_read_counts_optimized(sample_dict, total_reads, min_ambiguity, max_ambiguity, isoform_lengths, min_read_length, avg_read_length, max_read_length, std_dev=30):
    """
    Generate read counts for each isoform proportionate to their presence in the sample using multinomial distribution,
    and assign reads with ambiguity to multiple isoforms.

    Args:
        sample_dict (dict): Dictionary with isoform proportions.
        total_reads (int): Total number of reads to generate.
        min_ambiguity (int): Minimum number of isoforms a read can be ambiguously assigned to.
        max_ambiguity (int): Maximum number of isoforms a read can be ambiguously assigned to.
        isoform_lengths (dict): Dictionary with isoform lengths.
        min_read_length (int): Minimum read length to generate.
        avg_read_length (int): Average read length for the Gaussian distribution.
        max_read_length (int): Maximum read length to generate.
        std_dev (int): Standard deviation for read length Gaussian distribution.

    Returns:
        dict: Dictionary of reads and their assigned isoforms and probabilities.
    """
    read_dict = defaultdict(dict)
    read_probabilities = defaultdict(dict)

    # Step 1: List of isoforms and their proportions
    isoforms = list(sample_dict.keys())  # List of isoforms
    proportions = np.array(list(sample_dict.values()), dtype=np.float64)  # Proportions for each isoform

    # Step 4: Generate read hashes and assign isoforms with ambiguity, respecting proportions
    for read_index in range(total_reads):
        # Generate a unique hash for each read
        read_hash = f"SRR{random.randint(60000000, 70000000)}.{random.randint(1, 10)}"

        # Randomly assign this read to a number of isoforms (between min_ambiguity and max_ambiguity)
        num_isoforms_to_assign = random.randint(min_ambiguity, max_ambiguity)

        # Sample isoforms based on their proportions with no replacement (each isoform selected once)
        assigned_isoforms = np.random.choice(isoforms, size=num_isoforms_to_assign, replace=False, p=proportions)

        # Assign equal proportions for each isoform (assuming equal contribution from each isoform to the read)
        proportion_per_isoform = 1.0 / len(assigned_isoforms)
        read_dict[read_hash] = {isoform: proportion_per_isoform for isoform in assigned_isoforms}

        # Step 5: Generate a random read length from a Gaussian distribution
        read_length = np.clip(int(np.random.normal(loc=avg_read_length, scale=std_dev)), min_read_length, max_read_length)

        # Step 6: Calculate the probabilities for each isoform based on the read length and isoform length
        isoform_probs = {}
        for isoform in assigned_isoforms:
            isoform_length = isoform_lengths.get(isoform, 0)
            denominator = isoform_length - read_length + 1

            # Ensure the denominator is positive and valid
            if denominator > 0:
                isoform_probs[isoform] = 1.0 / denominator
            else:
                isoform_probs[isoform] = 1.0  # Fallback if the denominator is invalid (e.g., if the isoform length < read length)

        # Store the probabilities for this read
        read_probabilities[read_hash] = isoform_probs
        #print(f"{read_index}_DONE")

    return read_dict, read_probabilities

# Function to generate read counts proportionate to isoform proportions for multiple samples using multinomial distribution
def generate_read_counts(sample_dict, total_reads, min_ambiguity, max_ambiguity, isoform_lengths, min_read_length, avg_read_length, max_read_length, std_dev=30):
    """
    Generate read counts for each isoform proportionate to their presence in the sample using multinomial distribution.

    Args:
        sample_dict (dict): Dictionary with isoform proportions.
        total_reads (int): Total number of reads to generate.
        min_ambiguity (int): Minimum number of isoforms a read can be ambiguously assigned to.
        max_ambiguity (int): Maximum number of isoforms a read can be ambiguously assigned to.

    Returns:
        dict: Dictionary of reads and their assigned isoforms.
    """
    read_dict = defaultdict(dict)
    read_probabilities = defaultdict(dict)

    # Step 1: Use multinomial distribution to generate read counts based on isoform proportions
    isoforms = list(sample_dict.keys())
    proportions = list(sample_dict.values())

    # Multinomial distribution gives the number of reads assigned to each isoform
    read_counts = np.random.multinomial(total_reads, proportions)

    # Step 2: Create a list of isoforms with repetition proportionate to their read count
    isoform_list = []
    for isoform, count in zip(isoforms, read_counts):
        isoform_list.extend([isoform] * count)  # Repeat isoform proportionally to its count

    # Shuffle the isoform list to simulate random assignment of reads
    random.shuffle(isoform_list)

    # Step 3: Generate read hashes and assign isoforms with ambiguity
    i=0
    for read_index in range(total_reads):
        read_hash = generate_read_hash(f"read_{read_index}")

        # Randomly assign this read to between min_ambiguity and max_ambiguity isoforms, maintaining proportions
        num_isoforms_to_assign = random.randint(min_ambiguity, max_ambiguity)
        assigned_isoforms = random.sample(isoform_list, num_isoforms_to_assign)  # Sample from isoform_list, not set
        unique_assigned_isoforms = (set(assigned_isoforms))

        # Assign isoforms with equal proportion (assuming ambiguity is equally likely between assigned isoforms)
        # proportion_per_isoform = 1.0 / num_isoforms_to_assign
        proportion_per_isoform = 1.0 / len(unique_assigned_isoforms)
        read_dict[read_hash] = {isoform: proportion_per_isoform for isoform in unique_assigned_isoforms}

        # Generate lengths from a Gaussian distribution
        lengths = int(np.random.normal(loc=avg_read_length, scale=std_dev))
        
        # Generate a random read length between min_read_length and max_read_length
        read_length = max(min_read_length, min(max_read_length, lengths))

        # Step 4: Calculate the probabilities for each isoform
        isoform_probs = {}

        for isoform in unique_assigned_isoforms:
            isoform_length = isoform_lengths.get(isoform, 0)  # Get the isoform length
            if isoform_length > read_length:
                # Calculate the probability using the formula: 1 / (isoform_length - read_length + 1)
                denominator = (isoform_length - read_length + 1)
                # Check for potential division by zero
                if denominator == 0:
                        prob = np.nan  # This will trigger the NaN check below
                else:
                    prob = 1 / denominator
                # Check if prob is NaN or less than 1
                if np.isnan(prob) or prob < 0:
                    isoform_probs[isoform] = 1
                else:
                    isoform_probs[isoform] = prob
            else:
                isoform_probs[isoform] = 1
                i+=1
                print(f"len longer {i}")

        # Store the read's isoform probabilities
        read_probabilities[read_hash] = isoform_probs

    return read_dict, read_probabilities


# Function to generate isoform lengths based on a Gaussian distribution, with same lengths for all samples
def generate_isoform_lengths(isoform_proportion_samples, avg_length=2000, std_dev=100, min_length=300):
    """
    Generate a list of isoform lengths using a Gaussian distribution, ensuring the same lengths for all samples.
    
    Args:
        isoform_proportion_samples (dict): Dictionary of isoform proportions for multiple samples.
        avg_length (int): The avg length of the isoforms (default 1800 bp).
        std_dev (int): The standard deviation for the lengths (default 300 bp).
        min_length (int): The minimum length allowed for an isoform (default 100 bp).
    
    Returns:
        dict: A dictionary with isoform names as keys and their corresponding lengths as values (shared across all samples).
    """
    # Get a unique set of isoforms from all samples
    all_isoforms = set()
    for isoforms in isoform_proportion_samples.values():
        all_isoforms.update(isoforms.keys())
    
    # Number of unique isoforms
    num_isoforms = len(all_isoforms)
    
    # Generate lengths from a Gaussian distribution
    lengths = np.random.normal(loc=avg_length, scale=std_dev, size=num_isoforms)
    
    # Apply a minimum length cutoff to avoid negative or too short lengths
    lengths = np.maximum(lengths, min_length)
    
    # Convert to integers
    lengths = lengths.astype(int)
    
    # Create a dictionary with isoform names as keys and their lengths as values
    isoform_lengths = dict(zip(all_isoforms, lengths))
    
    return isoform_lengths


# Main function
def main():

    ###### GENERATE ONLY 2 SAMPLES AT A TIME ########

    start = time.time()

    isoform_number= int(150e3)     # How many isoforms per sample
    alpha_initial = 100e3      # Initial alpha summation
    num_samples = 2            # Specify the number of samples you want to generate
    simulation_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/round4'

    short_total_reads = int(2e7)   # total number of short reads
    short_min_ambiguity = int(3)    # short_min_ambiguity, short_max_ambiguity -> amount of ambiguity eg if the numbers are 3 and 10 that means read would be compatible with at least 3 to max 10 isoforms
    short_max_ambiguity = int(10)
    short_min_read_length = int(50)     # short read minimum length
    short_avg_read_length = int(150)    # average length
    short_max_read_length =  int(400)    # maximum length

    long_total_reads = int(2e6)   # total number of long reads
    long_min_ambiguity = int(1)    # long_min_ambiguity, long_max_ambiguity -> amount of ambiguity eg if the numbers are 3 and 10 that means read would be compatible with at least 3 to max 10 isoforms
    long_max_ambiguity = int(4)
    long_min_read_length = int(1500)     # long read minimum length
    long_avg_read_length = int(1800)    # average length
    long_max_read_length =  int(2000)    # maximum length
    
    ### DUMMY
    # isoform_number= int(150)     # How many isoforms per sample
    # alpha_initial = 100      # Initial alpha summation
    # num_samples = 2             # Specify the number of samples you want to generate
    # simulation_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/trial'
    
    # short_total_reads = int(20000)   # total number of short reads
    # short_min_ambiguity = int(10)    # short_min_ambiguity, short_max_ambiguity -> amount of ambiguity eg if the numbers are 3 and 10 that means read would be compatible with at least 3 to max 10 isoforms
    # short_max_ambiguity = int(20)
    # short_min_read_length = int(50)     # short read minimum length
    # short_avg_read_length = int(150)    # average length
    # short_max_read_length =  int(400)    # maximum length

    # long_total_reads = int(2000)   # total number of long reads
    # long_min_ambiguity = int(1)    # long_min_ambiguity, long_max_ambiguity -> amount of ambiguity eg if the numbers are 3 and 10 that means read would be compatible with at least 3 to max 10 isoforms
    # long_max_ambiguity = int(6)
    # long_min_read_length = int(1500)     # long read minimum length
    # long_avg_read_length = int(1800)    # average length
    # long_max_read_length =  int(2000)    # maximum length
    
    element_num = isoform_number 
    desired_corr=0.65
       
    read_name = 'SHORT' # (AT)

    # Parameters for the Dirichlet distribution
    alpha = assign_alpha(element_num, alpha_initial)

    # Step 1: Generate Dirichlet samples
    proportion_samples = generate_dirichlet_samples(alpha, desired_corr, num_samples)
    
    # Step 2: Generate hashed names for each element in the samples
    isoform_proportion_samples = create_sample_hashes(proportion_samples, element_num)
    
    # Generate isoform lengths for all samples (same lengths across samples)
    isoform_lengths = generate_isoform_lengths(isoform_proportion_samples)

    # Generate read counts and assignments for multiple samples
    corr = generate_read_counts_multiple_samples(isoform_proportion_samples, short_total_reads, isoform_lengths, simulation_dir, 
        short_min_ambiguity, short_max_ambiguity,  short_min_read_length, short_avg_read_length, short_max_read_length, read_name = 'SHORT')
    file_path = os.path.join(simulation_dir, f"{read_name}{timestamp}_spCorr_{corr:.2f}_isoformAbundance_groundTruth.pkl")
    # Save the read_dict as a pickle file
    with open(file_path, 'wb') as f:
        pickle.dump(isoform_proportion_samples, f)



    read_name = 'LONG' # (AT)
    proportion_samples_new = generate_dirichlet_samples(alpha, desired_corr, num_samples)
    
    # Step 5: Reuse the same isoform names (hashes) and lengths, but with new proportions
    isoform_proportion_samples_new = {}
    for sample_name, isoform_proportions in isoform_proportion_samples.items():
        isoform_proportion_samples_new[sample_name] = {}
        for idx, hash_name in enumerate(isoform_proportions.keys()):
            # Assign new proportions to the same hash names (isoform names)
            isoform_proportion_samples_new[sample_name][hash_name] = proportion_samples_new[int(sample_name[-1])-1][idx]
    
    corr_new = generate_read_counts_multiple_samples(isoform_proportion_samples_new, long_total_reads, isoform_lengths, simulation_dir, 
        long_min_ambiguity, long_max_ambiguity,  long_min_read_length, long_avg_read_length, long_max_read_length, read_name = 'LONG')
    
    file_path = os.path.join(simulation_dir, f"{read_name}{timestamp}_spCorr_{corr_new:.2f}_isoformAbundance_groundTruth.pkl")
    # Save the read_dict as a pickle file
    with open(file_path, 'wb') as f:
        pickle.dump(isoform_proportion_samples_new, f)

    end = time.time()
    interval = (end-start)/60
    print(f"total_time {interval} min")


# Run the main function
if __name__ == "__main__":
    main()
