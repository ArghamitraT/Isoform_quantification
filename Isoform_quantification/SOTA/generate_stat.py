import numpy as np
from scipy.stats import spearmanr, pearsonr
import pandas as pd
import datetime
import time
import os
import re
from collections import defaultdict
from scipy.interpolate import UnivariateSpline
#from scipy.integrate import trapz
import sys


def fraction_to_float(fraction_str):
    # Strip the input to remove any leading/trailing whitespace
    fraction_str = fraction_str.strip()

    # Try to convert the string directly to a float (works if there is no fraction)
    try:
        return float(fraction_str)
    except ValueError:
        # Split the string on the '/'
        fraction_parts = fraction_str.split('/')
        if len(fraction_parts) == 2:  # It's a fraction
            numerator = float(fraction_parts[0])
            denominator = float(fraction_parts[1])
            return numerator / denominator
        else:
            raise ValueError(f"Invalid fraction string: {fraction_str}")

def fraction_to_float_gen(value):
    try:
        return float(value)
    except ValueError:
        return None

def split_transcript_id(transcript):
    return transcript.split('.')[0]


def csv_tpm_processing(file_path1, file_path2, suffixes=('_quant', '_truth')):
    # Load the datasets
    our_quant = pd.read_csv(file_path1, sep="\t")
    ground_truth = pd.read_csv(file_path2, sep="\t")

    # Apply the function to the 'transcript_id' column of your DataFrame
    our_quant['target_id'] = our_quant['target_id'].apply(split_transcript_id)
    ground_truth['target_id'] = ground_truth['target_id'].apply(split_transcript_id)

    # Find common isoforms
    common_isoforms = pd.merge(our_quant, ground_truth, on='target_id', suffixes=suffixes)

    return common_isoforms

def spearman_corr_generic(file_path1, file_path2, type):
    common_isoforms = csv_tpm_processing(file_path1, file_path2)

    if 'tpm_quant' in common_isoforms and 'tpm_truth' in common_isoforms:
        common_isoforms['tpm_truth'] = common_isoforms['tpm_truth'].apply(fraction_to_float_gen)
        correlation, p_value = spearmanr(common_isoforms['tpm_quant'], common_isoforms['tpm_truth'])
        pcorrelation, p_value = pearsonr(common_isoforms['tpm_quant'], common_isoforms['tpm_truth'])

        return correlation, pcorrelation


def main():
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    

    """ ## CALL 'spearman_corr_generic' FUNC FIRST ## """
    data_type = 'long'
    sc, pc = spearman_corr_generic(file1, file2, type=data_type)
    #cv = calculate_im_acvc(file1, file2, type=data_type)

    print("SC: ", sc, " PC: ", pc)
    #print(cv)

# Example usage
if __name__ == "__main__":
    main()
