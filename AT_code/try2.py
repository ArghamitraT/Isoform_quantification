import pickle
import os
from scipy.stats import spearmanr, rankdata


simulation_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/round2/'
timestamp = '_2024_09_09__16_12_04'     # (AT) comment
file_name = 'SHORT_spCorr_0.89_isoformAbundance_groundTruth'
dir = os.path.join(simulation_dir, f"{file_name}{timestamp}.pkl")
with open(dir, 'rb') as file:
    short_read_theta = pickle.load(file)

file_name = 'LONG_spCorr_0.89_isoformAbundance_groundTruth'
dir = os.path.join(simulation_dir, f"{file_name}{timestamp}.pkl")
with open(dir, 'rb') as file:
    long_read_theta = pickle.load(file)

# Extract data from sample_1 and sample_2
sample_1_values = []
sample_2_values = []

# Make sure to compare only the common keys between sample_1 and sample_2
for key in long_read_theta['sample_1'].keys():
    if key in short_read_theta['sample_1']:
        sample_1_values.append(long_read_theta['sample_1'][key])
        sample_2_values.append(short_read_theta['sample_1'][key])

# Calculate the Spearman correlation
corr, p_value = spearmanr(sample_1_values, sample_2_values)

print(f"Spearman correlation: {corr}")
print(f"P-value: {p_value}")

print()


