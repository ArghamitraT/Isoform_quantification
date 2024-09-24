from scipy.stats import spearmanr
import pickle


main_directory = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/simulation/round8/'
theta_file_pairs = [ ['ds_100_num1_aln_01_short_theta.pkl', 'ds_10_num1_aln_01_long_theta.pkl'],
                    ['ds_100_num1_aln_02_short_theta.pkl', 'ds_10_num1_aln_02_long_theta.pkl'],
                    ['ds_100_num1_aln_11_short_theta.pkl', 'ds_10_num1_aln_11_long_theta.pkl'],
                    ['ds_100_num1_aln_12_short_theta.pkl', 'ds_10_num1_aln_12_long_theta.pkl']]


for file_pairs in theta_file_pairs:
    theta_S1_name = file_pairs[0]
    theta_S2_name = file_pairs[1]

    with open(main_directory+theta_S1_name, 'rb') as f:
        theta_S1 = pickle.load(f)
    with open(main_directory+theta_S2_name, 'rb') as f:
        theta_S2 = pickle.load(f)



    all_keys = set(theta_S1.keys()).union(set(theta_S2.keys()))
    values_S1 = [theta_S1.get(key, 0) for key in all_keys]
    values_S2 = [theta_S2.get(key, 0) for key in all_keys]
    corr, _ = spearmanr(values_S1, values_S2)
    print(f"file1 {theta_S1_name} file2 {theta_S2_name}")
    print(f"SP corr {corr}")
    # Get the keys (isoforms) for each counter
    keys_S1 = set(theta_S1.keys())
    keys_S2 = set(theta_S2.keys())

    # Find the intersection of the two sets (common isoforms)
    common_isoforms = keys_S1.intersection(keys_S2)
    print(f"common_isoform_num {len(common_isoforms)}")
print()