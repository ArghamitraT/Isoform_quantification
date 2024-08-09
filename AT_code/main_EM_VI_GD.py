from EM_VI_GD import Expec_Max
import os
import argparse


##main_folder = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/data/'

main_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/'
output_file_default = os.path.join(os.getcwd(), '../../files/results/single_run/files/output_PacIllu_VIGD_00000')

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description="Process BAM files and output results.")
parser.add_argument("--output_path", type=str, default=output_file_default,
                    help="Path for the output file. Default is '../../files/results/single_run/files/output_files/output_generic_VIGD_00000'.")
parser.add_argument("--sample1", type=str, default='ds_5_aln_02_long.bam', help="Sample1 (LR) file name.")
parser.add_argument("--sample2", type=str, default='aln_01_short.bam', help="Sample2 (SR) file name.")
parser.add_argument("--GD_lr", type=float, default=0.01, help="Learning rate for dirichlet gradient descent.")
parser.add_argument("--alpha_initial", type=float, default=10000, help="The fixed sum value of alpha")


# Parse the arguments
args = parser.parse_args()
# Assign the user-provided or default path to the output_file variable
output_file = args.output_path
sample1 = main_folder + args.sample1
sample2 = main_folder + args.sample2
GD_lr = args.GD_lr
alpha_initial = args.alpha_initial


## (AT) SIRV
# main_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/SIRV/'
# sample1 = main_folder + 'aln_E0.bam'
# sample2 = main_folder + 'aln_E2.bam'
# output_file = os.path.join(os.getcwd(), '../../files/results/single_run/files/output_files/output_SIRV_VIGD_00000')

# Print all the parameters
print("Output_file_path", output_file)
print("GD_lr", GD_lr)
print("alpha_initial", alpha_initial)

# call the main class ## (AT)
# Expec_Max (file_names=[sample2, sample2], count_file=output_file, GD_lr=GD_lr, alpha_initial=alpha_initial)
Expec_Max (file_names=[sample1, sample2], count_file=output_file, GD_lr=GD_lr, alpha_initial=alpha_initial)

# takes 100min/EM iteration

print("#########END###########")#