from EM_VI_GD import Expec_Max
import os
import argparse


main_folder = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/data/'
main_folder = '/gpfs/commons/home/atalukder/RNA_Splicing/data/'
sample1 = main_folder + 'SIRV/SIRV_Shree/alignments/aln_E0.bam'
sample2 = main_folder + 'SIRV/SIRV_Shree/alignments/aln_E2.bam'
output_file_default = os.path.join(os.getcwd(), '../../files/results/single_run/output_PacIllu_VIGD')
# main_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/'
# sample1 = main_folder + 'ds_5_aln_01_long.bam'
# sample2 = main_folder + 'aln_02_short.bam'
# output_file_default = os.path.join(os.getcwd(), '../../files/results/single_run/output_PacIllu_VIGD')

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description="Process BAM files and output results.")
parser.add_argument("--output_path", type=str, default=output_file_default,
                    help="Path for the output file. Default is '../../files/results/single_run/'.")

# Parse the arguments
args = parser.parse_args()
# Assign the user-provided or default path to the output_file variable
output_file = args.output_path


# call the main class
Expec_Max (file_names=[sample1, sample2], count_file=output_file)
# Expec_Max (file_names=[sample2], count_file=output_file)


print("#########END###########")