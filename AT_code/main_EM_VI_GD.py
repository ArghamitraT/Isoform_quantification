from EM_VI_GD import Expec_Max
import os
import argparse


##main_folder = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/data/'

main_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/'
output_file_default = os.path.join(os.getcwd(), '../../files/results/single_run/output_PacIllu_VIGD')

# Set up argparse to handle command-line arguments
parser = argparse.ArgumentParser(description="Process BAM files and output results.")
parser.add_argument("--output_path", type=str, default=output_file_default,
                    help="Path for the output file. Default is '../../files/results/single_run/'.")
parser.add_argument("--sample1", type=str, help="Sample1 file name.")
parser.add_argument("--sample2", type=str, help="Sample2 file name.")

## (AT)
# parser.add_argument("--sample1", type=str, default='ds_5_aln_01_long.bam', help="Sample1 file name.")
# parser.add_argument("--sample2", type=str, default='aln_02_short.bam', help="Sample2 file name.")

# Parse the arguments
args = parser.parse_args()
# Assign the user-provided or default path to the output_file variable
output_file = args.output_path
sample1 = main_folder + args.sample1
sample2 = main_folder + args.sample2


## (AT) SIRV
# main_folder = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/SIRV/'
# sample1 = main_folder + 'aln_E0.bam'
# sample2 = main_folder + 'aln_E2.bam'
# output_file = os.path.join(os.getcwd(), '../../files/results/single_run/output_SIRV_VIGD')

# call the main class
Expec_Max (file_names=[sample1, sample2], count_file=output_file)
# Expec_Max (file_names=[sample2], count_file=output_file)


print("#########END###########")#