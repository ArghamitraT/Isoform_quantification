from EM_MAP_GD import Expec_Max
import os

# input file names
#long_read_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF646VGI.bam')
# long_read_file = os.path.join(os.getcwd(), '../../data/NanoCount_data/ERR4352441_aligned_sorted_reads.bam')
# short_read_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF097KLY.bam')
# #output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF646VGI_LR_ENCFF097KLY_ShR__NC_output.tsv')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/output')


# SIRV_E0 = os.path.join(os.getcwd(), '../../data/SIRV/SIRV_Shree/alignments/aln_E0.bam')
# SIRV_E2 = os.path.join(os.getcwd(),'../../data/SIRV/SIRV_Shree/alignments/aln_E2.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/SIRV_output')

main_folder = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/CU_courses/Spring_24/CBMF4761/Project/RNA_Splicing/data/PacBio_data_Liz/'
SIRV_E0 = main_folder + 'day0-rep1/3-ClusterMap/mapped.bam'
#SIRV_E0 = main_folder + 'day0-rep1/Day0_1_illumina.bam'
SIRV_E2 = main_folder + 'day0-rep2/Day0_2_illumina.bam'
output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/Illumina_output')

# call the main class
# Expec_Max (file_names=[long_read_file, short_read_file], count_file=output_file)
Expec_Max (file_names=[SIRV_E0, SIRV_E2], count_file=output_file)


print("I am done")