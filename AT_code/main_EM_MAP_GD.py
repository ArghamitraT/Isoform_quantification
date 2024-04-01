from EM_MAP_GD import Expec_Max
import os

# input file names
#long_read_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF646VGI.bam')
long_read_file = os.path.join(os.getcwd(), '../../data/NanoCount_data/ERR4352441_aligned_sorted_reads.bam')
short_read_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF097KLY.bam')
output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF646VGI_LR_ENCFF097KLY_ShR__NC_output.tsv')

# call the main class
Expec_Max (file_names=[long_read_file, short_read_file], count_file=output_file)


print("I am done")