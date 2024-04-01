#from NanoCount_AT import NanoCount
from EM_MAP_GD import Expec_Max
# from NanoCount.common import jhelp, head
import os

# input_file = os.path.join(os.getcwd(), '../../data/NanoCount_data/ERR4352441_aligned_sorted_reads.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ERR4352441_NC_output.tsv')

# long_read_file = os.path.join(os.getcwd(), '../../data/PacBio_data_longWshort/ENCFF011BFA_aligned_reads.bam')
# short_read_file = os.path.join(os.getcwd(), '../../data/NanoCount_data/ERR4352441_aligned_sorted_reads.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF011BFA_ERR4352441_NC_output.tsv')

# short_read_file = os.path.join(os.getcwd(), '../../data/PacBio_data_longWshort/ENCFF417ZXW_aligned_reads.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF011BFA_lr_ENCFF417ZXW_shr_NC_output.tsv')

#long_read_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF646VGI.bam')
long_read_file = os.path.join(os.getcwd(), '../../data/NanoCount_data/ERR4352441_aligned_sorted_reads.bam')
short_read_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF097KLY.bam')
output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF646VGI_LR_ENCFF097KLY_ShR__NC_output.tsv')


# Expec_Max (short_read_file=short_read_file, long_read_file=long_read_file, count_file=output_file)
Expec_Max (file_names=[long_read_file, short_read_file], count_file=output_file)


print("I am done")