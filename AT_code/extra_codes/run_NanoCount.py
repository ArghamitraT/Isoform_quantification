from NanoCount_AT import NanoCount
from NanoCount.common import jhelp, head
import os

# input_file = os.path.join(os.getcwd(), '../../data/NanoCount_data/ERR4352441_aligned_sorted_reads.bam')
#output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ERR4352441_NC_output.tsv')

# input_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF011BFA_aligned_reads.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF011BFA_NC_output.tsv')

# input_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF417ZXW_aligned_reads.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF417ZXW_NC_output.tsv')

# input_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF646VGI.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF646VGI_aligned_reads_NC_output.tsv')

# input_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF097KLY.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF097KLY_aligned_reads_NC_output.tsv')


input_file = "pkl_files/Day0_2_illumina.bam"
output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/Illu_day0rep2_NC.tsv')

NanoCount (alignment_file=input_file, downsampled_files= "pkl_files/Day0_1_illumina.bam", count_file=output_file)


print("I am done")