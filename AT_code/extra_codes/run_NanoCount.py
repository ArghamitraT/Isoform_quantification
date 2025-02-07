from NanoCount import NanoCount
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


input_file = "/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/pklfiles/new_minimap_file_originalparsing/ds_100_num1_aln_21_long.bam"
output_file = os.path.join(os.getcwd(), '../../files/results/exprmntSingleRun_2024_00_00__00_00_00/files/output_files/nancount')

NanoCount (alignment_file=input_file, count_file=output_file)


print("I am done")