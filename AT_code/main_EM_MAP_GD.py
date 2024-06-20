from EM_VI_GD import Expec_Max
import os

main_folder = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/data/'
SIRV_E0 = main_folder + 'SIRV/SIRV_Shree/alignments/aln_E0.bam'
SIRV_E2 = main_folder + 'SIRV/SIRV_Shree/alignments/aln_E2.bam'
output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/output_SIRV_VIGD')


# call the main class
Expec_Max (file_names=[SIRV_E0, SIRV_E2], count_file=output_file)


print("I am done")