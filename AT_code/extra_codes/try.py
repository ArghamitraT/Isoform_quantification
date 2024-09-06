import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import pickle
import os


# script_dir='/gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code/output_dummy/output_dummy/out_main_000000__2024_07_13__22_52_01.txt'
# data_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmntDummy_2024_00_00__00_00_00/files/output_files/'
# os.system(f"cp {script_dir} {data_dir}")



pickle_file1_path = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/single_run/weights/allWeights_00000_2024_7_31_17_38_19.pkl'

with open(pickle_file1_path, 'rb') as f:
        our_quant = pickle.load(f)
print()

"""


/gpfs/commons/home/atalukder/miniconda3/envs/NanoCount_5/bin/python /gpfs/commons/home/atalukder/RNA_Splicing/code/AT_code/main_EM_VI_GD.py > /gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmntDummy_2024_00_00__00_00_00/files/out_main_000000__2024_07_13__22_52_00.txt 2>&1



"""


# # from NanoCount_AT import NanoCount
# # from NanoCount.common import jhelp, head
# import os

# # input_file = os.path.join(os.getcwd(), '../../data/NanoCount_data/ERR4352441_aligned_sorted_reads.bam')
# # output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ERR4352441_NC_output.tsv')

# # long_read_file = os.path.join(os.getcwd(), '../../data/PacBio_data_longWshort/ENCFF011BFA_aligned_reads.bam')
# # short_read_file = os.path.join(os.getcwd(), '../../data/NanoCount_data/ERR4352441_aligned_sorted_reads.bam')
# # output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF011BFA_ERR4352441_NC_output.tsv')

# # short_read_file = os.path.join(os.getcwd(), '../../data/PacBio_data_longWshort/ENCFF417ZXW_aligned_reads.bam')
# # output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF011BFA_lr_ENCFF417ZXW_shr_NC_output.tsv')

# long_read_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF646VGI.bam')
# short_read_file = os.path.join(os.getcwd(),'../../data/PacBio_data_longWshort/ENCFF097KLY.bam')
# output_file = os.path.join(os.getcwd(), '../../files/NanoCount_output/ENCFF646VGI_LR_ENCFF097KLY_ShR__NC_output.tsv')


# NanoCount (short_read_file=short_read_file, long_read_file=long_read_file, count_file=output_file)


# print("I am done")


# # # import json
# # #
# # # with open(
# # #         '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/CU_courses/Spring_24/CBMF4761/Project/RNA_Splicing/data/PacBio_data_Liz/day0-rep1/'
# # #         'day0_rep1.json', 'r') as file:
# # #     all_read_dicts = json.load(file)
# # #     print()
# # #
# # #
# # # # def dfs(buckt_id, series):
# # # #     if buckt_id == 4:
# # # #         print(series)
# # # #         t = get_total_cost(series)
# # # #         s_rate = get_s_rate(series)
# # # #         if t<=D:
# # # #             if s_rate > pre_rate:
# # # #                 target = series
# # # #                 pre_rate = s_rate
# # # #
# # # #         return
# # # #     for i in range(3):
# # # #         dfs(buckt_id+1, series+[i])
# # # #
# # # # dfs(0, [])
# #
# # import random
# #
# # def randomized_quick_sort(arr, low, high):
# #     if low < high:
# #         # Randomly select a pivot index and swap with the last element
# #         pivot_index = random.randint(low, high)
# #         arr[pivot_index], arr[high] = arr[high], arr[pivot_index]
# #
# #         # Partition the array
# #         pivot_index = partition(arr, low, high)
# #
# #         # Recursively sort the partitions
# #         randomized_quick_sort(arr, low, pivot_index - 1)
# #         randomized_quick_sort(arr, pivot_index + 1, high)
# #
# # def partition(arr, low, high):
# #     pivot = arr[high]
# #     i = low - 1
# #     for j in range(low, high):
# #         if arr[j] <= pivot:
# #             i += 1
# #             arr[i], arr[j] = arr[j], arr[i]
# #     arr[i + 1], arr[high] = arr[high], arr[i + 1]
# #     return i + 1
# #
# # # Example list of numbers
# # numbers = [10, 7, 8, 9, 1, 5]
# # print("Original array:", numbers)
# #
# # # Applying randomized quicksort
# # randomized_quick_sort(numbers, 0, len(numbers) - 1)
# # print("Sorted array:", numbers)


# import types

# # Assuming Read is a class or function you have defined somewhere
# def Read():
#     pass

# # Manually adding Read to __main__
# import sys
# sys.modules['__main__'].Read = Read

# # Now try to load your pickle file
# import pickle


# with open('pkl_files/PacBio_Day0_rep1.pkl', 'rb') as file:
#     all_read_dicts = pickle.load(file)
# print()