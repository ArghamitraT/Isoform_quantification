# # import json
# #
# # with open(
# #         '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/CU_courses/Spring_24/CBMF4761/Project/RNA_Splicing/data/PacBio_data_Liz/day0-rep1/'
# #         'day0_rep1.json', 'r') as file:
# #     all_read_dicts = json.load(file)
# #     print()
# #
# #
# # # def dfs(buckt_id, series):
# # #     if buckt_id == 4:
# # #         print(series)
# # #         t = get_total_cost(series)
# # #         s_rate = get_s_rate(series)
# # #         if t<=D:
# # #             if s_rate > pre_rate:
# # #                 target = series
# # #                 pre_rate = s_rate
# # #
# # #         return
# # #     for i in range(3):
# # #         dfs(buckt_id+1, series+[i])
# # #
# # # dfs(0, [])
#
# import random
#
# def randomized_quick_sort(arr, low, high):
#     if low < high:
#         # Randomly select a pivot index and swap with the last element
#         pivot_index = random.randint(low, high)
#         arr[pivot_index], arr[high] = arr[high], arr[pivot_index]
#
#         # Partition the array
#         pivot_index = partition(arr, low, high)
#
#         # Recursively sort the partitions
#         randomized_quick_sort(arr, low, pivot_index - 1)
#         randomized_quick_sort(arr, pivot_index + 1, high)
#
# def partition(arr, low, high):
#     pivot = arr[high]
#     i = low - 1
#     for j in range(low, high):
#         if arr[j] <= pivot:
#             i += 1
#             arr[i], arr[j] = arr[j], arr[i]
#     arr[i + 1], arr[high] = arr[high], arr[i + 1]
#     return i + 1
#
# # Example list of numbers
# numbers = [10, 7, 8, 9, 1, 5]
# print("Original array:", numbers)
#
# # Applying randomized quicksort
# randomized_quick_sort(numbers, 0, len(numbers) - 1)
# print("Sorted array:", numbers)


import types

# Assuming Read is a class or function you have defined somewhere
def Read():
    pass

# Manually adding Read to __main__
import sys
sys.modules['__main__'].Read = Read

# Now try to load your pickle file
import pickle


with open('pkl_files/PacBio_Day0_rep1.pkl', 'rb') as file:
    all_read_dicts = pickle.load(file)
print()