from collections import Counter
import numpy as np
import pickle
import matplotlib.pyplot as plt


def convert_text_to_counter(file_path):
    # Create an empty Counter object
    isoform_counter = Counter()

    # Open the file and read line by line
    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            # Skip the header line
            if i == 0:
                continue

            # Extract the relevant data
            parts = line.strip().split()
            # Assume that the first column is the isoform identifier and the second column is the count
            isoform_id = parts[0].strip().split(',')[0]
            count_fl = int(parts[0].strip().split(',')[1])  # Convert count to integer

            # isoform_id = parts[0]
            # count_fl = int(parts[1])  # Convert count to integer

            # Update the Counter object
            isoform_counter[isoform_id] = count_fl

    return isoform_counter

with open('pkl_files/theta_illumina_day0_rep1.pkl', 'rb') as file:
     theta_sample1 = pickle.load(file)
with open('pkl_files/theta_PacBio_day0_rep1.pkl', 'rb') as file:
    theta_sample2 = pickle.load(file)
# theta_sample1 = Counter({
#     'PB.1.1': 12,
#     'PB.1.2': 25,
#     'PB.1.3': 7,
#     'PB.1.4': 60,
#     'PB.1.5': 2,
#     'PB.1.6': 5,
#     'PB.1.7': 2,
#     'PB.1.8': 2,
#     'PB.1.9': 2,
#     'PB.1.10': 5,
#     'PB.1.11': 56
# })
# theta_sample2 = Counter({
#     'PB.1.100': 12,
#     'PB.1.2': 20,
#     'PB.1.3': 7,
#     'PB.1.4': 60,
#     'PB.1.51': 2,
#     'PB.1.67': 5,
#     'PB.1.7': 2,
#     'PB.1.8': 25,
#     'PB.1.9': 2,
#     'PB.1.10': 5,
#     'PB.1.11': 56
# })
# Find common isoforms between the two samples
#common_isoforms = theta_sample2 & theta_sample1
common_isoforms_dict = {isoform: theta_sample2[isoform] for isoform in theta_sample1 if isoform in theta_sample2}
common_isoforms = Counter(common_isoforms_dict)
# Function to plot the read distribution for a sample
# def plot_read_distribution(sample, title, subplot_number):
#     #plt.subplot(1, 3, subplot_number)
#     # plt.bar(sample.keys(), sample.values())
#     plt.bar(range(len(sample)), sample.values())
#     plt.title(title)
#     #plt.xticks(rotation=90)
#     plt.tight_layout()
#     plt.show()
def plot_read_distribution(sample, title, nameFig):
    # Assuming 'sample' is a dictionary and 'title' and 'nameFig' are defined
    max_key = max(sample, key=sample.get)
    max_value = sample[max_key]

    # Prepare the data
    sorted_items = sorted(sample.items(), key=lambda item: item[1])
    mid_index = len(sorted_items) // 2
    left_half = sorted_items[:mid_index]
    right_half = sorted_items[mid_index:]
    ordered_keys = [k for k, v in left_half] + [max_key] + [k for k, v in reversed(right_half)]
    ordered_values = [sample[k] for k in ordered_keys]

    # Set the figure size and bar width
    fig_width = 20  # You might need to adjust this based on your actual data
    bar_width = fig_width / len(ordered_values)  # Adjust bar width based on number of data points

    # Create the figure and plot the bars
    plt.figure(figsize=(fig_width, 5))
    plt.bar(range(len(ordered_values)), ordered_values, width=bar_width)

    # Set a logarithmic scale if the data varies widely
    plt.yscale('log')

    # Set the y-axis limit to slightly above the max value
    plt.ylim(0.1, max_value * 10)  # Adjust the lower limit if you have zero or negative values

    # Add a title and other plot settings
    plt.title(title)
    plt.tight_layout()

    # Save the figure to a file and show it
    plt.savefig('figures/' + nameFig)
    plt.show()


def plot_normalized_read_distribution(sample, title, nameFig):
    # Normalize the values
    values = np.array(list(sample.values()), dtype=float)
    min_value = np.min(values)
    max_value = np.max(values)
    normalized_values = (values - min_value) / (max_value - min_value)

    # Define 10 bins for normalized data which ranges from 0 to 1
    bins = np.linspace(0, 1, 11)

    # Calculate the count of isoforms in each bin
    hist, bin_edges = np.histogram(normalized_values, bins=bins)

    # Normalize the histogram counts for the y-axis
    normalized_hist = hist / float(hist.sum())

    # Plot the histogram with the bin counts
    fig, ax = plt.subplots(figsize=(10, 6))
    bars = ax.bar(range(len(normalized_hist)), normalized_hist, width=0.8, align='center')

    # Annotate each bar with the actual range of data values
    for bar, bin_count, bin_start, bin_end in zip(bars, hist, bin_edges[:-1], bin_edges[1:]):
        actual_start = int(bin_start * (max_value - min_value) + min_value)
        actual_end = int(bin_end * (max_value - min_value) + min_value)
        annotation_text = f'Reads #: {actual_start}-{actual_end}\nIsoforms #: {bin_count}'

        # Dynamically set the vertical position based on the height of the bar
        vertical_position = bar.get_height() if bar.get_height() > 0 else 0.1
        ax.annotate(annotation_text,
                    xy=(bar.get_x() + bar.get_width() / 2, vertical_position),
                    xytext=(0, 5),  # Adjust text position for visibility
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8, rotation=45, color='blue')

    # Set the title and labels
    plt.yscale('log')  # Log scale for y-axis
    plt.title(title)
    plt.xlabel('Normalized Value Bins')
    plt.ylabel('Isoforms # (normalized)')

    # Adjust the x-tick labels to match the bin edges and prevent overlap
    ax.set_xticks(range(len(bin_edges) - 1))
    ax.set_xticklabels([f"{int(edge):d}" for edge in bin_edges[:-1]], rotation=45)

    # Save the figure to a file and show it
    plt.tight_layout()
    plt.savefig(nameFig)
    plt.show()



# Plot read distribution for theta_sample1
# plot_read_distribution(theta_sample1, 'Read Distribution for Sample 1', 'theta_sample1.png')
#plot_normalized_read_distribution(theta_sample1, 'Read Distribution for Sample 1', 'theta_sample1.png')

# Plot read distribution for theta_sample2
# plot_normalized_read_distribution(theta_sample2, 'Read Distribution for Sample 2', 'theta_sample2.png')

# Plot overlapping read distribution for common isoforms
plot_normalized_read_distribution(common_isoforms, 'Overlapping Read Distribution', 'theta_common.png')

# Display the plot

file_path = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/CU_courses/Spring_24/CBMF4761/Project/RNA_Splicing/data/PacBio_data_Liz/day0-rep1/4-Collapse/'
# Replace 'path_to_your_text_file.txt' with the actual path to your text file
# counter_from_text_file = convert_text_to_counter(file_path+'collapsed_transcripts.flnc_count.txt')
