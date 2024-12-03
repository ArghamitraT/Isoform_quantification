"""
This script plots the result stats from the pickle files
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
import time
import pickle

def plot_GDloss_stats(final_result_dir, main_dir, experiment_file, fig_generic_name):
    print()
    # Step 1: Read all files starting with 'allEMstats'
    emstat_files = [f for f in os.listdir(final_result_dir) if f.startswith('allGDloss') and f.endswith('.pkl')]
    # Create a folder named 'figures'
    figure_dir = os.path.join(main_dir, experiment_file, 'figures')
    os.makedirs(figure_dir, exist_ok=True)
    # Regular expression to extract GDlr, AlphaInitial, and EMround from the filenames
    # file_pattern = re.compile(r'file1_(.*?)_file2_(.*?)_GDlr_(0\.\d+)_AlphaInitial_(\d+(\.\d+)?)_EMround_(\d+)_token_(\d+)_')
    file_pattern = re.compile(r'file1_(.*?)_GDlr_(0\.\d+)_AlphaInitial_(\d+(\.\d+)?)_EMround_(\d+)_token_(\d+)_')

    for file in emstat_files:
        file_path = os.path.join(final_result_dir, file)
        with open(file_path, 'rb') as file:
            data = pickle.load(file)
        
        # Extract GDlr, AlphaInitial, and EMround from filename
        match = file_pattern.search(file.name)
        if match:
            # filename1, filename2, GDlr, AlphaInitial, _, EMround, token = match.groups()  # Extract the desired groups (AT)
            filename2 = 'NA'
            filename1, GDlr, AlphaInitial, _, EMround, token = match.groups()  # Extract the desired groups
        else:
            print(f"Could not extract metadata from file name: {file}")
            continue
        # Convert GDlr and AlphaInitial to scientific notation
        GDlr = "{:.1e}".format(float(GDlr))
        AlphaInitial= "{:.1e}".format(float(AlphaInitial))
        
        plt.figure(figsize=(10, 6))
        plt.plot(data.index, data['GD_Loss'], marker='o', linestyle='-', color='b')
        plt.xlabel('Iterations')
        plt.ylabel('NegativeLogLikelihood')
        plt.title('Iterations vs GDLoss')
        plt.grid(True)
        # Print the single occurrence data at the bottom of the figure
        single_occurrence_data = {
            "Filename 1": filename1,
            "Filename 2": filename2,
            "GDlr": GDlr,
            "AlphaInitial": AlphaInitial,
            "EMround": EMround
        }
        single_occurrence_text = "\n".join([f"{key}: {value}" for key, value in single_occurrence_data.items()])
        
        # Print the single-occurrence text at the bottom of the figure
        plt.figtext(0.5, 0.01, single_occurrence_text, ha='center', fontsize=12, bbox={"facecolor": "orange", "alpha": 0.5, "pad": 5})
                
        # Adjust layout
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Leave space at the bottom for the text

        # Save the figure
        timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
        fig_name = f"{fig_generic_name}_file1_{filename1}_file2_{filename2}_GDlr_{GDlr}_AlphaInitial_{AlphaInitial}_EMround_{EMround}_token_{token}"
        plt.savefig(os.path.join(figure_dir, fig_name + timestamp + '.png'))
        plt.close()  # Close the figure after saving
    
    print("GD loss Plotting and saving completed.")



def plot_result_stats(final_result_dir, main_dir, experiment_file, columns_to_plot, fig_generic_name):

    # Step 1: Read all files starting with 'allEMstats'
    emstat_files = [f for f in os.listdir(final_result_dir) if f.startswith('allEMstats') and f.endswith('.csv')]

    # Create a folder named 'figures'
    figure_dir = os.path.join(main_dir, experiment_file, 'figures')
    os.makedirs(figure_dir, exist_ok=True)
    # Regular expression to extract GDlr, AlphaInitial, and EMround from the filenames
    # file_pattern = re.compile(r'file1_(.*?)_file2_(.*?)_GDlr_(0\.\d+)_AlphaInitial_(\d+(\.\d+)?)_EMround_(\d+)_token_(\d+)_')
    file_pattern = re.compile(r'file1_(.*?)_GDlr_(0\.\d+)_AlphaInitial_(\d+(\.\d+)?)_EMround_(\d+)_token_(\d+)_')


    # Step 3: Load all files, extract metadata, and plot specified columns
    for file in emstat_files:
        file_path = os.path.join(final_result_dir, file)
        
        # Extract GDlr, AlphaInitial, and EMround from filename
        match = file_pattern.search(file)
        if match:
            #filename1, filename2, GDlr, AlphaInitial, _, EMround, token = match.groups()  # Extract the desired groups (AT)
            filename2 = 'NA'
            filename1, GDlr, AlphaInitial, _, EMround, token = match.groups()  # Extract the desired groups
        else:
            print(f"Could not extract metadata from file name: {file}")
            continue
        # Convert GDlr and AlphaInitial to scientific notation
        GDlr = "{:.1e}".format(float(GDlr))
        AlphaInitial= "{:.1e}".format(float(AlphaInitial))
        
        # Read the CSV file
        df = pd.read_csv(file_path)
        
        # Determine the number of rows and columns for the subplot grid
        n_columns = 2  # Number of columns in the grid
        n_rows = (len(columns_to_plot) + 1) // n_columns  # Number of rows needed

        # Create a figure with subplots
        fig, axs = plt.subplots(n_rows, n_columns, figsize=(15, 5 * n_rows))

        # Plot each column against EM_loop in a subplot
        for i, column in enumerate(columns_to_plot):
            row = i // n_columns
            col = i % n_columns
            
            # Mask for finite and infinite values
            finite_mask = np.isfinite(df[column])
            inf_mask = np.isinf(df[column])

            # Plot finite values
            axs[row, col].plot(df['EM_loop'][finite_mask], df[column][finite_mask], marker='o', label=f'{column} (finite)')
            
            # Plot inf values with a distinct marker and color
            if inf_mask.any():
                axs[row, col].plot(df['EM_loop'][inf_mask], df[column][inf_mask], 'rx', label=f'{column} (inf)')
            
            # Set labels and title
            axs[row, col].set_xlabel('EM_loop', fontsize=16)
            axs[row, col].set_ylabel(column, fontsize=16)
            axs[row, col].set_title(f'{column} vs EM_loop', fontsize=18)
            axs[row, col].grid(True)
            #axs[row, col].legend()
            # Adjust tick parameters for the grid (both x and y axes)
            axs[row, col].tick_params(axis='both', which='major', labelsize=14)


        # Remove any empty subplots
        for j in range(len(columns_to_plot), n_rows * n_columns):
            fig.delaxes(axs.flatten()[j])
        
        # Print the single occurrence data at the bottom of the figure
        single_occurrence_data = {
            "Filename 1": filename1,
            "Filename 2": filename2,
            "GDlr": GDlr,
            "AlphaInitial": AlphaInitial,
            "EMround": EMround
        }
        single_occurrence_text = "\n".join([f"{key}: {value}" for key, value in single_occurrence_data.items()])
        
        # Print the single-occurrence text at the bottom of the figure
        plt.figtext(0.5, 0.01, single_occurrence_text, ha='center', fontsize=12, bbox={"facecolor": "orange", "alpha": 0.5, "pad": 5})
                
        # Adjust layout
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Leave space at the bottom for the text

        # Save the figure
        timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
        fig_name = f"{fig_generic_name}_file1_{filename1}_file2_{filename2}_GDlr_{GDlr}_AlphaInitial_{AlphaInitial}_EMround_{EMround}_token_{token}"
        plt.savefig(os.path.join(figure_dir, fig_name + timestamp + '.png'))
        plt.close()  # Close the figure after saving


    print("result stats Plotting and saving completed.")

# Main function to run the entire process
def main():
    experiment_file = 'exprmnt_2024_11_14__13_21_32'
    # "Different experiment setup, 1: for 1 sample, 2 for merged, 4 for multisample, 5 for merged multisample"
    experiment = 4


    main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/'
    
    if experiment == 1 or experiment==2:
        columns_to_plot = ['ELBO_sample_1', 'Convergence_sample_1', 'EM_convergence', 
                            'Spearman_corr_theta1_alpha', 'Alpha_summation', 'EM_loop']
    elif experiment == 4 :
        columns_to_plot = ['ELBO_sample_1', 'ELBO_sample_2', 'Convergence_sample_1', 'Convergence_sample_2',
                        'EM_convergence', 'Spearman_corr_theta1_theta2',
                        'Spearman_corr_theta1_alpha', 'Spearman_corr_theta2_alpha', 'Alpha_summation', 'EM_loop']

    
    final_result_dir = os.path.join(main_dir, experiment_file, 'weights/') 
    
    plot_result_stats(final_result_dir, main_dir, experiment_file, columns_to_plot, fig_generic_name="result_stats")
    plot_GDloss_stats(final_result_dir, main_dir, experiment_file, fig_generic_name="GDloss_stats")


# Run the main function
if __name__ == "__main__":
    main()


"""
great. now we generate read count from these isoforms. first short reads, for each sample
1. number of reads would be proportionate isoform proportions we have in sample dict. so if a isoform A is present in 80% and another B is 20% then a significant number of reads would be generated from isoform A then isoform B
2. we are generating short read counts, here how it is going to be:
a. generate a read hash, and for each read hash assign isoforms
b. as they are
"""