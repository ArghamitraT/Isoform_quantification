"""
This files is used parse the result statistics from the log file
"""

import re
import pandas as pd
import os
import time
import matplotlib.pyplot as plt
import numpy as np

# Function to find log file paths
def find_log_file_paths(main_dir, experiment_file):
    final_result_dir = os.path.join(main_dir, experiment_file, 'files/output_files/') 
    #final_result_dir = os.path.join(main_dir, experiment_file, 'files/problem/') #(AT)
    log_file_paths = []
    for file in os.listdir(final_result_dir):
        if file.startswith('out_main_') and not file.endswith('.csv'):
            log_file_paths.append(os.path.join(final_result_dir, file))
    return log_file_paths


def parse_log_file(log_file_path, multi_occurrence_vars, single_occurrence_vars):
    with open(log_file_path, 'r') as file:
        lines = file.readlines()

    # Initialize dictionaries to store regular expressions and corresponding data lists
    data = {var: [] for var in multi_occurrence_vars}
    single_occurrence_data = {var: None for var in single_occurrence_vars}

    # Update the regex to handle numbers or the word 'inf'
    regex_patterns_multi = {var: re.compile(f'{var} (-?\\d*\\.?\\d+|\\d+|inf|-inf)') for var in multi_occurrence_vars}
    regex_patterns_single = {var: re.compile(f'{var} (-?\\d*\\.?\\d+|\\d+|inf|-inf|/[^\\s]+)') for var in single_occurrence_vars}

    # Parse the file line by line
    for line in lines:
        for var, regex in regex_patterns_multi.items():
            match = regex.search(line)
            if match:
                value = match.group(1)
                # Convert the value to float or inf as necessary
                data[var].append(float(value) if value not in ['inf', '-inf'] else float(value))

        for var, regex in regex_patterns_single.items():
            match = regex.search(line)
            if match and single_occurrence_data[var] is None:
                value = match.group(1)
                # Check if the value is a file path
                if '/' in value:
                    # Extract the file name from the path
                    single_occurrence_data[var] = value.split('/')[-1]
                else:
                    # Convert the value to float if it's not a file path
                    single_occurrence_data[var] = float(value) if value not in ['inf', '-inf'] else float(value)

    # Create a DataFrame from the parsed data
    df = pd.DataFrame(data)

    return df, single_occurrence_data


# Function to save DataFrame to CSV
def save_dataframe_to_csv(df, log_file_path):
    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    output_file_path = log_file_path + timestamp + '.csv'
    df.to_csv(output_file_path, index=False)
    return output_file_path


def plot_results(df, single_occurrence_data, experiment_file, main_dir, name, columns_to_plot):
    # Create a folder named 'figures'
    figure_dir = os.path.join(main_dir, experiment_file, 'figures')
    os.makedirs(figure_dir, exist_ok=True)

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
        axs[row, col].set_xlabel('EM_loop')
        axs[row, col].set_ylabel(column)
        axs[row, col].set_title(f'{column} vs EM_loop')
        axs[row, col].grid(True)
        axs[row, col].legend()

    # Remove any empty subplots
    for j in range(len(columns_to_plot), n_rows * n_columns):
        fig.delaxes(axs.flatten()[j])

    # Print single-occurrence data on the figure
    single_occurrence_text = "\n".join([f"{key}: {value}" for key, value in single_occurrence_data.items()])
    plt.figtext(0.5, 0.01, single_occurrence_text, ha='center', fontsize=12, bbox={"facecolor": "orange", "alpha": 0.5, "pad": 5})

    # Adjust layout
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Leave space at the bottom for the text

    # Save the figure
    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    plt.savefig(os.path.join(figure_dir, name + timestamp + '.png'))


# Read the file and process lines
def process_file(input_file):
    # Initialize the base counter for iterations
    base_counter = 0

    # Initialize results list
    results = []

    # Initialize lists to store iterations and loss values for plotting
    iterations = []
    losses = []

    # Read the file and process lines
    with open(input_file, 'r') as file:
        for line in file:
            # Check if the line contains an iteration
            if 'GD_Iteration' in line:
                
                # Extract the iteration number
                iteration_number_str = re.findall(r'\d+', line)[0]  # Extract digits
                iteration_number = int(iteration_number_str)  # Convert to integer

                # Update the iteration number based on the base counter
                updated_iteration_number = base_counter + iteration_number

                # Append the updated iteration number to the iterations list for plotting
                iterations.append(updated_iteration_number)


                
            # Check if the line contains a loss value
            if 'GD_Current_Loss' in line:
                # Extract the loss value, including the negative sign if present
                loss_value = float(re.findall(r'-?\d+\.\d+', line)[0])

                # Append the loss value to the losses list for plotting
                losses.append(loss_value)


            # Check if the iteration number has reached 9 to update the base counter
            if 'GD_Iteration 9' in line:
                base_counter += 10

            # Append the processed line to results
            # results.append(line)

    print("File processed successfully.")
    return iterations, losses


def plot_iterations_vs_loss(iterations, losses, experiment_file, main_dir, name):
    figure_dir = os.path.join(main_dir, experiment_file, 'figures')

    plt.figure(figsize=(10, 6))
    plt.plot(iterations, losses, marker='o', linestyle='-')
    plt.xlabel('Iterations')
    plt.ylabel('Loss')
    plt.title('Iterations vs Loss')
    plt.grid(True)
    # Save the figure
    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    plt.savefig(os.path.join(figure_dir, name+'_Gradient' + timestamp + '.png'))
    #plt.show()


# Main function to run the entire process
def main():
    experiment_file = 'exprmnt_2024_08_10__02_05_36'
    main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/'

    # Find all log file paths
    log_file_paths = find_log_file_paths(main_dir, experiment_file)
    if not log_file_paths:
        print("No file found matching the format 'out_main_'")
        return

    # Parse each log file and combine the data
    combined_df = pd.DataFrame()
    for log_file_path in log_file_paths:

        multi_occurrence_vars = ['ELBO_sample_1', 'ELBO_sample_2', 'Convergence_sample_1', 'Convergence_sample_2',
                     'EM_convergence', 'Spearman_corr_theta1_theta2',
                     'Spearman_corr_theta1_alpha', 'Spearman_corr_theta2_alpha', 'alpha_summation ', 'EM_loop']
        
        single_occurrence_vars = ['GD_lr', 'alpha_initial', 'sample_1', 'sample_2']
    
        df, single_occurrence_data = parse_log_file(log_file_path, multi_occurrence_vars, single_occurrence_vars)
        
        # Save the parsed data to a CSV file
        #csv_file_path = save_dataframe_to_csv(df, log_file_path) #(AT)

        # Plot the results
        name = log_file_path.split('/')[-1]
        plot_results(df, single_occurrence_data, experiment_file, main_dir, name, multi_occurrence_vars)

        # (AT)
        iterations, losses = process_file(log_file_path)
        plot_iterations_vs_loss(iterations, losses, experiment_file, main_dir, log_file_path.split('/')[-1])


# Run the main function
if __name__ == "__main__":
    main()
