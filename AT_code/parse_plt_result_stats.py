import re
import pandas as pd
import os
import time
import matplotlib.pyplot as plt

# Function to find log file paths
def find_log_file_paths(main_dir, experiment_file):
    final_result_dir = os.path.join(main_dir, experiment_file, 'files/output_files/')
    log_file_paths = []
    for file in os.listdir(final_result_dir):
        if file.startswith('out_main_') and not file.endswith('.csv'):
            log_file_paths.append(os.path.join(final_result_dir, file))
    return log_file_paths

# Function to parse a single log file
def parse_log_file(log_file_path):
    with open(log_file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables to store parsed data
    em_loops, elbo_sample_1, elbo_sample_2, convergence_sample_1, convergence_sample_2, em_convergence = [], [], [], [], [], []

    # Regular expressions to match the required lines
    em_loop_re = re.compile(r'EM_loop (\d+)')
    elbo_sample_1_re = re.compile(r'ELBO_sample_1 ([\d.-]+)')
    elbo_sample_2_re = re.compile(r'ELBO_sample_2 ([\d.-]+)')
    convergence_sample_1_re = re.compile(r'Convergence_sample_1 ([\d.-]+)')
    convergence_sample_2_re = re.compile(r'Convergence_sample_2 ([\d.-]+)')
    em_convergence_re = re.compile(r'EM_convergence ([\d.-]+)')

    # Parse the file line by line
    for line in lines:
        em_loop_match = em_loop_re.search(line)
        if em_loop_match:
            em_loops.append(int(em_loop_match.group(1)))

        elbo_sample_1_match = elbo_sample_1_re.search(line)
        if elbo_sample_1_match:
            elbo_sample_1.append(float(elbo_sample_1_match.group(1)))

        elbo_sample_2_match = elbo_sample_2_re.search(line)
        if elbo_sample_2_match:
            elbo_sample_2.append(float(elbo_sample_2_match.group(1)))

        convergence_sample_1_match = convergence_sample_1_re.search(line)
        if convergence_sample_1_match:
            convergence_sample_1.append(float(convergence_sample_1_match.group(1)))

        convergence_sample_2_match = convergence_sample_2_re.search(line)
        if convergence_sample_2_match:
            convergence_sample_2.append(float(convergence_sample_2_match.group(1)))

        em_convergence_match = em_convergence_re.search(line)
        if em_convergence_match:
            em_convergence.append(float(em_convergence_match.group(1)))

    # Create a DataFrame from the parsed data
    data = {
        'EM_loop': em_loops,
        'ELBO_sample_1': elbo_sample_1,
        'ELBO_sample_2': elbo_sample_2, # (AT)
        'Convergence_sample_1': convergence_sample_1,
        'Convergence_sample_2': convergence_sample_2, # (AT)
        'EM_convergence': em_convergence
    }

    return pd.DataFrame(data)

# Function to save DataFrame to CSV
def save_dataframe_to_csv(df, log_file_path):
    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    output_file_path = log_file_path + timestamp + '.csv'
    df.to_csv(output_file_path, index=False)
    return output_file_path

# Function to plot results
def plot_results(df, experiment_file, main_dir, name):
    # Create a folder named 'figures'
    figure_dir = os.path.join(main_dir, experiment_file, 'figures')
    os.makedirs(figure_dir, exist_ok=True)

    # List of columns to plot
    # (AT)
    columns_to_plot = ['ELBO_sample_1', 'ELBO_sample_2', 'Convergence_sample_1', 'Convergence_sample_2', 'EM_convergence']
    #columns_to_plot = ['ELBO_sample_1', 'Convergence_sample_1', 'EM_convergence']

    # Create a figure with subplots
    fig, axs = plt.subplots(len(columns_to_plot), 1, figsize=(10, 15))

    # Plot each column against EM_loop in a subplot
    for i, column in enumerate(columns_to_plot):
        axs[i].plot(df['EM_loop'], df[column], marker='o')
        axs[i].set_xlabel('EM_loop')
        axs[i].set_ylabel(column)
        axs[i].set_title(f'{column} vs EM_loop')
        axs[i].grid(True)

    # Adjust layout
    plt.tight_layout()

    # Save the figure
    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    plt.savefig(os.path.join(figure_dir, name + timestamp + '.png'))

    # Show the plot
    plt.show()

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
                parts = line.split(' ')
                iteration_number = int(parts[1])

                # Update the iteration number based on the base counter
                parts[1] = str(base_counter + iteration_number)
                # line = ' '.join(parts)

                # Append the updated iteration number to the iterations list for plotting
                iterations.append(base_counter + iteration_number)

            # Check if the line contains a loss value
            if 'GD_Current_Loss' in line:
                parts = line.split(' = ')
                loss_value = float(parts[1])
                # Append the loss value to the losses list for plotting
                losses.append(loss_value)

            # Check if the iteration number has reached 9 to update the base counter
            if 'GD_Iteration 9' in line:
                base_counter += 10

            # Append the processed line to results
            # results.append(line)

    # Write the modified lines to a new file
    # with open(output_file, 'w') as file:
    #     file.writelines(results)

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
    plt.show()


# Main function to run the entire process
def main():
    experiment_file = 'exprmnt_2024_07_16__12_46_32'
    main_dir = '/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/RNA_Splicing/files/results/'

    # Find all log file paths
    log_file_paths = find_log_file_paths(main_dir, experiment_file)
    if not log_file_paths:
        print("No file found matching the format 'out_main_'")
        return

    # Parse each log file and combine the data
    combined_df = pd.DataFrame()
    for log_file_path in log_file_paths:
        iterations, losses = process_file(log_file_path)
        # (AT)
        # plot_iterations_vs_loss(iterations, losses, experiment_file, main_dir, log_file_path.split('/')[-1])

        df = parse_log_file(log_file_path)
        #combined_df = pd.concat([combined_df, df], ignore_index=True)

        # Save the parsed data to a CSV file
        csv_file_path = save_dataframe_to_csv(df, log_file_path)

        # Plot the results
        plot_results(df, experiment_file, main_dir, log_file_path.split('/')[-1])


# Run the main function
if __name__ == "__main__":
    main()
