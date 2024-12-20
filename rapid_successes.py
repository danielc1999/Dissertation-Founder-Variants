#!/usr/bin/python3

# This script is used to plot graphs for the accuracy and power metrics of RaPID
# with different window sizes, number of runs and number of successes.

import matplotlib.pyplot as plt
import numpy as np

# Read only the lines which contain the numerical data to plot the graph.
def read_lines(file_path):
    '''
    Read specific lines from a file and extract numerical data.

    :param file_path: The path to the input file containing IBD benchmark result
    
    :return: A list containing numerical data from specific lines in the file
    '''
    lines_to_read = [3, 6, 9, 12, 15, 18, 21]
    with open(file_path, "r") as file:
        lines = [
            line.strip().split("\t")[1:] 
            for i, line in enumerate(file) if i + 1 in lines_to_read
        ]
    return lines

# Define the value of the genotype error rate. This value can be changed.
error_rate = "0.01"

# Define the list of threshold values
thresholds = [0.9]

# Define the file paths for the four data files
file_paths = [
    f"./{error_rate}_error/w3s4.max.IBD_BM_Result.txt",
    f"./{error_rate}_error/w3s7.max.IBD_BM_Result.txt",
    f"./{error_rate}_error/w5s4.max.IBD_BM_Result.txt",
    f"./{error_rate}_error/w5s7.max.IBD_BM_Result.txt",
    f"./{error_rate}_error/w3s2.max.IBD_BM_Result.txt",
    f"./{error_rate}_error/w5s2.max.IBD_BM_Result.txt"
]

# Headings for the 7 accuracy and power metric subplots
headings = [
    "Accuracy", "Length Accuracy", "Length Discrepancy", "Recall", 
    "Power", "Accumulative Recall", "Accumulative Power"
]

# Simplified names for plot legend
legend_names = ["w3s4", "w3s7", "w5s4", "w5s7", "w3s2", "w5s2"]

# Corresponding markers and linestyles, one for each tool
markers = ["o", "^", "s", "x", "*", "v"]
linestyles = ["-", "--", "-.", ":", "--", "dotted"]

for threshold in thresholds:
    # Set up the subplots in a single figure, arranged in a 2x4 grid
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    
    # Flatten the 2D array of axes into a 1D array for easier indexing
    axes = axes.flatten()

    # Create a plot for each heading
    for i, heading in enumerate(headings):
        for j, file_path in enumerate(file_paths):
            file_path = file_path.format(threshold)
            data = read_lines(file_path)
            
            # Set x-axis ticks from 2 to 7 according to the centimorgan bins
            x_positions = np.arange(2, 8)
            data_values = [
                float(value) if float(value) >= 0 else None for value in data[i]
            ]

            if heading == "Length Discrepancy":
                # Set y-axis ticks from 1 to 16 with 1 intervals
                axes[i].set_yticks(np.arange(0, 16, 1))
                axes[i].set_ylim(0, 16) # Set y-axis limits to 1 to 16.
            else:
                # Set y-axis ticks from 0.0 to 1.2 with 0.2 intervals
                axes[i].set_yticks(np.arange(0.0, 1.2, 0.2))
                axes[i].set_ylim(0.0, 1.2)  # Set y-axis limits to 0.0 to 1.2
            
            # Assign the axes for the plots and a marker for each tool
            axes[i].plot(
                x_positions, data_values, label = legend_names[j],
                marker = markers[j], linestyle = linestyles[j],
                markerfacecolor = "none"
            )

            axes[i].set_title(heading)  # Set the plot title
            axes[i].set_xlabel("Centimorgan Bins (cM)") # Set the x-axis label
            
            # Set the y-axis label.
            if heading == "Length Discrepancy":
                axes[i].set_ylabel("Centimorgan Length (cM)")
            else:
                axes[i].set_ylabel("Precentage Decimal (% / 100)")
                
    # Remove the empty subplot that is generated
    fig.delaxes(axes[-1])
    
    # Create legend for the figure.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc = "lower right", fontsize = "14")
    
    # Add a title to the figure.
    plt.suptitle(
        "IBD Benchmark at " + str(threshold * 100) + "% Threshold "
        "and " + error_rate + "% Genotype Error Rate", fontsize = 16
    )
    
    plt.tight_layout()  # Adjust layout for better spacing
    plt.savefig(f"./{error_rate}_error/{threshold}.png") # Save the plots
    plt.close()  # Close the plots to save memory

