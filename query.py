#!/usr/bin/python3

# Used to modify RaPID-Query .ibd files to the hap-IBD parser format to be used
# in the IBD benchmark tool.

# Open the RefinedIBD .ibd file.
with open("query.0.1.ibd", "r") as input_file:
    # Open a new output file.
    with open("modified_query.0.1.ibd", "w") as output_file:
        for index, line in enumerate(input_file):
            if index == 0:
                continue  # Skip the first line
            # Split the tab-separated line into individual values.
            values = line.split(",")

            # Assign values to variables based on the desired indices.
            id1 = 0             # First sample identifier
            id2 = 2             # Second sample identifier
            hap1 = 1            # First sample haplotype index
            hap2 = 3            # Second sample haplotype index
            physical_start = 4  # Starting genomic position
            physical_end = 5    # Ending genomic position
            length = 6          # Length in centiMorgans
            start_index = 7     # Starting index position
            end_index = 8       # Ending index position

            # Write the formatted line to the output file according to the
            # Hap-IBD parser format of the IBD benchmark tool.
            output_file.write(
                f"20\t{values[id1]}\t{values[id2]}\t{values[hap1]}\t"
                f"{values[hap2]}\t{values[physical_start]}\t"
                f"{values[physical_end]}\t{values[length]}\t"
                f"{values[start_index]}\t{values[end_index]}"
            )

