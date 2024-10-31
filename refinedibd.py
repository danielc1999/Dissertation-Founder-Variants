#!/usr/bin/python3

# Used to modify RefinedIBD .ibd files to the hap-IBD parser format to be used
# in the IBD benchmark tool.

# Open the RefinedIBD .ibd file.
with open("refinedibd_eur.0.01.ibd", "r") as input_file:
    # Open a new output file.
    with open("modified_refinedibd_eur.0.01.ibd", "w") as output_file:
        for line in input_file:
            # Split the tab-separated line into individual values.
            values = line.split("\t")

            # Assign values to variables based on the desired indices.
            id1 = 0         # First sample identifier
            id2 = 2         # Second sample identifier
            hap1 = 1        # First sample haplotype index
            hap2 = 3        # Second sample haplotype index
            chromosome = 4  # Chromosome number
            start = 5       # Starting genomic position
            end = 6         # Ending genomic position
            length = 8      # Length of IBD segment in centiMorgan

            # Write the formatted line to the output file according to the
            # Hap-IBD parser format of the IBD benchmark tool.
            output_file.write(
                f"{values[id1]}\t{values[hap1]}\t{values[id2]}\t{values[hap2]}"
                f"\t{values[chromosome]}\t{values[start]}\t{values[end]}\t"
                f"{values[length]}"
            )

