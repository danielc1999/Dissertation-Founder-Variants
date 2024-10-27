#!/usr/bin/python3

# Used to modify RefinedIBD .ibd files to the Hap-IBD parser format to be used
# in the IBD benchmark tool.

# Open the RefinedIBD .ibd file.
with open("refinedibd_eur.0.01.ibd", "r") as input_file:
    # Open a new output file.
    with open("modified_refinedibd_eur.0.01.ibd", "w") as output_file:
        for line in input_file:
            # Split the tab-separated line into individual values.
            values = line.split("\t")

            # Assign values to variables based on the desired indices.
            ID1_idx = 0
            ID2_idx = 2
            Hap1_idx = 1
            Hap2_idx = 3
            chromosome = 4
            PhyStart_idx = 5
            PhyEnd_idx = 6
            Len_idx = 8

            # Write the formatted line to the output file according to the
            # Hap-IBD parser format of the IBD benchmark tool.
            output_file.write(
        f"{values[ID1_idx]}\t{values[Hap1_idx]}\t{values[ID2_idx]}\t"
        f"{values[Hap2_idx]}\t{values[chromosome]}\t{values[PhyStart_idx]}\t"
        f"{values[PhyEnd_idx]}\t{values[Len_idx]}"
        )

