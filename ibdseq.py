#!/usr/bin/python3

# This script calculates the centiMorgan distance of IBD segments for IBDSeq's 
# .ibd files by using a genetic map. It also changes the format of the output
# to match Hap-IBD's parser for the IBD Benchmark tool.

import bisect

def find_closest_genetic_map_position(
    genomic_position, genomic_positions, map_lines
):
    '''
    Find the closest genetic map position for a given genomic position.
    
    This function takes a basepair position and finds the closest position from
    the genetic map file.
    
    :param genomic_position: The genomic basepair position to search for.
    :param map_lines: The respective closest genomic basepair position.
    
    :return: Returns the closest basepair positions.
    '''
    index = bisect.bisect_left(genomic_positions, genomic_position)

    if index == 0:
        closest_position = map_lines[0]
    elif index == len(map_lines):
        closest_position = map_lines[-1]
    else:
        left_position = map_lines[index - 1]
        right_position = map_lines[index]
        closest_position = left_position if abs(
                int(left_position.split("\t")[1]) - genomic_position
                ) < abs(
                    int(right_position.split("\t")[1]) - genomic_position
                        ) else right_position

    closest_position = closest_position.split("\t")
    
    return int(closest_position[1]), float(closest_position[3])

# Load the genetic map file
with open("genetic_map_GRCh37_chr20.txt", "r") as map_file:
    next(map_file)  # Skip the header line
    map_lines = map_file.readlines()    # Read the rest of the file

# Create genomic positions list outside the loop
genomic_positions = [int(line.split("\t")[1]) for line in map_lines]

# Load the IBD file
with open("ibdseq_eur.0.1.ibd", "r") as ibd_file:
    # Write the results to the output file immediately after calculation
    with open("modified_ibdseq_eur.0.1.ibd", "w") as output_file:
        for line in ibd_file:
            # Split the line into individual values
            values = line.split("\t")

            # Assign values to variables based on the desired indices
            id1, hap1, id2, hap2, chromosome, start, end = 0, 1, 2, 3, 4, 5, 6

            # Extract the relevant values
            id1, hap1, id2, hap2, chromosome, start, end = (
                values[id1], values[hap1], values[id2], values[hap2], 
                values[chromosome], values[start], values[end]
            )

            # Find the closest genetic map positions for the start and
            # end coordinates
            genomic_start = int(start)
            genomic_end = int(end)

            closest_start_position, closest_start_map = (
                find_closest_genetic_map_position(
                    genomic_start, genomic_positions, map_lines
                    )
            )

            closest_end_position, closest_end_map = (
                find_closest_genetic_map_position(
                    genomic_end, genomic_positions, map_lines
                    )
            )

            # Calculate centimorgan distance for the IBD segments
            centimorgan_distance = round(closest_end_map - closest_start_map, 3)

            # Write the formatted line to the output file immediately
            output_file.write(
                f"{id1}\t{hap1}\t{id2}\t{hap2}\t{chromosome}\t{start}"
                f"\t{end}\t{centimorgan_distance}\n"
            )

