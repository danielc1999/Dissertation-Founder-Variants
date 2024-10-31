#!/usr/bin/python3

# This script is used to read variant information and search for IBD segments
# that go through this position. The start and end basepair positions of the
# segments will be used to plot horizontal bar plots, heatmaps and variant 
# frameworks.

import random
import pysam
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors
import sys
import os
import datetime
from collections import defaultdict, Counter

# Function to read variants information from a text file
def read_variants(filename):
    '''
    Read the lines of a text file to extract variant data.

    :param filename: The path to the input file containing the variant list.

    :return: A list with the chromosome number, basepair position, reference
    allele and alternate allele of the variants.
    '''
    variants = []
    with open(filename, "r") as file:
        for line in file:
            chr_num, bp_pos, ref_allele, alt_allele = line.strip().split()
            variants.append((chr_num, bp_pos, ref_allele, alt_allele))
    return variants

# Function to find individuals with the variant from the respective VCF
def find_variant_individuals(vcf_file, chr_num, bp_pos, ref_allele, alt_allele):
    '''
    Takes the variant data and goes to the respective VCF to identify the
    individuals that have the variant.

    :param vcf_file: The path to the respective VCF file.
    :param chr_num: The chromosome number of the variant.
    :param bp_pos: The basepair position of the variant.
    :param ref_allele: The reference allele of the variant.
    :param alt_allele: The alternate allele of the variant.

    :return: Four lists of individuals with the variant, homozygous reference
    samples, heterozygous samples and homozygous alternate sampels.
    '''
    individuals_with_variant = set()
    homo_ref_samples = set()
    hetero_samples = set()
    homo_alt_samples = set()
    bp_pos = int(bp_pos)
    chr_num = "chr" + chr_num   # Prepend "chr" to the chromosome number

    with open(vcf_file, "r") as vcf:
        header = []
        for line in vcf:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    header = line.strip().split("\t")
                continue

            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])

            if chrom != chr_num:
                continue    # Skip the chromosome if it does not match
            if pos > bp_pos:
                break   # Stop processing if base pair position goes over
            if chrom == chr_num and pos == bp_pos:
                if fields[3] == ref_allele and fields[4] == alt_allele:
                    for i, sample_info in enumerate(fields[9:]):
                        # Correctly map the sample names to genotype fields
                        sample = header[9 + i]
                        genotype = sample_info.split(":")[0]

                        # Group individuals based on genotype
                        if genotype == "0|0":  # Homozygous Reference
                            homo_ref_samples.add(sample)
                        elif genotype == "1|1":  # Homozygous Alternate
                            homo_alt_samples.add(sample)
                            individuals_with_variant.add(sample)
                        elif genotype in ["0|1", "1|0"]:  # Heterozygous
                            hetero_samples.add(sample)
                            individuals_with_variant.add(sample)

                break   # Break the loop once the variant is processed

    return individuals_with_variant, homo_ref_samples, hetero_samples, \
    homo_alt_samples

# Function to find segments with the variant from IBD file
def find_segments_with_variant(ibd_file, individuals_with_variant, bp_pos):
    '''
    Takes the variant basepair position and goes to the respective IBD file to
    extract the segments that pass through that position. One of the pair of
    individuals with the segment needs to have the variant for the IBD segment
    to be extracted.

    :param ibd_file: The path to the respective IBD file containing IBD segments
    :param individuals_with_variant:
    :param bp_pos: The basepair position of the variant.

    :return: A list containing the start and end positions of the IBD segment,
    together with the pair of individuals sharing that segment.
    '''
    segments_with_variant = []

    bp_pos = int(bp_pos)

    with open(ibd_file, "r") as ibd:
        for line in ibd:
            fields = line.strip().split()

            individual1 = fields[1]
            individual2 = fields[2]
            start_pos = int(fields[5])
            end_pos = int(fields[6])

            # Check if the position is within the segment range
            if start_pos <= bp_pos <= end_pos:
                # Check if at least 1 individual in the segment has the variant
                if (individual1 in individuals_with_variant or 
                        individual2 in individuals_with_variant):
                    segments_with_variant.append(line.strip())

    return segments_with_variant

def read_data_and_find_common_positions(segments_with_variant):
    '''
    Takes all of the IBD segments that contain the variant of interest and finds
    the common range of basepair positions which goes through all of them.

    :param segments_with_variant: A tuple containing the start and end positions
    of the IBD segment, together with the pair of individuals sharing it.

    :return: A list containing segment information, a list containing the start
    basepair position of each segment and a list containing the end basepair 
    position of each segment.
    '''
    data = []
    for line in segments_with_variant:
        entry = line.strip().split()
        entry[5] = int(entry[5])  # Start base pair positions
        entry[6] = int(entry[6])  # End base pair positions
        data.append(entry)  # Stores the entire line

    # Extract common base positions present in all individuals
    start_positions = np.array([entry[5] for entry in data])
    end_positions = np.array([entry[6] for entry in data])
    common_positions = np.arange(start_positions.max(), end_positions.min() + 1)

    if common_positions.size == 0:
        print("No common positions found.")
        sys.exit(1)

    start = common_positions[0]
    end = common_positions[-1]
    print("Common Position Range: Start Position:", start, "End Position:", end)

    return data, start, end

def plot_horizontal_bars(
    data, start, end, chr_num, bp_pos, ref_allele, alt_allele, chunk_index, 
    homo_ref_samples, hetero_samples, homo_alt_samples
):
    '''
    Plots horizontal bars to represent the IBD segments that go through the
    basepair position of the variants and highlights the common baspair
    positions between the segments.

    :param data: Information about the segment.
    :param start: The start basepair position of each segment.
    :param end: The end basepair position of each segment.
    :param chr_num: The chromosome number of the variant.
    :param bp_pos: The basepair position of the variant.
    :param ref_allele: The reference allele of the variant.
    :param alt_allele: The alternate allele of the variant.
    :param chunk_index: Index used to divide the plots into multiple figures for
    better readability. Each chunk corresponds to a separate figure.
    :param homo_ref_samples: List of homozygous reference individuals.
    :param hetero_samples: List of heterozygous individuals.
    :param homo_alt_samples: List of homozygous alternate individuals.

    :return: None
    '''
    fig, ax1 = plt.subplots(figsize=(20, 10))
    sample_labels = [entry[1] for entry in data]
    sample_labels_opposite = [entry[2] for entry in data]

    for i, entry in enumerate(data):
        start_position, end_position = entry[5], entry[6]
        individual_1 = entry[1]
        individual_2 = entry[2]

        # Check the zygosity of the individuals
        is_homo_ref_1 = individual_1 in homo_ref_samples
        is_homo_alt_1 = individual_1 in homo_alt_samples

        is_hetero_1 = individual_1 in hetero_samples
        is_hetero_2 = individual_2 in hetero_samples

        is_homo_ref_2 = individual_2 in homo_ref_samples
        is_homo_alt_2 = individual_2 in homo_alt_samples

        # Determine color based on individual statuses
        if is_homo_alt_1 and is_homo_alt_2:
            ax1.barh(
                i, end_position - start_position + 1, left = start_position, 
                height = 0.6, color = "red"
            )

        elif is_homo_ref_1 and is_homo_ref_2:
            ax1.barh(
                i, end_position - start_position + 1, left = start_position, 
                height = 0.6, color = "blue"
            )

        elif is_hetero_1 and is_hetero_2:
            ax1.barh(
                i, end_position - start_position + 1, left = start_position, 
                height = 0.6, color = "yellow"
            )

        elif (is_homo_alt_1 and is_homo_ref_2) or (is_homo_ref_1 and is_homo_alt_2):
            ax1.barh(
                i - 0.15, end_position - start_position + 1, 
                left = start_position, height = 0.3, color = "red"
            )  # Top half
            ax1.barh(
                i + 0.15, end_position - start_position + 1, 
                left = start_position, height = 0.3, color = "blue"
            )   # Bottom half

        elif (is_homo_alt_1 and is_hetero_2) or (is_hetero_1 and is_homo_alt_2):
            ax1.barh(
                i - 0.15, end_position - start_position + 1, 
                left = start_position, height = 0.3, color = "red"
            )  # Top half
            ax1.barh(
                i + 0.15, end_position - start_position + 1, 
                left = start_position, height = 0.3, color = "yellow"
            )   # Bottom half

        elif (is_hetero_1 and is_homo_ref_2) or (is_homo_ref_1 and is_hetero_2):
            ax1.barh(
                i - 0.15, end_position - start_position + 1, 
                left = start_position, height = 0.3, color = "yellow"
            )  # Top half
            ax1.barh(
                i + 0.15, end_position - start_position + 1, 
                left = start_position, height = 0.3, color = "blue"
            )   # Bottom half

    ax1.axvspan(start, end + 1, facecolor = "orange", alpha = 0.5)

    # Convert bp_pos to an integer
    bp_pos_int = int(bp_pos)

    # Set a consistent x-axis range of 500,000 basepairs each side of bp_pos
    x_min = bp_pos_int - 500000
    x_max = bp_pos_int + 500000
    ax1.set_xlim(x_min, x_max)

    # Global x-axis setup
    ax1.set_yticks(np.arange(len(data)))
    ax1.set_yticklabels(sample_labels, fontsize = 14)
    ax1.set_xlabel("Base Positions", fontsize = 16)
    ax1.set_ylabel("First Sample Identifier", fontsize = 16)

    ax2 = ax1.secondary_yaxis("right")
    ax2.set_yticks(np.arange(len(data)))
    ax2.set_yticklabels(sample_labels_opposite, fontsize = 14)
    ax2.set_ylabel("Second Sample Identifier", fontsize = 16)

    # Create a legend for the bar colors
    legend_handles = [
    plt.Line2D([0], [0], color = "orange", lw = 4, label = "Common Positions"),
    plt.Line2D([0], [0], color = "red", lw = 4, label = "Homozygous Alternate"),
    plt.Line2D([0], [0], color = "yellow", lw = 4, label = "Heterozygous"),
    plt.Line2D([0], [0], color = "blue", lw = 4, label = "Homozygous Reference")
    ]

    # Add the color legend
    ax1.legend(handles = legend_handles, loc = "upper right")

    # Set the title for the plot
    plt.title(
        f"IBD Segments for for chr{chr_num}, {bp_pos}, {ref_allele} > {alt_allele}", 
        fontsize = 18
    )

    plt.grid(True)  # Apply a grid to the plot
    plt.tight_layout()  # Apply a tight layout to the plot

    # Ensure the directory exists before saving the plot
    output_dir = os.path.join(
        ".", f"chr{chr_num}_{bp_pos}_{ref_allele}_{alt_allele}"
    )
    os.makedirs(output_dir, exist_ok = True)

    # Save the plot as a .png with the chromosome number, variant basepair
    # position, and respective reference and alternate alleles
    plt.savefig(os.path.join(output_dir, f"boxplot_{chunk_index}.png"))
    plt.close()

    print("Horizontal Bar Plot Generated.")

    return None

def create_microarray_heatmap(
    vcf_file, start, end, chr_num, bp_pos, ref_allele, alt_allele
):
    '''
    Plots microarray heatmaps for all of the individuals in the dataset
    showcasing the common basepair positions of the segments.

    :param vcf_file: The path to the VCF file.
    :param start: The start basepair position of each segment.
    :param end: The end basepair position of each segment.
    :param chr_num: The chromosome number of the variant.
    :param bp_pos: The basepair position of the variant.
    :param ref_allele: The reference allele of the variant.
    :param alt_allele: The alternate allele of the variant.

    :return: Positions, genotypes, references, alternates, basepair position
    of the variant, homozygous alternate samples, heterozygous samples, 
    homozygous reference samples, sorted samples, sorted sample indices, variant
    chromosome number, variant reference allele, variant alternate
    allele and output directory path.
    '''
    # Read all sample IDs from the VCF file header and store their
    # respective indices, genotypes, reference and alternate alleles
    #  at each position.
    sample_names = []
    sample_indices = []
    genotypes = []
    positions = []
    references = []
    alternates = []

    # Dictionary to count the occurrences of each position
    position_count = {}

    with open(vcf_file, "r") as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                sample_names = header[9:]
                sample_indices = list(range(9, len(header)))
            elif not line.startswith("#"):
                fields = line.strip().split("\t")
                pos = int(fields[1])
                if start <= pos <= end:
                    # Count occurrences of each position
                    if pos in position_count:
                        position_count[pos] += 1
                    else:
                        position_count[pos] = 1
                    genotypes.append([fields[i] for i in sample_indices])
                    positions.append(pos)
                    references.append(fields[3])
                    alternates.append(fields[4])

    # Initialize the dictionaries for grouping
    homo_ref_samples = []
    hetero_samples = []
    homo_alt_samples = []

    # Group samples based on their genotypes at bp_pos
    for i, sample_name in enumerate(sample_names):
        for j, pos in enumerate(positions):
            if pos == int(bp_pos):
                genotype = genotypes[j][i].split(":")[0]
                if genotype == "0|0":
                    homo_ref_samples.append(sample_name)
                elif genotype == "1|1":
                    homo_alt_samples.append(sample_name)
                elif genotype in ("0|1", "1|0"):
                    hetero_samples.append(sample_name)
                break

    # Identify positions that occur more than once (for indels)
    duplicate_positions = {
        pos for pos, count in position_count.items() if count > 1
    }

    # Filter out variants with duplicate positions
    filtered_positions = []
    filtered_references = []
    filtered_alternates = []
    filtered_genotypes = []

    # Track already processed positions to avoid repeats
    processed_positions = set()

    for pos, ref, alt, geno in zip(positions, references, alternates, genotypes):
        if pos not in duplicate_positions:
            # Kepp the the position if it is not duplicated
            filtered_positions.append(pos)
            filtered_references.append(ref)
            filtered_alternates.append(alt)
            filtered_genotypes.append(geno)
            processed_positions.add(pos)

        else:
            # Dictionary to store inserion counts
            insertion_counts = defaultdict(int)
            # Dictionary to store deletion counts
            deletion_counts = defaultdict(int)

            # Find a sample that is homozygous for the alternate allele (1|1)
            for sample_name, genotype in zip(sample_names, geno):
                if sample_name in homo_alt_samples and pos not in processed_positions:
                    # Use the sample name to locate the corresponding BAM file
                    bam_file_path = (
                            f"/data/targetid/dna/results/dnaseq_results/BAM/"
                            f"{sample_name}.post.sortNMMD.bam"
                    )
                    samfile = pysam.AlignmentFile(bam_file_path, "rb")

                    for pileupcolumn in samfile.pileup(
                                            f"chr{chr_num}", pos, pos + 1
                    ):
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                read = pileupread.alignment.query_sequence
                                indel_pos = pileupread.query_position

                                # Insertion: check if alt is longer than ref
                                if len(alt) > len(ref):
                                    indel_chars = pileupread.indel

                                    # Ensure we are within bounds of the read
                                    if (
                                        indel_chars > 0 and
                                        indel_pos + 1 + indel_chars <= len(read)
                                    ):
                                        inserted_alleles = read[
                                            indel_pos + 1:
                                            indel_pos + 1 + indel_chars
                                        ]

                                        # Increment count for this inserted allele
                                        insertion_counts[inserted_alleles] += 1

                                # Deletion: check if ref is longer than alt
                                elif len(ref) > len(alt):
                                    indel_chars = pileupread.indel

                                    # Ensure the deletion is within bounds
                                    if indel_chars < 0 and indel_pos + 1 <= len(read):
                                        deleted_alleles = read[
                                            indel_pos + 1:
                                            indel_pos + 1 + abs(indel_chars)
                                        ]

                                        # Increment count for this deleted allele
                                        deletion_counts[deleted_alleles] += 1

                    samfile.close()

            # Find the most frequent insertion
            if insertion_counts:
                most_frequent_insertion = max(
                    insertion_counts, key = insertion_counts.get
                )

                # Check if the most frequent insertion matches the ALT allele
                if most_frequent_insertion == alt[1:]:
                    # Keep the variant
                    filtered_positions.append(pos)
                    filtered_references.append(ref)
                    filtered_alternates.append(alt)
                    filtered_genotypes.append(geno)
                    # Mark this position as processed
                    processed_positions.add(pos)

            # Find the most frequent deletion length
            if deletion_counts:
                most_frequent_deletion = max(
                    deletion_counts, key = deletion_counts.get
                )

                # Check if the most frequent deletion length matches the
                # expected deletion length
                if most_frequent_deletion == ref[len(alt):]:
                    # Keep the variant
                    filtered_positions.append(pos)
                    filtered_references.append(ref)
                    filtered_alternates.append(alt)
                    filtered_genotypes.append(geno)
                    # Mark this position as processed
                    processed_positions.add(pos)

    # If no positions remain after filtering, return empty results
    if not filtered_positions:
        return positions, references, alternates, genotypes, sample_names, ""

    positions = filtered_positions
    references = filtered_references
    alternates = filtered_alternates
    genotypes = filtered_genotypes

    color_mapping = {"0|0": 0, "1|1": 1, "0|1": 2, "1|0": 3}
    # Colours correspond to the colour mapping
    custom_colors = ["blue", "red", "yellow", "orange"]
    other_color = "gray" # Other for any genotypes with missing data

    numerical_data = np.array(
   [[color_mapping.get(gt.split(":")[0], 4) for gt in row] for row in genotypes]
    )

    # Create the colourmap
    cmap = matplotlib.colors.ListedColormap(custom_colors + [other_color])

    # Create the heatmap
    fig, ax = plt.subplots(figsize = (20, 20))

    # Sort the samples based on the group
    sorted_samples = homo_alt_samples + hetero_samples + homo_ref_samples
    sorted_indices = [sample_names.index(sample) for sample in sorted_samples]

    sorted_numerical_data = numerical_data[:, sorted_indices].T

    cax = ax.matshow(
        sorted_numerical_data, cmap = cmap, aspect = "auto", vmin = 0, vmax = 4
    )

    # Set x-ticks for sample IDs
    num_samples = len(sorted_samples)
    ax.set_yticks(np.arange(num_samples))
    ax.set_yticklabels(sorted_samples, fontsize = 10)
    ax.set_ylabel("Sample IDs")

    # Set x-ticks for base pair positions
    if len(positions) <= 65:
        ax.set_xticks(np.arange(len(positions)))
        ax.set_xticklabels(positions, rotation = 90, fontsize = 20)
    else:
        # Calculate the position interval for when many positions are present
        interval = len(positions) // 65  
        tick_positions = positions[::interval]  # Sample positions at intervals
        # Set x-ticks at intervals
        ax.set_xticks(np.arange(0, len(positions), interval))  
        ax.set_xticklabels(tick_positions, rotation = 90, fontsize = 20)

    ax.set_xlabel("Base Pair Positions")

    # Set heatmap title
    ax.set_title(
        f"Heatmap for chr{chr_num}, {bp_pos}, {ref_allele} > {alt_allele}",
        fontsize = 22
    )

    # Create a colour legend for the heatmap
    legend_labels = ["0|0", "1|1", "0|1", "1|0"]
    legend_handles = [
        plt.Line2D([0], [0], color = custom_colors[0], lw = 4),
        plt.Line2D([0], [0], color = custom_colors[1], lw = 4),
        plt.Line2D([0], [0], color = custom_colors[2], lw = 4),
        plt.Line2D([0], [0], color = custom_colors[3], lw = 4),
    ]

    ax.legend(
        legend_handles, legend_labels, title = "Genotype",
        loc = "upper left", bbox_to_anchor = (1, 1),
        fontsize = 12, title_fontsize = 14,
    )
    plt.tight_layout()

    # Add vertical lines for group boundaries
    boundaries = [
        len(homo_alt_samples), len(homo_alt_samples) + len(hetero_samples)
    ]

    for boundary in boundaries:
        ax.axhline(
            boundary - 0.5, linestyle = "--", color = "black", linewidth = 1
        )

    # Ensure the directory exists before saving the heatmap
    output_dir = os.path.join(
        ".", f"chr{chr_num}_{bp_pos}_{ref_allele}_{alt_allele}"
    )
    os.makedirs(output_dir, exist_ok = True)

    # Save the heatmap as a .png
    plt.savefig(
      os.path.join(
        output_dir, 
            f"microarray_heatmap_chr{chr_num}_{bp_pos}_{ref_allele}_{alt_allele}.png", 
        )
    )
    plt.close()

    return positions, genotypes, references, alternates, bp_pos, \
    homo_alt_samples, hetero_samples, homo_ref_samples, sorted_samples, \
    sorted_indices, chr_num, ref_allele, alt_allele, output_dir

def create_mini_microarray(
    positions, genotypes, bp_pos, homo_alt_samples, hetero_samples,  
    homo_ref_samples, sorted_samples, sorted_indices, chr_num, ref_allele, 
    alt_allele
):
    '''
    Plots a smaller heatmap for all individuals with the variant and a random 50
    homozygous refrence individuals from the dataset. This showcases the common
    basepair positions of the IBD segments.

    :param positions: A list of basepair positions from the VCF dataset,
    representing the IBD segment.
    :param genotypes: A list of sample genotypes corresponding to each position.
    :param bp_pos: THe basepair positon of the variant.
    :param homo_alt_samples: A list of homozygous alternate sampels for the
    variant.
    :param hetero_samples: A list of heterozygous samples for the variant.
    :param homo_ref_samples: A list of homozygous reference samples for the
    variant.
    :param sorted_samples: A list of samples sorted according to their zygosity.
    :param sorted_indices: A list of sample indices sorted according to the
    sorted_samples variable.
    :param chr_num: The chromosome number of the variant.
    :param ref_allele: The reference allele of the variant.
    :param alt_allele: The alternate allele of the variant.

    :return: None.
    '''

    color_mapping = {"0|0": 0, "1|1": 1, "0|1": 2, "1|0": 3}
    # Colours correspond to the colour mapping
    custom_colors = ["blue", "red", "yellow", "orange"]
    other_color = "gray" # Other for any genotypes with missing data

    numerical_data = np.array(
   [[color_mapping.get(gt.split(":")[0], 4) for gt in row] for row in genotypes]
    )
    sorted_numerical_data = numerical_data[:, sorted_indices].T

    # Filter out the basepair positions where all individuals are 0|0
    non_zero_positions = np.any(sorted_numerical_data != 0, axis = 0)
    filtered_positions = np.array(positions)[non_zero_positions]
    filtered_numerical_data = sorted_numerical_data[:, non_zero_positions]

    # Find the index of bp_pos in filtered_positions
    bp_index = np.where(filtered_positions == int(bp_pos))[0][0]

    # Define limits for 40 positions before and after
    start = max(0, bp_index - 40)
    end = min(len(filtered_positions), bp_index + 40 + 1)

    # Slice the filtered positions and genotypes for the desired range
    filtered_positions = filtered_positions[start:end]
    filtered_numerical_data = filtered_numerical_data[:, start:end]

    # Limit the number of homozygous reference samples to 50 and chosen randomly
    if len(homo_ref_samples) > 50:
        homo_ref_samples = random.sample(homo_ref_samples, 50)

    # Combine the chosen homo_ref_samples with all hetero_samples and
    # homo_alt_samples
    combined_samples = homo_alt_samples + hetero_samples + homo_ref_samples

    # Filter the numerical data for the combined samples
    filtered_numerical_data = filtered_numerical_data[
        [s in combined_samples for s in sorted_samples], :
    ]

    # Create the heatmap
    cmap = matplotlib.colors.ListedColormap(custom_colors + [other_color])
    fig, ax = plt.subplots(figsize = (20, 20))

    cax = ax.matshow(
       filtered_numerical_data, cmap = cmap, aspect = "auto", vmin = 0, vmax = 4
    )

    # Update y-ticks to reflect the combined samples
    ax.set_yticks(np.arange(len(combined_samples)))
    ax.set_yticklabels(combined_samples, fontsize = 10)
    ax.set_ylabel("Sample IDs")

    # Use the filtered positions and numerical data in plotting
    ax.set_xticks(np.arange(len(filtered_positions)))
    ax.set_xticklabels(filtered_positions, rotation = 90, fontsize = 20)
    ax.set_xlabel("Base Pair Positions")

    # Set heatmap title
    ax.set_title(
        f"Mini Heatmap for chr{chr_num}, {bp_pos}, {ref_allele} > {alt_allele}", 
        fontsize = 22
    )

    # Create a colour legend for the heatmap
    legend_labels = ["0|0", "1|1", "0|1", "1|0"]
    legend_handles = [
        plt.Line2D([0], [0], color = custom_colors[0], lw = 4),
        plt.Line2D([0], [0], color = custom_colors[1], lw = 4),
        plt.Line2D([0], [0], color = custom_colors[2], lw = 4),
        plt.Line2D([0], [0], color = custom_colors[3], lw = 4),
    ]

    ax.legend(
        legend_handles, legend_labels, title = "Genotype",
        loc = "upper left", bbox_to_anchor = (1, 1)
    )
    plt.tight_layout()

    # Add vertical lines for group boundaries
    boundaries = [
        len(homo_alt_samples), len(homo_alt_samples) + len(hetero_samples)
    ]

    for boundary in boundaries:
        ax.axhline(
            boundary - 0.5, linestyle = "--", color = "black", linewidth = 1
        )

    # Ensure the directory exists before saving the heatmap
    output_dir = os.path.join(
        ".", f"chr{chr_num}_{bp_pos}_{ref_allele}_{alt_allele}"
    )
    os.makedirs(output_dir, exist_ok = True)

    # Save the heatmap as a .png
    plt.savefig(
        os.path.join(
            output_dir, 
          f"mini_microarray_chr{chr_num}_{bp_pos}_{ref_allele}_{alt_allele}.png"
        )
    )
    plt.close()

    return None

def plot_reference_vs_variant_allele(
    variant, positions, references, alternates, genotypes,
    homo_alt_samples, hetero_samples, homo_ref_samples, output_dir
):
    '''
    Plots the variant framework of the individuals that have the variant. It
    also filters out incorrect indels present within the dataset.

    :param variant: The chromosome number, basepair position, reference allele
    and alternate allele of the variant.
    :param positions: A list of positions for the IBD segment.
    :param references: A list of reference alleles for every position in the IBD
    segment.
    :param alternates: A list of alternate alleles for every position in the IBD
    segment.
    :param genotypes: A list of genotypes for every position in the IBD segment.
    :param homo_alt_samples: A list of homozygous alternate sampels for the
    variant.
    :param hetero_samples: A list of heterozygous samples for the variant.
    :param homo_ref_samples: A lit of homozygous reference samples for the
    variant.
    :param output_dir: The output directory file path.

    :return: None.
    '''
    chr_num, bp_pos, ref_allele, alt_allele = variant
    sample_names = homo_alt_samples + hetero_samples + homo_ref_samples

    # Count occurrences of each position
    position_count = Counter(positions)

    # Identify positions that occur more than once (for indels)
    duplicate_positions = {
        pos for pos, count in position_count.items() if count > 1
    }

    # Filter out variants with duplicate positions
    filtered_positions = []
    filtered_references = []
    filtered_alternates = []
    filtered_genotypes = []

    # Track already processed positions to avoid repeats
    processed_positions = set()

    for pos, ref, alt, geno in zip(positions, references, alternates, genotypes):
        if pos not in duplicate_positions:
            # If the position is not duplicated, keep it
            filtered_positions.append(pos)
            filtered_references.append(ref)
            filtered_alternates.append(alt)
            filtered_genotypes.append(geno)
            processed_positions.add(pos)

        else:
            # Dictionary to store inserion counts
            insertion_counts = defaultdict(int) 
            # Dictionary to store deletion counts
            deletion_counts = defaultdict(int)

            # Find a sample that is homozygous for the alternate allele (1|1)
            for sample_name, genotype in zip(sample_names, geno):
                if genotype == "1|1" and pos not in processed_positions:
                    # Use the sample name to locate the corresponding BAM file
                    bam_file_path = (
                        f"/data/targetid/dna/results/dnaseq_results/BAM/"
                        f"{sample_name}.post.sortNMMD.bam"
                    )
                    samfile = pysam.AlignmentFile(bam_file_path, "rb")

                    # Dictionary to store insertion counts
                    insertion_counts = defaultdict(int)
                    # Dictionary to store deletion counts
                    deletion_counts = defaultdict(int)

                    for pileupcolumn in samfile.pileup(
                            f"chr{chr_num}", pos, pos + 1
                    ):
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                read = pileupread.alignment.query_sequence
                                indel_pos = pileupread.query_position
                                cigar = pileupread.alignment.cigarstring

                                # Insertion: check if alt is longer than ref
                                if len(alt) > len(ref):
                                    indel_chars = pileupread.indel

                                    # Ensure we are within bounds of the read
                                    if (
                                        indel_chars > 0 and 
                                        indel_pos + 1 + indel_chars <= len(read)
                                    ):
                                        inserted_alleles = read[
                                            indel_pos + 1:
                                            indel_pos + 1 + indel_chars
                                        ]

                                        # Increment count for inserted allele
                                        insertion_counts[inserted_alleles] += 1

                                # Deletion: check if ref is longer than alt
                                elif len(ref) > len(alt):
                                    indel_chars = pileupread.indel

                                    # Ensure the deletion is within bounds
                                    if (
                                        indel_chars < 0 and 
                                        indel_pos + 1 <= len(read)
                                    ):
                                        deleted_alleles = read[
                                            indel_pos + 1:
                                            indel_pos + 1 + abs(indel_chars)
                                        ]

                                        # Increment count for deleted allele
                                        deletion_counts[deleted_alleles] += 1

                    samfile.close()

            # Find the most frequent insertion
            if insertion_counts:
                most_frequent_insertion = max(
                    insertion_counts, key = insertion_counts.get
                )

                # Check if the most frequent insertion matches the ALT allele
                if most_frequent_insertion == alt[1:]:
                    # Keep the variant
                    filtered_positions.append(pos)
                    filtered_references.append(ref)
                    filtered_alternates.append(alt)
                    filtered_genotypes.append(geno)
                    # Mark this position as processed
                    processed_positions.add(pos)

            # Find the most frequent deletion length
            if deletion_counts:
                most_frequent_deletion = max(
                    deletion_counts, key = deletion_counts.get
                )

                # Check if the most frequent deletion length matches the
                # expected deletion length
                if most_frequent_deletion == ref[len(alt):]:
                    # Keep the variant
                    filtered_positions.append(pos)
                    filtered_references.append(ref)
                    filtered_alternates.append(alt)
                    filtered_genotypes.append(geno)
                    # Mark this position as processed
                    processed_positions.add(pos)

    # If no positions remain after filtering, return empty results
    if not filtered_positions:
        return positions, references, alternates, genotypes, sample_names, ""

    positions = filtered_positions
    # Truncate references and alternates when assigning them
    references = [
        ref[:5] + "..." if len(ref) > 5 else ref
        for ref in filtered_references
    ]

    alternates = [
        alt[:5] + "..." if len(alt) > 5 else alt
        for alt in filtered_alternates
    ]
    genotypes = filtered_genotypes

    # Find the index of the variant position
    variant_index = positions.index(int(bp_pos))

    # Filter genotypes for individuals with "1" allele
    individuals_with_variant = {
     sample_names[i] for i, g in enumerate(genotypes[variant_index]) if "1" in g
    }

    variant_alleles_by_position = defaultdict(lambda: defaultdict(int))

    for i, sample in enumerate(sample_names):
        if sample in individuals_with_variant:
            for pos_idx in range(len(positions)):
                alleles = genotypes[pos_idx][i].split("|")
                if "1" in alleles:
                    variant_alleles_by_position[positions[pos_idx]][alternates[pos_idx]] += 1

    # Identify alleles present in the majority of individuals (90% of them)
    # Value can be changed
    required_count = int(len(individuals_with_variant) * 0.9)
    variant_framework = defaultdict(list)

    for pos, allele_counts in variant_alleles_by_position.items():
        for allele, count in allele_counts.items():
            if count >= required_count:
                variant_framework[pos] = allele
                break

    # Filter out positions where the most common allele matches the reference 
    # allele
    positions_common = []
    references_common = []
    common_alleles = []
    for pos in positions:
        if pos in variant_framework:
            reference_allele = references[positions.index(pos)]
            common_allele = variant_framework[pos]
            if reference_allele != common_allele:
                positions_common.append(pos)
                references_common.append(reference_allele)
                common_alleles.append(common_allele)

    # Determine the number of chunks
    window_size = 100
    num_chunks = (
        len(positions_common) // window_size + 
        (1 if len(positions_common) % window_size != 0 else 0)
    )

    # Plotting the reference sequence and common allele sequence in chunks
    for chunk_idx in range(num_chunks):
        start_idx = chunk_idx * window_size
        end_idx = min((chunk_idx + 1) * window_size, len(positions_common))
        chunk_positions = positions_common[start_idx:end_idx]
        chunk_references = references_common[start_idx:end_idx]
        chunk_common_alleles = common_alleles[start_idx:end_idx]

        # Adjust positions to avoid overlapping
        adjusted_positions = list(range(start_idx, end_idx))

        fig, ax = plt.subplots(figsize = (20, 5), constrained_layout = True)

        # Plot reference sequence
        ax.scatter(
            adjusted_positions, [0.75]*len(chunk_positions), s = 10, 
            label = "Reference Sequence", color = "blue"
        )

        # Plot common allele sequence without color differentiation
        for adj_pos, alt_allele in zip(adjusted_positions, chunk_common_alleles):
            ax.scatter(adj_pos, 0.25, s = 10, color = "red")
            ax.annotate(
                alt_allele, xy = (adj_pos, 0.25), xytext = (-5, 5), 
                textcoords = "offset points", ha = "center", 
                rotation = 90, fontsize = 11
            )

        # Annotate reference alleles
        for adj_pos, reference_allele in zip(adjusted_positions, chunk_references):
            ax.annotate(
                reference_allele, xy = (adj_pos, 0.75), xytext = (-5, 5), 
                textcoords = "offset points", ha = "center", 
                rotation = 90, fontsize = 11
            )

        # Add a straight line passing through both scatter plot lines
        ax.plot(
            adjusted_positions, [0.75]*len(adjusted_positions), linestyle = "--"
        )
        ax.plot(
            adjusted_positions, [0.25]*len(adjusted_positions), linestyle = "--"
        )

        # Customize plot
        ax.set_xticks(adjusted_positions)
        ax.set_xticklabels(
            [str(pos) for pos in chunk_positions], rotation = 90, fontsize = 14
        )
        ax.set_yticks([0.25, 0.75])
        ax.set_yticklabels(["Variant Allele", "Reference"], fontsize = 14)
        ax.set_xlabel("Base Pair Positions", fontsize = 16)
        ax.set_title(
            f"Framework for variant chr{chr_num}, {bp_pos}, "
            f"{variant[2]} > {variant[3]} ({chunk_idx + 1})", 
            fontsize = 18, pad = 20
        )

        # Remove the upper line of the box representing the plot
        ax.spines["top"].set_visible(False)

        # Save plot to file
        plt.savefig(os.path.join(
            output_dir, f"reference_vs_common_allele_chunk_{chunk_idx + 1}.png")
        )
        plt.close()

    print("Common variant sequence constructed.")

    # Save results to files
    current_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Write common variant sequence to a separate txt file
    with open(os.path.join(output_dir, "variant_framework.txt"), "w") as f:
        f.write(f"File generated on: {current_time}\n")
        f.write(
            f"Variant: chr{chr_num}, {bp_pos}, {variant[2]} > {variant[3]}\n"
        )

        f.write("\nVariant framework for individuals with the variant:\n")
        f.write(
            f"{'Bp Position':<15}{'Reference Allele':<20}{'Variant Allele'}\n"
        )
        f.write("-" * 50 + "\n")

        for pos, ref_allele, common_allele in zip(
                    positions_common, filtered_references, common_alleles
        ):
            # Writes the alleles in full rather than using the short version
            full_ref_allele = filtered_references[positions.index(pos)]
            full_common_allele = variant_framework[pos]
            f.write(f"{pos:<15}{full_ref_allele:<20}{full_common_allele}\n")

    print(
    f"Variant framework sequences have been written to files in {output_dir}.\n"
    )

    return None

# Main function
def main():
    """
    The main function of the script.

    This function combines the processing of genomic data, including loading
    VCF files, finding IBD segments and generating figures. The latter includes
    horizontal bar plots, heatmaps and variant frameworks. The outputs are saved 
    to a specified directory.

    :return: None
    """

    variants = read_variants("variants.txt")

    for variant in variants:
        chr_num, bp_pos, ref_allele, alt_allele = variant

        print(
            f"Processing variant: chr{chr_num}, {bp_pos}, "
            f"{ref_allele}, {alt_allele}"
        )

        # Adjusts the VCF file path based on the chromosome number
        vcf_file = f"./chromosome{chr_num}/chr{chr_num}.phased.recode.vcf"

        individuals_with_variant, homo_ref_samples, hetero_samples, \
        homo_alt_samples = find_variant_individuals(
            vcf_file, chr_num, bp_pos, ref_allele, alt_allele
        )

        # Print the count of individuals with the variant
        count = len(individuals_with_variant)
        print(f"Number of individuals with the variant: {count}")

        if not individuals_with_variant:
            print(
                f"No individuals with the variant chr{chr_num}, {bp_pos}, "
                f"{ref_allele}, {alt_allele} found.\n"
            )
            continue    # Move to the next variant if no individuals have it

        # Adjusts the IBD segment file path based on the chromosome number
        ibd_file = f"./chromosome{chr_num}/chr{chr_num}_0.5cM.max"
        segments_with_variant = find_segments_with_variant(
            ibd_file, individuals_with_variant, bp_pos
        )

        if not segments_with_variant:
            print(
                f"No segments with the variant chr{chr_num}, {bp_pos}, "
                f"{ref_allele}, {alt_allele} found.\n"
            )
            continue  # Move to the next variant if no segments are found

        data, start, end = read_data_and_find_common_positions(
            segments_with_variant
        )

        # Call the function for the creation of the microarray heatmap
        create_microarray_heatmap(
            vcf_file, start, end, chr_num, bp_pos, ref_allele, alt_allele
        )

        mini_array = False

        if mini_array:
            positions, genotypes, references, alternates, bp_pos, \
            homo_alt_samples, hetero_samples, homo_ref_samples, sorted_samples, \
            sorted_indices, chr_num, ref_allele, alt_allele = \
            create_microarray_heatmap(
                vcf_file, start, end, chr_num, bp_pos, ref_allele, alt_allele
            )

            create_mini_microarray(
                positions, genotypes, bp_pos, homo_alt_samples, hetero_samples, 
                homo_ref_samples, sorted_samples, sorted_indices, chr_num, \
                ref_allele, alt_allele
            )

        # Define the threshold and split data into chunks if necessary
        threshold = 53  # Value can be adjusted

        # Divide the data into chunks to process into the heatmap function
        chunks = [data[i:i + threshold] for i in range(0, len(data), threshold)]
        for i, chunk in enumerate(chunks):
            plot_horizontal_bars(
                chunk, start, end, chr_num, bp_pos, ref_allele, alt_allele, i, \
                homo_ref_samples, hetero_samples, homo_alt_samples
            )

        positions, genotypes, references, alternates, bp_pos, homo_alt_samples, \
        hetero_samples, homo_ref_samples, sorted_samples, sorted_indices, \
        chr_num, ref_allele, alt_allele, output_dir = \
        create_microarray_heatmap(
            vcf_file, start, end, chr_num, bp_pos, ref_allele, alt_allele
        )

        # Call the function for the creation of the reference vs common variant
        # allele framework sequence
        plot_reference_vs_variant_allele(
            variant, positions, references, alternates, genotypes, 
            homo_alt_samples, hetero_samples, homo_ref_samples,output_dir
        )

    return None

if __name__ == "__main__":
    main()
