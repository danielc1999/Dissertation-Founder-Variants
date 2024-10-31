#!/bin/bash

# This script is used to automatically generate IBD segments using RaPID for each chromosome. The chromosome number needs to be specified by the user while executing the script.

# Check if chromosome number is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <chromosome_number>"
    exit 1
fi

# Set the chromosome number from the first argument
chromosome=$1

# Filter by chromosome using VCFtools.
vcftools --vcf unphased_merged_relabelled_VCF.recode.vcf --chr chr${chromosome} --recode --out chr${chromosome}.unphased
echo "Filtering by chromosome and individuals complete."

# Perform Beagle phasing using plink genetic maps.
java -jar beagle.01Mar24.d36.jar gt=chr${chromosome}.unphased.recode.vcf map=plink_maps_GRCh38/plink.chr${chromosome}.GRCh38.map out=chr${chromosome}.phased.recode nthreads=20 window=10
echo "Phasing complete."

# Change to RaPID's directory.
cd RaPID/

#Perform genetic map filtering and interpolation of loci for RaPID using RaPID's genetic maps.
python filter_mapping_file.py genetic_maps/genetic_map_GRCh38_chr${chromosome}.txt filtered_chr${chromosome}_map.txt
python interpolate_loci.py filtered_chr${chromosome}_map.txt ../chr${chromosome}.phased.recode.vcf.gz chr${chromosome}_map
echo "Interpolation of loci complete."

# Perform RaPID IBD detection using its optimal parameters. The centiMorgan threshold parameter (-d) can be adjusted.
./RaPID_v.1.7 -w 5 -r 10 -s 2 -d 0.5 -i ../chr${chromosome}.phased.recode.vcf.gz -g chr${chromosome}_map -o ../chr${chromosome}.phased
echo "RaPID IBD detection complete."

# Sort the generated files into respective directories on a chromosomal basis.
cd
mkdir chromosome${chromosome}

cd chr${chromosome}.phased/
gunzip results.max.gz
mv results.max chr${chromosome}.max
mv chr${chromosome}.max ../chromosome${chromosome}/

cd
mv chr${chromosome}.unphased.recode.vcf chromosome${chromosome}
mv chr${chromosome}.unphased.log chromosome${chromosome}
mv chr${chromosome}.phased.recode.vcf.gz chromosome${chromosome}
mv chr${chromosome}.phased.recode.log chromosome${chromosome}
mv ./RaPID/filtered_chr${chromosome}_map.txt chromosome${chromosome}
mv ./RaPID/chr${chromosome}_map chromosome${chromosome}

rm -rf chr${chromosome}.phased
