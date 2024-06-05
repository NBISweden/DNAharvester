#!/usr/bin/env python3

from sys import argv
from itertools import islice
import gzip

def parse_haplo(haplo_file):
    # Define a list of sample names
    samples = ["SF-1_42", "AI-1_119", "BL-11_9", "KU-23_3",
               "MO34_ZA", "HD-3", "RY-3", "EL012", "Sheep", "domesticgoat"]

    # Create an empty dictionary to store the sequences for each sample
    sample_dict = {}
    for i in samples:
        sample_dict[i] = ""

    # Open the haplo_file using gzip and iterate over each line
    with gzip.open(haplo_file, 'rt') as f1:
        next(f1)  # Skip the header line
        for line in f1:
            splitted = line.strip().split("\t")
            # Splits each line by tabs and extracts chrom, pos, and alleles.
            # Alleles is formed by concatenating elements in the line starting from the fourth element onwards (splitted[3:]).
            chrom, pos, alleles = splitted[0], splitted[1], "".join(splitted[3:])
            # Filters out lines where alleles contains 'N'.
            if "N" not in alleles:
                # Checks if alleles contains exactly two unique values.
                if len(set(alleles)) == 2:
                    # Appends the respective allele (character) to the corresponding sample's sequence in sample_dict.
                    for i in range(len(alleles)):
                        focal_sample = samples[i]
                        sample_dict[focal_sample] += alleles[i]

    # Print the sequences for each sample
    for key, value in sample_dict.items():
        print(">" + key)
        print(value)

if __name__ == "__main__":
    # Get the filename from command line arguments
    filename = argv[1]
    # Call the parse_haplo function with the filename
    parse_haplo(filename)
