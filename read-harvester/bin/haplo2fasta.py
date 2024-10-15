#!/usr/bin/env python3

from sys import argv
from itertools import islice
import gzip

def parse_haplo(haplo_file, samplename):
    # Create an empty dictionary to store the sequence with a samplename
    sample_dict = {}
    sample_dict[samplename] = ""

    # Open the haplo_file using gzip and iterate over each line
    with gzip.open(haplo_file, 'rt') as f:
        next(f)  # Skip the header line
        for line in f:
            splitted = line.strip().split("\t")
            # Splits each line by tabs and extracts chrom, pos, and alleles.
            # Allele is extracted from column 3
            chrom, pos, allele = splitted[0], splitted[1], splitted[3]
            # Appends the respective allele (character) to the sample's sequence in sample_dict.
            sample_dict[samplename] += allele

    # Print the sequences for each sample
    for key, value in sample_dict.items():
        print(">" + key)
        print(value)

if __name__ == "__main__":
    # Get filename and sample name from command line arguments
    filename = argv[1]
    sample = argv[2]
    # Call the parse_haplo function with filename and sample name
    parse_haplo(filename, sample)
