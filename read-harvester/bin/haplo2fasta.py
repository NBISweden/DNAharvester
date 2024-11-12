#!/usr/bin/env python3

from sys import argv
from itertools import islice
import gzip

"""
This script parses a *.haplo.gz file produced with ANGSD 
doHaploCall for one sample and converts it to FASTA format. 
"""

def parse_haplo(haplo_file, fasta_fai):
    # Open the fasta_fai and store chromosome names and lengths in a dictionary
    chr_dict = {}
    with open(fasta_fai) as fai:
        for line in fai:
            splitted = line.strip().split("\t")
            chrom, length = splitted[0], splitted[1]
            if chrom not in chr_dict:
                chr_dict[chrom] = length
    # Create an empty dictionary to store the sequence with a chromosome/scaffold name
    sequence_dict = {}
    # Open the haplo_file using gzip and iterate over each line
    with gzip.open(haplo_file, 'rt') as f:
        next(f)  # Skip the header line
        for line in f:
            splitted = line.strip().split("\t")
            # Splits each line by tabs and extracts chromosome/scaffold name and allele for ind0.
            # Allele is extracted from column 4
            chrom, pos, allele = splitted[0], splitted[1], splitted[3]
            if chrom not in sequence_dict:
                sequence_dict[chrom] = ""
            # Appends the respective allele (character) to the chromosome's/scaffold's sequence in sequence_dict.
            sequence_dict[chrom] += allele

    # Print the sequences for each chromosome/scaffold
    for key, value in sequence_dict.items():
        print(">" + key)
        print(value)

if __name__ == "__main__":
    # Get *haplo.gz and *.fasta.fai filenames from command line arguments
    filename = argv[1]
    fastaindex = argv[2]
    # Call the parse_haplo function with filename
    parse_haplo(filename, fastaindex)
