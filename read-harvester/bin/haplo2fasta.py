#!/usr/bin/env python3

from sys import argv
from itertools import islice
import gzip
import pandas as pd

"""
This script parses a *.haplo.gz file produced with ANGSD 
doHaploCall for one sample and the *.fasta.fai file of the 
reference genome that was used to generate the *.haplo.gz 
file and converts the *.haplo.gz to FASTA format, filling 
any missing sites with "N". 
"""

def parse_haplo(haplo_file):
    # Open haplo_file using gzip and iterate over each line
    # Create a list to store the data from haplo_file
    haplo_data = []
    with gzip.open(haplo_file, 'rt') as f:
        next(f)  # Skip the header line
        for line in f:
            splitted = line.strip().split("\t")
            # Splits each line by tabs and extracts chromosome/scaffold name, position, and allele for ind0.
            # Allele is extracted from column 4
            chrom, position, allele = splitted[0], int(splitted[1]), splitted[3]
            # Append the data to the haplo_data list
            haplo_data.append([chrom, position, allele])
    # Convert the list to a pandas dataframe
    haplo_df = pd.DataFrame(haplo_data, columns=['chrom', 'position', 'allele'])
    return haplo_df

def parse_fai(fasta_fai):
    # Open the fasta_fai and store chromosome names and lengths in a dataframe
    # Create a list to store the chromosome names and reference genome positions
    fai_data = []
    with open(fasta_fai) as fai:
        for line in fai:
            splitted = line.strip().split("\t")
            chrom, length = splitted[0], int(splitted[1])
            for position in range(1, length + 1):
                fai_data.append([chrom, position])
    # Convert the list to a pandas dataframe
    fai_df = pd.DataFrame(fai_data, columns=['chrom', 'position'])
    return fai_df

def merge_fai_haplo_df(fai_df, haplo_df):
    # Merge the dataframes by chromosome name and position, using all positions from fai_df
    merged_df = fai_df.merge(haplo_df, on=['chrom', 'position'], how='left')
    # Fill any positions with a missing allele from haplo_df with N
    merged_df['allele'] = merged_df['allele'].fillna('N')
    return merged_df

def print_to_fasta(merged_df):
    # Print the merged_df in fasta format
    for chrom, group in merged_df.groupby('chrom'):
        alleles = ''.join(group['allele'])
        print(f">{chrom}\n{alleles}")

if __name__ == "__main__":
    # Get *haplo.gz and *.fasta.fai filenames from command line arguments
    haplo = argv[1]
    fai = argv[2]
    # Call the functions with haplo and fai
    fai_haplo_df = merge_fai_haplo_df(parse_fai(fai), parse_haplo(haplo))
    print_to_fasta(fai_haplo_df)