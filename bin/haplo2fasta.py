#!/usr/bin/env python3

from sys import argv
import gzip
import pandas as pd

"""
This script parses a *.haplo.gz file produced with ANGSD 
doHaploCall and converts it to FASTA format, filling any 
missing sites with "N". 

Input is a *.haplo.gz file for one chromosome and one sample 
and the reference genome *.fai index file.
"""

def parse_haplo(haplo_file):
    # Open haplo_file using gzip and iterate over each line
    # Create a list to store the data from haplo_file
    haplo_data = []
    with gzip.open(haplo_file, 'rt') as f:
        for line in f:
            splitted = line.strip().split("\t")
            # Splits each line by tabs and extracts chromosome/scaffold name, position, and allele for ind0.
            # Allele is extracted from column 4
            chrom, position, allele = splitted[0], int(splitted[1]), splitted[3]
            if splitted[1] == "pos":
                continue
            else:
                # Append the data to the haplo_data list
                haplo_data.append([chrom, position, allele])
    # Convert the list to a pandas dataframe
    haplo_df = pd.DataFrame(haplo_data, columns=['chrom', 'position', 'allele'])
    return haplo_df

def parse_fai(fasta_fai):
    # Open the fasta_fai and store chromosome names and lengths in a dataframe
    fai_data = []
    with open(fasta_fai) as fai:
        for line in fai:
            splitted = line.strip().split("\t")
            chrom, length = splitted[0], int(splitted[1])
            fai_data.append([chrom, length])
    # Convert the list to a pandas dataframe
    fai_df = pd.DataFrame(fai_data, columns=['chrom', 'length'])
    return fai_df

def extract_chromosome_positions(haplo_df, fai_df):
    chrom_list = haplo_df['chrom'].unique().tolist()
    # Check that the *.haplo.gz file contained data from only one chromosome
    # and if yes, create a dataframe with chromosome name and all positions
    # on that chromosome
    if len(chrom_list) == 1:
        genomefile_data = []
        chrom = chrom_list[0]
        chrom_fai_df = fai_df[fai_df['chrom'].str.contains(chrom, na=False)]
        length = chrom_fai_df.iloc[0, 1]
        for position in range(1, length + 1):
            genomefile_data.append([chrom, position])
        # Convert the list to a pandas dataframe
        genomefile_df = pd.DataFrame(genomefile_data, columns=['chrom', 'position'])
    else:
        raise ValueError("The *.haplo.gz file contains data from more than one chromosome.\nPlease provide a *.haplo.gz file for one chromosome as input to haplo2fasta.py.")
    return genomefile_df

def merge_genomefile_haplo_df(genomefile_df, haplo_df):
    # Merge the dataframes by chromosome name and position, using all positions from genomefile_df
    merged_df = genomefile_df.merge(haplo_df, on=['chrom', 'position'], how='left')
    # Fill any positions with a missing allele from haplo_df with N
    merged_df['allele'] = merged_df['allele'].fillna('N')
    return merged_df

def print_to_fasta(merged_df):
    # Print the merged_df in fasta format
    for chrom, group in merged_df.groupby('chrom'):
        alleles = ''.join(group['allele'])
        print(f">{chrom}\n{alleles}")

if __name__ == "__main__":
    # Read in the *haplo.gz and *.fasta.fai files as dataframes
    haplo_df = parse_haplo(argv[1])
    fai_df = parse_fai(argv[2])
    # Call the functions with haplo_df and fai_df
    genomefile_df = extract_chromosome_positions(haplo_df, fai_df)
    fai_haplo_df = merge_genomefile_haplo_df(genomefile_df, haplo_df)
    print_to_fasta(fai_haplo_df)