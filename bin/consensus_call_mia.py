#!/usr/bin/env python

import sys
import argparse
import os

'''
Script taken from https://github.com/aersoares81/mia-helper-scripts
Original author: André Soares
License: MIT License
'''

def parse_arguments():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(description='This script will create a consensus FASTA file from the output of MIA.'
    'It will also filter sequences based on coverage, percentage agreement, and quality. Please create a .41 file using ma.'
    'For more details, check https://github.com/aersoares81/mia-helper-scripts')
    parser.add_argument('-c', '--coverage', type=int, default=3, help='Coverage cutoff; default = 3')
    parser.add_argument('-I', '--id', type=str, default=None, help='ID to assign; default is the name of the output FASTA file (without the .fasta extension)')
    parser.add_argument('-p', '--percent', type=float, default=0.66, help='Percentage of agreement of the base called; default = 0.66')
    parser.add_argument('-q', '--quality', type=int, default=40, help='Aggregate score quality cutoff; default = 40')
    parser.add_argument('-m', '--input', type=str, required=True, help='Input .41 file')
    parser.add_argument('-o', '--output', type=str, default=None, help='Output .fasta file')
    return parser.parse_args()

def filter_sequences(args):
    """
    Filters sequences based on coverage, percentage agreement, and quality.
    """
    seq = ''
    with open(args.input, 'r') as file:
        for line in file:
            el = line.strip().split()
            if el[2] != '-':
                if int(el[3]) >= args.coverage and float(el[14]) >= args.percent and int(el[13]) >= args.quality:
                    seq += el[2]
                else:
                    seq += 'N'
    return seq

def print_fasta(seq_id, seq, output_file):
    """
    Creates the output FASTA file.
    """
    with open(output_file, 'w') as file:
        file.write('>' + seq_id + '\n')
        for i in range(0, len(seq), 60):
            file.write(seq[i:i+60] + '\n')

def main():
    args = parse_arguments()
    seq = filter_sequences(args)

    if args.output is None:
        filename, _ = os.path.splitext(args.input)
        if ".maln" in filename:
            filename = filename.split(".maln", 1)[0]
        args.output = filename + ".fasta"

    if args.id is None:
        args.id = os.path.splitext(os.path.basename(args.output))[0]

    print_fasta(args.id, seq, args.output)

if __name__ == "__main__":
    main()