#!/usr/bin/env python
from __future__ import division
from Bio import SeqIO
from sys import argv

"""
Author: Tom van der Valk
tom.vandervalk@ebc.uu.se

Script modifications: Verena Kutschera, Bilal Sharif

Script identifies CpG sites in a fasta reference
output is in BED-format
"""
fasta_in= argv[1]
bed_out= argv[2]

def find_cpg(reference,outputname):
    outputfile = open(outputname, "w")
    fasta_sequences = SeqIO.parse(open(reference),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq).upper()
        for i in range(len(sequence)-1):
            position = i
            if (sequence[i] + sequence[i+1]) == "CG":
                outputfile.write(name + "\t" + str(position) + "\t" + str(position + 2) + "\n")

    outputfile.close()

if __name__ == "__main__":
    find_cpg(fasta_in, bed_out)