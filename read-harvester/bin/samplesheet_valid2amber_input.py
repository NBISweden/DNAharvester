#!/usr/bin/env python


"""
Provide a command line tool to transform a tabular samplesheet 
to an input file listing a bam file name for AMBER to be used
as input in a Nextflow pipeline, i.e. without file path.

Requires a validated tabular samplesheet as processed by 
check_samplesheet.py and the bam file name.
"""

import csv
import sys


def convert_samplesheet(file_in, bam):
    """
    Convert the validated, tabular samplesheet for nf-core pipelines into one tab-delimited 
    samplesheet per sample for AMBER.

    Args:
        file_in: Tabular samplesheet that was previously processed by check_samplesheet.py
        bam: BAM file for which the AMBER samplesheet is being generated
    """
    with open(file_in, "r") as in_handle:
        reader = csv.DictReader(in_handle)
        fieldnames = ["sample", "bam"]
        for row in reader:
            if row["sample"] + ".bam" == str(bam):
                file_out = row["sample"] + ".amber.tsv"
                with open(file_out, "w") as out_handle:
                    writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
                    new_row = {"sample": row["sample"], "bam": row["sample"] + ".bam"}
                    writer.writerow(new_row)

convert_samplesheet(sys.argv[1], sys.argv[2])