#!/usr/bin/env python


"""
Creates an input file listing a bam file name for AMBER to
as input in a Nextflow pipeline, i.e. without file path.

Requires the bam file name.
"""

import csv
import sys


def create_amber_samplesheet(bam):
    """
    Create one tab-delimited samplesheet per BAM file for AMBER.

    Args:
        bam: BAM file for which the AMBER samplesheet is being generated
    """
    with open(file_out, "w") as out_handle:
        sample = str(bam).strip(".bam")
        writer = csv.DictWriter(out_handle, fieldnames=["sample", "bam"], delimiter="\t")
        new_row = {"sample": sample, "bam": str(bam)}
        writer.writerow(new_row)

create_amber_samplesheet(sys.argv[1])