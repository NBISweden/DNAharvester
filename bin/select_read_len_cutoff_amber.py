#!/usr/bin/env python3

import sys
from statistics import mean
from kneed import KneeLocator

"""
Author:		    Bilal Sharif
Contact: 	    bilal.bioinfo@gmail.com
Usage:          read_len_cutoff_amber.py <amber_output.txt> <curve> <direction>
Description:    This script reads mismatch rates from an AMBER output file and determines the read length cutoff
                using knee location.
"""


### Checking the input arguments
if len(sys.argv) != 4:
    print("Usage: read_len_cutoff_amber.py <amber_output.txt> <curve> <direction>")
    sys.exit(1)


### Setting up initial variables
filein = sys.argv[1]
read_lengths = []
mismatch_rates = []
curve = sys.argv[2]
direction = sys.argv[3]


### Check if the curve is valid
if curve not in ["convex", "concave"]:
    print("Error: Invalid curve type. Please use 'convex' or 'concave'.")
    sys.exit(1)
### Check if the direction is valid
if direction not in ["increasing", "decreasing"]:
    print("Error: Invalid direction. Please use 'increasing' or 'decreasing'.")
    sys.exit(1)


############### Read the read length and mismatch rate from the input file ###############
data = False
with open(filein, "r") as f1:
    for line in f1:
        if "MISMATCH_RATE" in line:
            data = True
            continue  ## skip the header line
        if data:
            if line.startswith("-"):
                break
            line = line.strip()
            columns = line.split("\t")
            try:
                read_lengths.append(int(columns[0]))
                mismatch_rates.append(float(columns[1]))
            except:
                print("Error: Unexpected data format in the input file")
                sys.exit(1)
if not data:
    print("Error: No mismatch rate data found in the input file. Please check the input file")
    sys.exit(1)


############### Determine the read length cutoff using knee locator ###############

## subset the read lengths and mismatch rates to those less than 40
subset_indices = [i for i, rl in enumerate(read_lengths) if rl < 40]
read_lengths_sub = [read_lengths[i] for i in subset_indices]
mismatch_rates_sub = [mismatch_rates[i] for i in subset_indices]

# Sort both lists in decreasing order of read length (walk backward)
sorted_pairs = sorted(zip(read_lengths_sub, mismatch_rates_sub), reverse=True)
read_lengths_sub_sorted, mismatch_rates_sub_sorted = zip(*sorted_pairs)

knee = KneeLocator(read_lengths_sub_sorted, mismatch_rates_sub_sorted, curve=curve, direction=direction)

# Check if a knee point was found
if knee.knee:
    print(f"Selected read length cutoff: {knee.knee}")
else:
    print("Error: No suitable read length cutoff found.")
    sys.exit(1)



