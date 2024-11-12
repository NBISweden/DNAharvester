#!/usr/bin/env python3

import sys
from statistics import mean

"""
Author:		Bilal Sharif
Contact: 	bilal.bioinfo@gmail.com
Usage:      read_len_cutoff_amber.py <amber_output.txt> <tolerance>
"""

if len(sys.argv) != 3:
    print("Usage: read_len_cutoff_amber.py <amber_output.txt>")
    sys.exit(1)

if float(sys.argv[2]) <= 0 or float(sys.argv[2]) >= 1:
    print("Error: The tolerance value should be a float between 0 and 1")
    sys.exit(1)

### Setting up initial variables
filein = sys.argv[1]
read_lengths = []
mismatch_rates = []
tolerance = 1 + float(sys.argv[2])



### Read the mismatch rates from the input file
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




### Estimate the read length cutoff based on the average mismatch rate and tolerance threshold

## Calculate the average mismatch rate using the mismatch rates of reads with lengths between 40 and 60 bp
avg_mismatch_rate = mean([mismatch_rates[read_lengths.index(i)] for i in range(40, 61)])

cutoff_length = None
warning = False
for i in range(39, read_lengths[0]-1, -1):  # Walk backward from 39 to the smallest read length
    mismatch_rate = mismatch_rates[read_lengths.index(i)]
    tolerance_threshold = avg_mismatch_rate * tolerance
    if mismatch_rate > tolerance_threshold:
        cutoff_length = i+1 ## the cutoff length is the read length of the last read with a mismatch rate below the tolerance threshold
        if not mismatch_rates[read_lengths.index(i-1)] > tolerance_threshold: ## if the mismatch rate of the previous read is not higher than the tolerance threshold, probably not necessary!!!!
            warning = True
        break
    else:
        ## updating the average mismatch rate and tolerance threshold
        avg_mismatch_rate = mean([avg_mismatch_rate, mismatch_rate])

## Output the result
if cutoff_length:
    print(f"Selected read length cutoff: {cutoff_length}")
    if warning:
        print(f"\nWARNING: The mismatch rate of read length {i} is higher than the tolerance threshold, but the mismatch rate of read length {i-1} is not. You might need to second look at the plot")
else:
    print(f"No significant increase in mismatch rate was found until read length {read_lengths[0]-1}. Check if short reads are already filtered out")


