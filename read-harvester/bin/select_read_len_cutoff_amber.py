#!/usr/bin/env python3

import sys
from statistics import mean

"""
Author:		Bilal Sharif
Contact: 	bilal.bioinfo@gmail.com
Usage:      read_len_cutoff_amber.py <amber_output.txt> <tolerance>
"""

### Checking the input arguments
if len(sys.argv) != 3:
    print("Usage: read_len_cutoff_amber.py <amber_output.txt> <tolerance>")
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

###### Estimate the read length cutoff based on the average mismatch rate and tolerance threshold

# If the shortest read length is >= 40, use it as the cutoff
if read_lengths[0] >= 40:
    print(f"Selected read length cutoff: {read_lengths[0]}")
    print("The shortest read length is already >= 40. No walk-back step required.")
    sys.exit(0)


## Calculate the average mismatch rate using lengths between 40 and 60 if available
filtered_mismatch_rates = [mismatch_rates[read_lengths.index(i)] for i in range(40, 61) if i in read_lengths]

if not filtered_mismatch_rates:
    print("Error: No mismatch rate data available for lengths between 40 and 60.")
    sys.exit(1)

avg_mismatch_rate = mean(filtered_mismatch_rates)


### Walk back from 39 to the smallest read length
cutoff_length = None
warning = False
for i in range(39, read_lengths[0] - 1, -1):  # Walk backward from 39 to the smallest read length
    if i not in read_lengths:
        continue  # Skip missing lengths
    mismatch_rate = mismatch_rates[read_lengths.index(i)]
    tolerance_threshold = avg_mismatch_rate * tolerance
    if mismatch_rate > tolerance_threshold: ## only setting higher here as it works best with bwa aligner - not for bowtie2
        cutoff_length = i + 1  # The cutoff length is the last read length with a mismatch rate below the tolerance threshold
        if i - 1 in read_lengths and not mismatch_rates[read_lengths.index(i - 1)] > tolerance_threshold:
            warning = True
        break
    else:
        # Update the average mismatch rate and tolerance threshold
        avg_mismatch_rate = mean([avg_mismatch_rate, mismatch_rate])

# Output the result
if cutoff_length:
    print(f"Selected read length cutoff: {cutoff_length}")
    if warning:
        print(
            f"\nWARNING: The mismatch rate of read length {i} is higher than the tolerance threshold, "
            f"but the mismatch rate of read length {i-1} is not. You might wamt to visually inspect the amber plot."
            )
else:
    print(
        f"Selected read length cutoff: {read_lengths[0]}\n"
        f"{read_lengths[0]} is also the shortest read in the file. "
        f"Check if short reads are already filtered out because no incease in mismatch rate was observed."
    )

