#!/usr/bin/env python

### Author: Bilal Sharif

from sys import argv
import matplotlib.pyplot as plt

if len(argv) != 4:
    print("Usage: plot_preseq.py <preseq_output_file> <steps_to_plot> <output_file>")
    exit(1)


filein = argv[1]
try:
    steps_to_plot = int(argv[2])
except ValueError:
    print("Error: steps_to_plot must be an integer")
    exit(1)
fileout = argv[3]

reads = []
uniqs = []

i = 0
with open(filein, 'r') as f1:
    next(f1)
    for line in f1:
        fields = line.strip().split()
        reads.append(int(float(fields[0])))
        uniqs.append(float(fields[1]))
        i += 1
        if i > steps_to_plot: ## beak after steps_to_plot to make the plot more readable
            break

plt.plot(reads, uniqs)
plt.xlabel('Number of Reads')
plt.ylabel('Unique Reads')
plt.title('Number of Unique Reads vs Number of Reads')
plt.savefig(fileout)
