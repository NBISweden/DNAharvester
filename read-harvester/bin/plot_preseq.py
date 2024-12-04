#!/usr/bin/env python

### Author: Bilal Sharif

from sys import argv
import matplotlib.pyplot as plt

filein = argv[1]
fileout = argv[2]

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
        if i > 250: ## plotting only the first 250 points to avoid overplotting. increase this number if you want to plot more points
            break

plt.plot(reads, uniqs)
plt.xlabel('Number of Reads')
plt.ylabel('Unique Reads')
plt.title('Number of Unique Reads vs Number of Reads')
plt.savefig(fileout)
