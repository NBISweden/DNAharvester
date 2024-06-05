#!/usr/bin/python
from sys import argv
from itertools import izip, islice
import gzip

def parse_haplo(haplo_file):

  samples = ["SF-1_42","AI-1_119","BL-11_9","KU-23_3",
        "MO34_ZA","HD-3","RY-3","EL012","Sheep","domesticgoat"]

  sample_dict	= {}
  for	i in samples:
    sample_dict[i] = ""

  with gzip.open(haplo_file) as f1:
    next(f1)
    for line in f1:
      splitted = line.strip().split("\t")
      chrom,pos,alleles = splitted[0],splitted[1],"".join(splitted[3:])
      if "N" not in alleles:
        if len(set(alleles)) == 2:
          for	i in range(len(alleles)):
            focal_sample = samples[i]
            sample_dict[focal_sample] += alleles[i]

  for key,value in sample_dict.items():
    print ">" + key
    print value

if __name__ == "__main__":
  filename = argv[1]
  parse_haplo(filename)