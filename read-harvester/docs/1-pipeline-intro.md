# Readharvester pipeline

Welcome to the readharvester documentation!

This Nextflow pipeline analyzes whole-genome sequencing data 
from highly degraded ancient DNA.

The pipeline takes raw FASTQ files with paired-end read data 
as input, merges them with fastp and maps them to a reference 
genome with bwa aln and ancient DNA-specific parameters. It 
runs AMBER on each BAM file and determines the minimum read 
length per sample at which reads can be confidently mapped to 
the reference by analysing the AMBER output. PCR duplicates 
are removed per PCR/library and from the merged BAM files. 
Reads are realigned and variants are randomly sampled to 
generate pseudo-haploid VCF and fasta files. 


## Contact

For bug reports, comments or suggestions, please open an issue here on GitHub.

> Please note that response times are longer around July and the Christmas holidays.

