# Data-processing steps of the pipeline

## Subworkflows

The following subworkflows can be run subsequently to one 
another, so each subworkflow of the pipeline depends on the 
subworkflows before. 

### Input check

- Process the sample sheet
- Run FastQC and MultiQC on the raw fastq files

### Fastq processing

- Run FastP on the raw fastq files for paired-end read merging, 
default quality trimming and adapter removal (parameters: 
`--merge --correction --overlap_len_require 15 --overlap_diff_limit 1 --detect_adapter_for_pe`)

### Processed fastq quality check

- Run FastQC and MultiQC on the processed fastq files

### Mapping

- Run bwa index on the reference genome
- Run bwa aln with ancient DNA specific parameters to map  
the merged reads to the reference genome
- Convert SAI to BAM with bwa samse
- Index the BAM files with samtools index

### Raw BAM quality check

- Run samtools flagstat, mapDamage2 (without rescaling) and 
MultiQC on the BAM files 
- Run AMBER on a subset of mapped reads and identify a 
minimum read length from the AMBER output 

### BAM file processing

- Run samtools faidx on the reference genome
- Merge BAM files per PCA/library with samtools merge and 
index the merged BAM files with samtools index 
- Remove PCR duplicates with a custom script and 
index the deduplicated BAM files with samtools index 
- Merge the BAM files per sample with samtools merge and 
index the merged BAM files with samtools index 
- Remove PCR duplicates with a custom script and 
index the deduplicated BAM files with samtools index 
- Run Picard CreateSequenceDictionary on the reference genome 
- Run GATK RealignerTargetCreator and IndelRealigner 

### Processed BAM quality check

- Run samtools flagstat and MultiQC on the BAM files 
after merging per PCR/library 
- Run samtools flagstat and MultiQC on the BAM files 
after removing duplicates from each PCR/library BAM file 
- Run samtools flagstat and MultiQC on the BAM files 
after merging per sample 
- Run samtools flagstat and MultiQC on the BAM files 
after removing duplicates from each sample BAM file 
- Run samtools flagstat and MultiQC on the BAM files 
after indel realignment
- Run samtools depth to calculate mean depth across 
the reference genome per BAM file after indel realignment 

### Random sampling BAM

- Run ANGSD doHaploCall 1 on the processed BAM files to 
randomly sample one base per site 
- Convert the ANGSD `*.haplo.gz` files to plink format with 
ANGSD haploToPlink 
- Convert plink to VCF format with plink 
- Convert the ANGSD `*.haplo.gz` files to fasta format with 
a custom script, coding reference genome sites without data 
as "N" 
