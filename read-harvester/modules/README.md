# Process definitions for tools used in the workflow

Currently planned modules:

- Fastp
- bwa index
- bwa aln
- bwa samse
- FastQC
- MultiQC
- MapDamage2
- AMBER
- samtools view -q minMQ
- custom code: identify minimum read length from AMBER output
- bcftools mpileup
- custom code: pseudohaploidization
- bcftools call
- Reference genome bed file
- Reference genome genomefile (bedtools format)
- CpG-site identification from reference genome
- CpG-site2bed
- custom code?: C>T transition identification
- C>T transitions2bed
- custom code?: G>A transition identification
- G>A transitions2bed
- Repeatmodeler?
- Repeatmasker?
- repeats2bed?
- bedtools intersect (bed files vs. VCF files)
- report: damage at read ends
- report: depurination-induced fragmentation patterns
- report: summary tables and figures of key data and results
- custom code for competetive mapping subworkflow