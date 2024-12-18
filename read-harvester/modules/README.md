# Process definitions for tools used in the workflow

Implemented modules:

- samplesheet_check
- fastp
- bwa index
- bwa aln
- bwa samse
- FastQC
- MultiQC
- MapDamage2
- amber
- amber create_amber_samplesheet
- samremovedup
- samtools faidx
- samtools merge
- samtools index

Planned modules:

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