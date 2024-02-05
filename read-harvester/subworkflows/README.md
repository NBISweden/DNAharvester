# Custom workflows for different stages of the main analysis

## Implemented subworkflows:

merge_filter_reads
- fastp (paired-end and single-strand libraries, minimum read length from config file or from subworkflow min_read_length_identification)

mapping
- bwa index
- bwa aln (aDNA specific parameters, paired-end and single-strand libraries, config file)
- bwa samse

raw_processed_mapped_reads_qc
- FastQC
- MapDamage2
- AMBER
- MultiQC (including fastp reports)

merge_dedup_realign_bams
- samtools faidx
- samtools merge to merge bam files per library index (sample_index_lane.bam --> sample_index.bam)
- samremovedup on merged bam files
- samtools merge to merge bam files per sample (sample_index.bam --> sample.bam)
- samremovedup on merged bam files
- samtools index
- picard createsequencedictionary
- GATK realignertargetcreator
- GATK indelrealigner

## Planned subworkflows:

bam_filtering
- samtools view -q minMQ (config file)

min_read_length_identification
- Custom code: identify min. read length from AMBER output (mismatch/read length plot)

low_depth_variant_calling
- bcftools mpileup (generate mpileup file)
- pseudohaploidization (draw 1 random read per site from mpileup file)
- bcftools call (pseudohaploid mpileup file > VCF file)

variant_filtering
- Reference genome bed file
- Reference genome genomefile (for BEDtools)
- CpG-site identification from reference genome (USER and non-USER-treated data)
- CpG-sites2bed (USER and non-USER-treated data)
- C>T transition identification (non-USER-treated data)
- C>T transitions2bed (non-USER-treated data)
- G>A transition identification (non-USER-treated and paired-end data)
- G>A transitions2bed (non-USER-treated and paired-end data)
- Repeatmodeler?
- Repeatmasker?
- Repeats2bed?
- BEDtools intersect (CpG, transition [and Repeat-]bed files vs. vcf files)

competetive_mapping
- Map to concatenated genome (reference + taxa identified to be present in sample via eMeta pipeline)

report
- Damage at read ends (prior to CpG filtering)
- Depurination-induced fragmentation patterns (non-USER treated data)
- Summary tables and figures of key data and results