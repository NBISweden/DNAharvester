#! /usr/bin/env nextflow

include { SAMTOOLS_UNMAPPED_READS                           } from '../../../modules/local/samtools/samtools_unmapped_reads.nf'
include { BWA_INDEX as MS_BWA_INDEX                         } from '../../../modules/local/bwa/bwa_index.nf'
include { BWA_ALN as MS_BWA_ALN                             } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_MEM as MS_BWA_MEM                             } from '../../../modules/local/bwa/bwa_mem.nf'
include { SPLIT_FASTQ as MS_SPLIT_FASTQ                     } from '../../../modules/local/awk/split_fastq.nf'
include { BWA_ALN as MS_BWA_ALN_SHORT                       } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_MEM as MS_BWA_MEM_LONG                        } from '../../../modules/local/bwa/bwa_mem.nf'
include { SAMTOOLS_MERGE as MS_BWA_ALN_MEM_MERGE            } from '../../../modules/local/samtools/samtools_merge.nf'
include { BOWTIE2_BUILD as MS_BOWTIE2_BUILD                 } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BOWTIE2 as MS_BOWTIE2                             } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_MERGE as MS_MERGED_UNMERGED_READS_BAM    } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX as MS_RAW_BAM_INDEX                } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW_MQ as MS_SAMTOOLS_VIEW_MQ           } from '../../../modules/local/samtools/samtools_view_mq.nf'
include { SAMTOOLS_MERGE as MS_SAMTOOLS_MERGE_SAMPLE        } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMREMOVEDUP as MS_SAMREMOVEDUP_SAMPLE            } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as MS_SAMREMOVEDUP_SAMPLE_INDEX    } from '../../../modules/nf-core/samtools/index/main'
include { FILTERBAM as MS_FILTERBAM                         } from '../../../modules/local/stats_output/filterbam.nf'
include { FILTERBAM_PLOT as MS_FILTERBAM_PLOT               } from '../../../modules/local/microbial_screening/filterbam_plot.nf'


workflow MICROBIAL_SCREENING {
    take:
    raw_bam
    ms_reference

    main:
    ch_versions             = Channel.empty()

    // get unmapped reads from the raw BAM file
    SAMTOOLS_UNMAPPED_READS ( raw_bam )
    ch_versions             = ch_versions.mix(SAMTOOLS_UNMAPPED_READS.out.versions)
    ch_unmapped_reads       = SAMTOOLS_UNMAPPED_READS.out.unmapped_fastq

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Index the reference database if it is not already indexed
    ////////////////////////////////////////////////////////////////////////////////////////////////

    if (params.ms_mapping_tool == 'bwa-aln' || params.ms_mapping_tool == 'bwa-mem' || params.ms_mapping_tool == 'bwa-aln-mem') {
        // Build the BWA index
        MS_BWA_INDEX (ms_reference, file(params.ms_reference_database).getParent())
        ch_versions             = ch_versions.mix(MS_BWA_INDEX.out.versions)
        ch_ms_reference_index   = MS_BWA_INDEX.out.index_dir
    } else if (params.ms_mapping_tool == 'bowtie2') {
        // Build the Bowtie2 index
        MS_BOWTIE2_BUILD (ms_reference, file(params.ms_reference_database).getParent())
        ch_versions             = ch_versions.mix(MS_BOWTIE2_BUILD.out.versions)
        ch_ms_reference_index   = MS_BOWTIE2_BUILD.out.index_dir
    } else {
        error "Invalid mapping tool specified for Microbial Screening: ${params.ms_mapping_tool}. Use 'bwa-aln', 'bwa-mem', 'bwa-aln-mem', or 'bowtie2'."
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Map the unmapped reads to the microbial reference database
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // BWA ALN
    if (params.ms_mapping_tool == 'bwa-aln' ) {
        MS_BWA_ALN ( ch_unmapped_reads, ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BWA_ALN.out.versions)
        ch_raw_bam          = MS_BWA_ALN.out.bam
    }
    // BWA MEM
    else if (params.ms_mapping_tool == 'bwa-mem') {
        MS_BWA_MEM ( ch_unmapped_reads, ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BWA_MEM.out.versions)
        ch_raw_bam          = MS_BWA_MEM.out.bam
    }
    else if (params.ms_mapping_tool == 'bwa-aln-mem') {
        //split fastq
        MS_SPLIT_FASTQ ( ch_unmapped_reads )
        ch_versions         = ch_versions.mix(MS_SPLIT_FASTQ.out.versions)
        //align short reads with BWA ALN
        MS_BWA_ALN_SHORT ( MS_SPLIT_FASTQ.out.short_reads, ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BWA_ALN_SHORT.out.versions)
        //align long reads with BWA MEM
        MS_BWA_MEM_LONG ( MS_SPLIT_FASTQ.out.long_reads,  ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BWA_MEM_LONG.out.versions)
        //merge BAMs from short and long reads
        ch_raw_bam_aln_mem  = MS_BWA_ALN_SHORT.out.bam.join(MS_BWA_MEM_LONG.out.bam)
            .map { meta, file1, file2 -> [meta, [file1, file2]] }
        MS_BWA_ALN_MEM_MERGE ( ch_raw_bam_aln_mem, reference )
        ch_versions         = ch_versions.mix(MS_BWA_ALN_MEM_MERGE.out.versions)
        ch_raw_bam          = MS_BWA_ALN_MEM_MERGE.out.bam
    }
    // BOWTIE2
    else if (params.ms_mapping_tool == 'bowtie2') {
        MS_BOWTIE2 ( ch_unmapped_reads, ch_ms_reference_index )
        ch_versions         = ch_versions.mix(MS_BOWTIE2.out.versions)
        ch_raw_bam          = MS_BOWTIE2.out.bam
    }
    else {
        error "Invalid mapping tool specified for Taxonomic Classification: ${params.ms_mapping_tool}. Use 'bwa-aln', 'bwa-mem', 'bwa-aln-mem' or 'bowtie2'."
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Merge the mapped unmerged reads if provided
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // if both merged reads and unmerged reads are mapped, merge the BAM files per sample
    def ch_merged_raw_bam = null
    if (params.merge_reads.toBoolean() && params.keep_unmerged_reads.toBoolean()) {
        // Group the unmerged reads BAMs by sample ID (removing the '-unmerged' suffix)
        ch_raw_bam_grouped = ch_raw_bam.map { meta, bam ->
                def new_meta = meta.clone()
                new_meta.id = meta.id.replace("-unmerged", "")
                tuple(new_meta.id, new_meta, bam)
            }
            .groupTuple()
            .map { id, metas, bams ->
                // pick one meta; choose the one with single_end == false if present
                def final_meta = metas.find { !it.single_end } ?: metas[0]
                // remove grouping id, return meta + joined bam list
                tuple(final_meta, bams)
            }

        MS_MERGED_UNMERGED_READS_BAM ( ch_raw_bam_grouped, reference )
        ch_versions = ch_versions.mix(MS_MERGED_UNMERGED_READS_BAM.out.versions)
        ch_merged_raw_bam = MS_MERGED_UNMERGED_READS_BAM.out.bam
    }
    // Use merged BAM if available, otherwise use raw BAM
    ch_raw_bam_for_processing = ch_merged_raw_bam ?: ch_raw_bam

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 4. BAM processing: index, filter, merge per sample, remove duplicates
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Index the BAM file
    MS_RAW_BAM_INDEX ( ch_raw_bam_for_processing )
    ch_versions             = ch_versions.mix ( MS_RAW_BAM_INDEX.out.versions )

    // mapping quality filter
    MS_SAMTOOLS_VIEW_MQ ( ch_raw_bam_for_processing )
    ch_versions             = ch_versions.mix ( MS_SAMTOOLS_VIEW_MQ.out.versions )

    // Prepare channel to merge BAM files per sample
    ch_ms_bams_per_sample  = MS_SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
        // update only the 'id' field in meta, keep all other fields
        [['id': meta.id.split("_")[0]], bam]
        }.groupTuple()

    // Merge BAM files per sample
    MS_SAMTOOLS_MERGE_SAMPLE ( ch_ms_bams_per_sample, ms_reference )
    ch_versions             = ch_versions.mix(MS_SAMTOOLS_MERGE_SAMPLE.out.versions)

    // Remove PCR duplicates
    MS_SAMREMOVEDUP_SAMPLE ( MS_SAMTOOLS_MERGE_SAMPLE.out.bam, ms_reference )
    ch_versions             = ch_versions.mix(MS_SAMREMOVEDUP_SAMPLE.out.versions)
    MS_SAMREMOVEDUP_SAMPLE_INDEX ( MS_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions             = ch_versions.mix(MS_SAMREMOVEDUP_SAMPLE_INDEX.out.versions)
    ch_bam_bai_final        = MS_SAMREMOVEDUP_SAMPLE.out.dedup.join(MS_SAMREMOVEDUP_SAMPLE_INDEX.out.bai)

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 5. Generate stats and plots
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // run filterBAM to generate stats
    MS_FILTERBAM ( ch_bam_bai_final )
    ch_versions             = ch_versions.mix ( MS_FILTERBAM.out.versions )

    // generate plot from filterBAM output
    MS_FILTERBAM_PLOT ( ch_bam_bai_final, MS_FILTERBAM.out.filterBAM_stats )
    ch_versions             = ch_versions.mix ( MS_FILTERBAM_PLOT.out.versions )

    ////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    unmapped_fastq             = SAMTOOLS_UNMAPPED_READS.out.unmapped_fastq  // channel: [ val(meta), [ fastq ] ]
    raw_bam                    = ch_raw_bam_for_processing                   // channel: [ val(meta), [ bam ] ]
    raw_bam_bai                = MS_RAW_BAM_INDEX.out.bai                    // channel: [ val(meta), [ bai ] ]
    mq_filtered_bam            = MS_SAMTOOLS_VIEW_MQ.out.bam                 // channel: [ val(meta), [ mq_filtered_bam ] ]
    merged_bam_sample          = MS_SAMTOOLS_MERGE_SAMPLE.out.bam            // channel: [ val(meta), [ merged_bam_sample ] ]
    dedup_bam_sample           = MS_SAMREMOVEDUP_SAMPLE.out.dedup            // channel: [ val(meta), [ dedup_bam_sample ] ]
    dedup_bam_sample_bai       = MS_SAMREMOVEDUP_SAMPLE_INDEX.out.bai        // channel: [ val(meta), [ dedup_bam_sample_bai ] ]
    filterBAM_stats            = MS_FILTERBAM.out.filterBAM_stats            // channel: [ val(meta), [ filterBAM.csv ] ]
    screening_plot             = MS_FILTERBAM_PLOT.out.ms_plot               // channel: [ val(meta), [ pathogen_screening_plot.pdf ] ]
    versions                   = ch_versions                                 // channel: [ versions.yml ]
}