#! /usr/bin/env nextflow

include { BOWTIE2_BUILD as TC_BOWTIE2_BUILD                 } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BOWTIE2 as TC_BOWTIE2                             } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_MERGE as TC_MERGED_UNMERGED_READS_BAM    } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX as TC_RAW_BAM_INDEX                } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW_MQ as TC_SAMTOOLS_VIEW_MQ           } from '../../../modules/local/samtools/samtools_view_mq.nf'
include { SAMTOOLS_MERGE as TC_SAMTOOLS_MERGE_SAMPLE        } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMREMOVEDUP as TC_SAMREMOVEDUP_SAMPLE            } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as TC_SAMREMOVEDUP_SAMPLE_INDEX    } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS as TC_SAMTOOLS_IDXSTATS         } from '../../../modules/local/samtools/samtools_idxstats.nf'
include { SORT_IDXSTATS as TC_SORT_IDXSTATS                 } from '../../../modules/local/samtools/sort_idxstats.nf'

workflow TAXONOMIC_CLASSIFICATION {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads
    workflow_name

    main:
    ch_versions         = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Index the reference database if it is not already indexed
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Build the Bowtie2 index
    TC_BOWTIE2_BUILD (reference, file(params.tc_reference_database).getParent())
    ch_versions         = ch_versions.mix(TC_BOWTIE2_BUILD.out.versions)
    ch_bowtie2_index    = TC_BOWTIE2_BUILD.out.index_dir

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Mapping the reads to the reference database
    ////////////////////////////////////////////////////////////////////////////////////////////////

    TC_BOWTIE2 ( reads, ch_bowtie2_index )
    ch_versions         = ch_versions.mix(TC_BOWTIE2.out.versions)
    ch_raw_bam          = TC_BOWTIE2.out.bam

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

        TC_MERGED_UNMERGED_READS_BAM ( ch_raw_bam_grouped, reference )
        ch_versions = ch_versions.mix(TC_MERGED_UNMERGED_READS_BAM.out.versions)
        ch_merged_raw_bam = TC_MERGED_UNMERGED_READS_BAM.out.bam
    }
    // Use merged BAM if available, otherwise use raw BAM
    ch_raw_bam_for_processing = ch_merged_raw_bam ?: ch_raw_bam

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 4. BAM processing
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Index the BAM file
    TC_RAW_BAM_INDEX ( ch_raw_bam_for_processing )
    ch_versions             = ch_versions.mix ( TC_RAW_BAM_INDEX.out.versions )

    // Filter the BAM file by mapping quality
    TC_SAMTOOLS_VIEW_MQ ( ch_raw_bam_for_processing )
    ch_versions         = ch_versions.mix(TC_SAMTOOLS_VIEW_MQ.out.versions)


    // Prepare BAM files for merging per sample
    ch_bam_sample_to_merge = TC_SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] // remove library_id and lane for merging all bams per sample_id
        new_meta.remove('single_end') // remove single_end info from meta as the same sample can have both single-end and paired-end data
        new_meta.remove('read_group') // remove read_group info from meta as the same sample can have multiple read groups
        new_meta.remove('library_type') // remove library_type info from meta as the same sample can have both double-stranded and single-stranded libraries
        [ new_meta, bam ]
    }.groupTuple()

    // Merge BAM files per sample
    TC_SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, reference )
    ch_versions         = ch_versions.mix(TC_SAMTOOLS_MERGE_SAMPLE.out.versions)

    // Remove duplicates from the merged BAM files per sample
    TC_SAMREMOVEDUP_SAMPLE ( TC_SAMTOOLS_MERGE_SAMPLE.out.bam, reference )
    ch_versions         = ch_versions.mix(TC_SAMREMOVEDUP_SAMPLE.out.versions)
    // Index the BAM file after removing duplicates
    TC_SAMREMOVEDUP_SAMPLE_INDEX ( TC_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions         = ch_versions.mix(TC_SAMREMOVEDUP_SAMPLE_INDEX.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 5. Generate idxstats
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Generate idxstats for the final BAM file
    TC_SAMTOOLS_IDXSTATS ( TC_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions         = ch_versions.mix(TC_SAMTOOLS_IDXSTATS.out.versions)

    // Sort idxstats files
    TC_SORT_IDXSTATS ( TC_SAMTOOLS_IDXSTATS.out.idxstats )
    ch_versions         = ch_versions.mix(TC_SORT_IDXSTATS.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    raw_bam                 = ch_raw_bam_for_processing                   // channel: [ val(meta), [ bam ] ]
    raw_bam_bai             = TC_RAW_BAM_INDEX.out.bai                    // channel: [ val(meta), [ raw_bam_bai ] ]
    mq_filtered_bam         = TC_SAMTOOLS_VIEW_MQ.out.bam                 // channel: [ val(meta), [ mq_filtered_bam ] ]
    merged_bam_sample       = TC_SAMTOOLS_MERGE_SAMPLE.out.bam            // channel: [ val(meta), [ merged_bam_sample ] ]
    dedup_bam_sample        = TC_SAMREMOVEDUP_SAMPLE.out.dedup            // channel: [ val(meta), [ dedup_bam_sample ] ]
    dedup_bam_sample_bai    = TC_SAMREMOVEDUP_SAMPLE_INDEX.out.bai        // channel: [ val(meta), [ dedup_bam_sample_bai ] ]
    idxstats                = TC_SORT_IDXSTATS.out.sorted_idxstats        // channel: [ val(meta), [ idxstats ] ]
    versions                = ch_versions                                 // channel: [ versions.yml ]
}