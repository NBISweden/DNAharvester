#! /usr/bin/env nextflow

include { BOWTIE2_BUILD                                     } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BOWTIE2 as MS_BOWTIE2                             } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_INDEX as MS_BOWTIE2_INDEX                } from '../../../modules/nf-core/samtools/index/main'

// Mapping quality filter
include { SAMTOOLS_VIEW_MQ as MS_SAMTOOLS_VIEW_MQ           } from '../../../modules/local/samtools/samtools_view_mq.nf'
include { SAMTOOLS_INDEX as MS_SAMTOOLS_VIEW_MQ_INDEX       } from '../../../modules/nf-core/samtools/index/main'

// Merge BAM files per library/PCR
include { SAMTOOLS_MERGE as MS_SAMTOOLS_MERGE_LIB           } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX as MS_SAMTOOLS_MERGE_LIB_INDEX     } from '../../../modules/nf-core/samtools/index/main'

// Remove duplicates from BAM files merged per library/PCR
include { SAMREMOVEDUP as MS_SAMREMOVEDUP_LIB               } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as MS_SAMREMOVEDUP_LIB_INDEX       } from '../../../modules/nf-core/samtools/index/main'

// Merge BAM files per sample
include { SAMTOOLS_MERGE as MS_SAMTOOLS_MERGE_SAMPLE        } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX as MS_SAMTOOLS_MERGE_SAMPLE_INDEX  } from '../../../modules/nf-core/samtools/index/main'

// Remove duplicates from BAM files merged per sample
include { SAMREMOVEDUP as MS_SAMREMOVEDUP_SAMPLE            } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as MS_SAMREMOVEDUP_SAMPLE_INDEX    } from '../../../modules/nf-core/samtools/index/main'

// SAMTOOLS_IDXSTATS
include { SAMTOOLS_IDXSTATS as MS_SAMTOOLS_IDXSTATS         } from '../../../modules/local/samtools/samtools_idxstats.nf'


workflow MYSTERY_SAMPLE {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads

    main:
    ch_versions         = Channel.empty()


    // Index the reference genome if it is not already indexed
    BOWTIE2_BUILD (reference, file(params.ms_reference_database).getParent())
    ch_versions             = ch_versions.mix(BOWTIE2_BUILD.out.versions)


    // Map the reads to the reference genome
    MS_BOWTIE2 ( reads, BOWTIE2_BUILD.out.index_dir )
    ch_versions             = ch_versions.mix(MS_BOWTIE2.out.versions)
    ch_bam                  = MS_BOWTIE2.out.bam
    // Index the BAM file
    MS_BOWTIE2_INDEX ( ch_bam )
    ch_versions             = ch_versions.mix ( MS_BOWTIE2_INDEX.out.versions )


    // Filter the BAM file by mapping quality
    MS_SAMTOOLS_VIEW_MQ ( ch_bam )
    ch_versions         = ch_versions.mix(MS_SAMTOOLS_VIEW_MQ.out.versions)
    // Index the filtered BAM file
    MS_SAMTOOLS_VIEW_MQ_INDEX ( MS_SAMTOOLS_VIEW_MQ.out.bam )
    ch_versions         = ch_versions.mix(MS_SAMTOOLS_VIEW_MQ_INDEX.out.versions)


    // Prepare the BAM files for merging by updating the metadata
    ch_bam_lib_to_merge = MS_SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
        // update only the 'id' field in meta, keep all other fields
        [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1], 'library_type': meta.library_type, 'single_end': meta.single_end], bam]
    }.groupTuple()


    // Merge BAM files per library/PCR
    MS_SAMTOOLS_MERGE_LIB ( ch_bam_lib_to_merge, reference )
    ch_versions         = ch_versions.mix(MS_SAMTOOLS_MERGE_LIB.out.versions)
    // Index the merged BAM file per library/PCR
    MS_SAMTOOLS_MERGE_LIB_INDEX ( MS_SAMTOOLS_MERGE_LIB.out.bam )
    ch_versions         = ch_versions.mix(MS_SAMTOOLS_MERGE_LIB_INDEX.out.versions)


    // Remove duplicates from the merged BAM files per library/PCR
    MS_SAMREMOVEDUP_LIB ( MS_SAMTOOLS_MERGE_LIB.out.bam, reference )
    ch_versions         = ch_versions.mix(MS_SAMREMOVEDUP_LIB.out.versions)
    // Index the BAM file after removing duplicates
    MS_SAMREMOVEDUP_LIB_INDEX ( MS_SAMREMOVEDUP_LIB.out.dedup )
    ch_versions         = ch_versions.mix(MS_SAMREMOVEDUP_LIB_INDEX.out.versions)


    // Prepare the BAM files for merging by updating the metadata
    ch_bam_sample_to_merge = MS_SAMREMOVEDUP_LIB.out.dedup.map { meta, bam ->
            // update only the 'id' field in meta, keep all other fields
            [['id': meta.id.split("_")[0]], bam]
        }.groupTuple()


    // Merge BAM files per sample
    MS_SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, reference )
    ch_versions         = ch_versions.mix(MS_SAMTOOLS_MERGE_SAMPLE.out.versions)
    // Index the merged BAM file per sample
    MS_SAMTOOLS_MERGE_SAMPLE_INDEX ( MS_SAMTOOLS_MERGE_SAMPLE.out.bam )
    ch_versions         = ch_versions.mix(MS_SAMTOOLS_MERGE_SAMPLE_INDEX.out.versions)


    // Remove duplicates from the merged BAM files per sample
    MS_SAMREMOVEDUP_SAMPLE ( MS_SAMTOOLS_MERGE_SAMPLE.out.bam, reference )
    ch_versions         = ch_versions.mix(MS_SAMREMOVEDUP_SAMPLE.out.versions)
    // Index the BAM file after removing duplicates
    MS_SAMREMOVEDUP_SAMPLE_INDEX ( MS_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions         = ch_versions.mix(MS_SAMREMOVEDUP_SAMPLE_INDEX.out.versions)


    // Generate idxstats for the final BAM file
    MS_SAMTOOLS_IDXSTATS ( MS_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions         = ch_versions.mix(MS_SAMTOOLS_IDXSTATS.out.versions)


    emit:
    reference_index         = BOWTIE2_BUILD.out.index_dir                 // channel: path(index)
    raw_bam                 = ch_bam                                      // channel: [ val(meta), [ bam ] ]
    raw_bam_bai             = MS_BOWTIE2_INDEX.out.bai                    // channel: [ val(meta), [ raw_bam_bai ] ]
    mq_filtered_bam         = MS_SAMTOOLS_VIEW_MQ.out.bam                 // channel: [ val(meta), [ mq_filtered_bam ] ]
    mq_filtered_bam_bai     = MS_SAMTOOLS_VIEW_MQ_INDEX.out.bai           // channel: [ val(meta), [ mq_filtered_bam_bai ] ]
    merged_bam_lib          = MS_SAMTOOLS_MERGE_LIB.out.bam               // channel: [ val(meta), [ merged_bam_lib ] ]
    merged_bam_lib_bai      = MS_SAMTOOLS_MERGE_LIB_INDEX.out.bai         // channel: [ val(meta), [ merged_bam_lib_bai ] ]
    dedup_bam_lib           = MS_SAMREMOVEDUP_LIB.out.dedup               // channel: [ val(meta), [ dedup_bam_lib ] ]
    dedup_bam_lib_bai       = MS_SAMREMOVEDUP_LIB_INDEX.out.bai           // channel: [ val(meta), [ dedup_bam_lib_bai ] ]
    merged_bam_sample       = MS_SAMTOOLS_MERGE_SAMPLE.out.bam            // channel: [ val(meta), [ merged_bam_sample ] ]
    merged_bam_sample_bai   = MS_SAMTOOLS_MERGE_SAMPLE_INDEX.out.bai      // channel: [ val(meta), [ merged_bam_sample_bai ] ]
    dedup_bam_sample        = MS_SAMREMOVEDUP_SAMPLE.out.dedup            // channel: [ val(meta), [ dedup_bam_sample ] ]
    dedup_bam_sample_bai    = MS_SAMREMOVEDUP_SAMPLE_INDEX.out.bai        // channel: [ val(meta), [ dedup_bam_sample_bai ] ]
    idxstats                = MS_SAMTOOLS_IDXSTATS.out.idxstats           // channel: [ val(meta), [ idxstats ] ]
    versions                = ch_versions                                 // channel: [ versions.yml ]
}