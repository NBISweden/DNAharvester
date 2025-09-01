#! /usr/bin/env nextflow

include { BOWTIE2_BUILD                                     } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BOWTIE2 as SI_BOWTIE2                             } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_INDEX as SI_BOWTIE2_INDEX                } from '../../../modules/nf-core/samtools/index/main'

// Mapping quality filter
include { SAMTOOLS_VIEW_MQ as SI_SAMTOOLS_VIEW_MQ           } from '../../../modules/local/samtools/samtools_view_mq.nf'
include { SAMTOOLS_INDEX as SI_SAMTOOLS_VIEW_MQ_INDEX       } from '../../../modules/nf-core/samtools/index/main'

// Merge BAM files per library/PCR
include { SAMTOOLS_MERGE as SI_SAMTOOLS_MERGE_LIB           } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX as SI_SAMTOOLS_MERGE_LIB_INDEX     } from '../../../modules/nf-core/samtools/index/main'

// Remove duplicates from BAM files merged per library/PCR
include { SAMREMOVEDUP as SI_SAMREMOVEDUP_LIB               } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as SI_SAMREMOVEDUP_LIB_INDEX       } from '../../../modules/nf-core/samtools/index/main'

// Merge BAM files per sample
include { SAMTOOLS_MERGE as SI_SAMTOOLS_MERGE_SAMPLE        } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX as SI_SAMTOOLS_MERGE_SAMPLE_INDEX  } from '../../../modules/nf-core/samtools/index/main'

// Remove duplicates from BAM files merged per sample
include { SAMREMOVEDUP as SI_SAMREMOVEDUP_SAMPLE            } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as SI_SAMREMOVEDUP_SAMPLE_INDEX    } from '../../../modules/nf-core/samtools/index/main'

// SAMTOOLS_IDXSTATS
include { SAMTOOLS_IDXSTATS as SI_SAMTOOLS_IDXSTATS         } from '../../../modules/local/samtools/samtools_idxstats.nf'
include { MERGE_IDXSTATS as SI_MERGE_IDXSTATS               } from '../../../modules/local/merge_idxstats/merge_idxstats.nf'


workflow SPECIES_IDENTIFICATION {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads

    main:
    ch_versions         = Channel.empty()
    // define workflow name to be used in summary stats output
    def workflow_name = params.workflow_run_name ?: workflow.runName


    // Index the reference genome if it is not already indexed
    BOWTIE2_BUILD (reference, file(params.si_reference_database).getParent())
    ch_versions             = ch_versions.mix(BOWTIE2_BUILD.out.versions)


    // Map the reads to the reference genome
    SI_BOWTIE2 ( reads, BOWTIE2_BUILD.out.index_dir )
    ch_versions             = ch_versions.mix(SI_BOWTIE2.out.versions)
    ch_bam                  = SI_BOWTIE2.out.bam
    // Index the BAM file
    SI_BOWTIE2_INDEX ( ch_bam )
    ch_versions             = ch_versions.mix ( SI_BOWTIE2_INDEX.out.versions )


    // Filter the BAM file by mapping quality
    SI_SAMTOOLS_VIEW_MQ ( ch_bam )
    ch_versions         = ch_versions.mix(SI_SAMTOOLS_VIEW_MQ.out.versions)
    // Index the filtered BAM file
    SI_SAMTOOLS_VIEW_MQ_INDEX ( SI_SAMTOOLS_VIEW_MQ.out.bam )
    ch_versions         = ch_versions.mix(SI_SAMTOOLS_VIEW_MQ_INDEX.out.versions)


    // Prepare the BAM files for merging by updating the metadata
    ch_bam_lib_to_merge = SI_SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
        // update only the 'id' field in meta, keep all other fields
        [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1], 'library_type': meta.library_type, 'single_end': meta.single_end], bam]
    }.groupTuple()


    // Merge BAM files per library/PCR
    SI_SAMTOOLS_MERGE_LIB ( ch_bam_lib_to_merge, reference )
    ch_versions         = ch_versions.mix(SI_SAMTOOLS_MERGE_LIB.out.versions)
    // Index the merged BAM file per library/PCR
    SI_SAMTOOLS_MERGE_LIB_INDEX ( SI_SAMTOOLS_MERGE_LIB.out.bam )
    ch_versions         = ch_versions.mix(SI_SAMTOOLS_MERGE_LIB_INDEX.out.versions)


    // Remove duplicates from the merged BAM files per library/PCR
    SI_SAMREMOVEDUP_LIB ( SI_SAMTOOLS_MERGE_LIB.out.bam, reference )
    ch_versions         = ch_versions.mix(SI_SAMREMOVEDUP_LIB.out.versions)
    // Index the BAM file after removing duplicates
    SI_SAMREMOVEDUP_LIB_INDEX ( SI_SAMREMOVEDUP_LIB.out.dedup )
    ch_versions         = ch_versions.mix(SI_SAMREMOVEDUP_LIB_INDEX.out.versions)


    // Prepare the BAM files for merging by updating the metadata
    ch_bam_sample_to_merge = SI_SAMREMOVEDUP_LIB.out.dedup.map { meta, bam ->
            // update only the 'id' field in meta, keep all other fields
            [['id': meta.id.split("_")[0]], bam]
        }.groupTuple()


    // Merge BAM files per sample
    SI_SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, reference )
    ch_versions         = ch_versions.mix(SI_SAMTOOLS_MERGE_SAMPLE.out.versions)
    // Index the merged BAM file per sample
    SI_SAMTOOLS_MERGE_SAMPLE_INDEX ( SI_SAMTOOLS_MERGE_SAMPLE.out.bam )
    ch_versions         = ch_versions.mix(SI_SAMTOOLS_MERGE_SAMPLE_INDEX.out.versions)


    // Remove duplicates from the merged BAM files per sample
    SI_SAMREMOVEDUP_SAMPLE ( SI_SAMTOOLS_MERGE_SAMPLE.out.bam, reference )
    ch_versions         = ch_versions.mix(SI_SAMREMOVEDUP_SAMPLE.out.versions)
    // Index the BAM file after removing duplicates
    SI_SAMREMOVEDUP_SAMPLE_INDEX ( SI_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions         = ch_versions.mix(SI_SAMREMOVEDUP_SAMPLE_INDEX.out.versions)


    // Generate idxstats for the final BAM file
    SI_SAMTOOLS_IDXSTATS ( SI_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions         = ch_versions.mix(SI_SAMTOOLS_IDXSTATS.out.versions)

    // Merge idxstats from all samples
    ch_merged_input = SI_SAMTOOLS_IDXSTATS.out.idxstats
        .map { it -> [['id': "$workflow_name"], it[1]] }
        .groupTuple()

    SI_MERGE_IDXSTATS ( ch_merged_input )

    emit:
    reference_index         = BOWTIE2_BUILD.out.index_dir                 // channel: path(index)
    raw_bam                 = ch_bam                                      // channel: [ val(meta), [ bam ] ]
    raw_bam_bai             = SI_BOWTIE2_INDEX.out.bai                    // channel: [ val(meta), [ raw_bam_bai ] ]
    mq_filtered_bam         = SI_SAMTOOLS_VIEW_MQ.out.bam                 // channel: [ val(meta), [ mq_filtered_bam ] ]
    mq_filtered_bam_bai     = SI_SAMTOOLS_VIEW_MQ_INDEX.out.bai           // channel: [ val(meta), [ mq_filtered_bam_bai ] ]
    merged_bam_lib          = SI_SAMTOOLS_MERGE_LIB.out.bam               // channel: [ val(meta), [ merged_bam_lib ] ]
    merged_bam_lib_bai      = SI_SAMTOOLS_MERGE_LIB_INDEX.out.bai         // channel: [ val(meta), [ merged_bam_lib_bai ] ]
    dedup_bam_lib           = SI_SAMREMOVEDUP_LIB.out.dedup               // channel: [ val(meta), [ dedup_bam_lib ] ]
    dedup_bam_lib_bai       = SI_SAMREMOVEDUP_LIB_INDEX.out.bai           // channel: [ val(meta), [ dedup_bam_lib_bai ] ]
    merged_bam_sample       = SI_SAMTOOLS_MERGE_SAMPLE.out.bam            // channel: [ val(meta), [ merged_bam_sample ] ]
    merged_bam_sample_bai   = SI_SAMTOOLS_MERGE_SAMPLE_INDEX.out.bai      // channel: [ val(meta), [ merged_bam_sample_bai ] ]
    dedup_bam_sample        = SI_SAMREMOVEDUP_SAMPLE.out.dedup            // channel: [ val(meta), [ dedup_bam_sample ] ]
    dedup_bam_sample_bai    = SI_SAMREMOVEDUP_SAMPLE_INDEX.out.bai        // channel: [ val(meta), [ dedup_bam_sample_bai ] ]
    idxstats                = SI_SAMTOOLS_IDXSTATS.out.idxstats           // channel: [ val(meta), [ idxstats ] ]
    merged_idxstats         = SI_MERGE_IDXSTATS.out.merged_idxstats       // channel: [ val(meta), [ merged_idxstats ] ]
    versions                = ch_versions                                 // channel: [ versions.yml ]
}