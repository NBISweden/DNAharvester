#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX                          } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_INDEX  } from '../../../modules/nf-core/samtools/merge/main'
include { SAMREMOVEDUP as SAMREMOVEDUP_INDEX      } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_SAMPLE } from '../../../modules/nf-core/samtools/merge/main'
include { SAMREMOVEDUP as SAMREMOVEDUP_SAMPLE     } from '../../../modules/local/samremovedup/main'

workflow MERGE_DEDUP_BAMS {
    take:
    reference
    bam

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_bam_index_to_merge = bam.map {
        meta, bam -> [ ['id':meta.id.split("_")[0] + "_" + meta.id.split("_")[1]], bam ]
        }
        .groupTuple()
    ch_bam_index_to_merge.view()

    SAMTOOLS_MERGE_INDEX ( ch_bam_index_to_merge, reference, SAMTOOLS_FAIDX.out.fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_INDEX.out.versions)

    SAMREMOVEDUP_INDEX ( SAMTOOLS_MERGE_INDEX.out.bam )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_INDEX.out.versions)

    ch_bam_sample_to_merge = SAMTOOLS_MERGE_INDEX.out.bam.map {
        meta, bam -> [ ['id':meta.id.split("_")[0]], bam ]
        }
        .groupTuple()
    ch_bam_sample_to_merge.view()

    SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, reference, SAMTOOLS_FAIDX.out.fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE_SAMPLE.out.versions)

    SAMREMOVEDUP_SAMPLE ( SAMTOOLS_MERGE_SAMPLE.out.bam )
    ch_versions = ch_versions.mix(SAMREMOVEDUP_SAMPLE.out.versions)

    emit:
    fai               = SAMTOOLS_FAIDX.out.fai                         // channel: path(index)
    merged_bam_index  = SAMTOOLS_MERGE_INDEX.out.bam                   // channel: [ val(meta), [ bam ] ]
    dedup_index       = SAMREMOVEDUP_INDEX.out.dedup                   // channel: [ val(meta), [ bam ] ]
    merged_bam_sample = SAMTOOLS_MERGE_SAMPLE.out.bam                  // channel: [ val(meta), [ bam ] ]
    dedup_sample      = SAMREMOVEDUP_SAMPLE.out.dedup                  // channel: [ val(meta), [ bam ] ]
    versions          = ch_versions                                    // channel: [ versions.yml ]
}