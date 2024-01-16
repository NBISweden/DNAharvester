#! /usr/bin/env nextflow

include { SAMREMOVEDUP   } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'


workflow MERGE_DEDUP_BAMS {
    take:
    reference
    bams

    main:
    ch_versions = Channel.empty()

    SAMREMOVEDUP ( bams )
    ch_versions = ch_versions.mix(SAMREMOVEDUP.out.versions)

    SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    SAMTOOLS_MERGE ( SAMREMOVEDUP.out.dedup, reference, SAMTOOLS_FAIDX.out.fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    emit:
    fai              = SAMTOOLS_FAIDX.out.fai                         // channel: path(index)
    dedup            = SAMREMOVEDUP.out.dedup                         // channel: [ val(meta), [ sai ] ]
    merged_bam       = SAMTOOLS_MERGE.out.bam                         // channel: [ val(meta), [ bam ] ]
    versions         = ch_versions                                    // channel: [ versions.yml ]
}