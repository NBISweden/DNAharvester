#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'
include { SAMREMOVEDUP   } from '../../../modules/local/samremovedup/main'


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

    SAMTOOLS_MERGE ( ch_bam_index_to_merge, reference, SAMTOOLS_FAIDX.out.fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    SAMREMOVEDUP ( SAMTOOLS_MERGE.out.bam )
    ch_versions = ch_versions.mix(SAMREMOVEDUP.out.versions)

    emit:
    fai              = SAMTOOLS_FAIDX.out.fai                         // channel: path(index)
    merged_bam       = SAMTOOLS_MERGE.out.bam                         // channel: [ val(meta), [ bam ] ]
    dedup            = SAMREMOVEDUP.out.dedup                         // channel: [ val(meta), [ bam ] ]
    versions         = ch_versions                                    // channel: [ versions.yml ]
}