#! /usr/bin/env nextflow

include { SAMREMOVEDUP   } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_FAIDX } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_MERGE } from '../../../modules/nf-core/samtools/merge/main'


workflow MERGE_DEDUP_BAMS {
    take:
    reference
    bam

    main:
    ch_versions = Channel.empty()

    SAMREMOVEDUP ( bam )
    ch_versions = ch_versions.mix(SAMREMOVEDUP.out.versions)

    SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_bam_index_to_merge = SAMREMOVEDUP.out.dedup.map {
        meta, bam -> [ ['id':meta.id.split("_")[0] + "_" + meta.id.split("_")[1]], bam ]
        }
        .groupTuple()
    ch_bam_index_to_merge.view()

    SAMTOOLS_MERGE ( ch_bam_index_to_merge, reference, SAMTOOLS_FAIDX.out.fai )
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

    emit:
    fai              = SAMTOOLS_FAIDX.out.fai                         // channel: path(index)
    dedup            = SAMREMOVEDUP.out.dedup                         // channel: [ val(meta), [ sai ] ]
    merged_bam       = SAMTOOLS_MERGE.out.bam                         // channel: [ val(meta), [ bam ] ]
    versions         = ch_versions                                    // channel: [ versions.yml ]
}