#! /usr/bin/env nextflow

include { FASTP } from '../../modules/local/fastp/main.nf'

workflow MERGE_FILTER_READS {
    take:
    reads

    main:
    FASTP ( reads )

    emit:
    reads_merged                  // channel: [ val(meta), [ reads ] ]
    reads_unmerged                // channel: [ val(meta), [ reads ] ]
    versions = FASTP.out.versions // channel: [ versions.yml ]
}