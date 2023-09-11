#! /usr/bin/env nextflow

include { FASTP } from '../../../modules/local/fastp/pairedend.nf'

workflow MERGE_FILTER_READS {
    take:
    reads

    main:
    FASTP ( reads )

    emit:
    reads          = FASTP.out.reads                         // channel: [ val(meta), [ reads ] ]. Merged paired-end reads.
    reads_unmerged = FASTP.out.reads_unmerged                // channel: [ val(meta), [ reads ] ]
    versions       = FASTP.out.versions                      // channel: [ versions.yml ]
}