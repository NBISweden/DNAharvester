#! /usr/bin/env nextflow

include { FASTP } from '../../../modules/local/fastp/pairedend.nf'

workflow FASTQ_PROCESSING {
    take:
    reads

    main:
    FASTP ( reads )

    emit:
    reads          = FASTP.out.reads                         // channel: [ val(meta), [ reads ] ]. Merged paired-end reads.
    json           = FASTP.out.json                          // channel: [ val(meta), [ reads ] ]
    fastp_log      = FASTP.out.log                           // channel: [ val(meta), [ reads ] ]
    reads_unmerged = FASTP.out.reads_unmerged                // channel: [ val(meta), [ reads ] ]
    versions       = FASTP.out.versions                      // channel: [ versions.yml ]
}