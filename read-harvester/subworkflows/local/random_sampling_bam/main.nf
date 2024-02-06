#! /usr/bin/env nextflow

include { ANGSD_DOHAPLOCALL } from '../../../modules/local/angsd/dohaplocall/main'

workflow RANDOM_SAMPLING_BAM {
    take:
    bam  // list of meta, bam, bai

    main:

    // Remove *.bai from input channel
    ch_bam_for_angsd_dohaplocall = bam.map {
        meta, bam, bai -> [ meta, bam ] }

    ANGSD_DOHAPLOCALL ( ch_bam_for_angsd_dohaplocall )

    emit:
    haplo          = ANGSD_DOHAPLOCALL.out.haplo                         // channel: [ val(meta), haplofile ]
    versions       = ANGSD_DOHAPLOCALL.out.versions                      // channel: [ versions.yml ]
}