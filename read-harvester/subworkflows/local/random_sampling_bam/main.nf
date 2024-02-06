#! /usr/bin/env nextflow

include { ANGSD_DOHAPLOCALL  } from '../../../modules/local/angsd/dohaplocall/main'
include { ANGSD_HAPLOTOPLINK } from '../../../modules/local/angsd/haplotoplink/main'

workflow RANDOM_SAMPLING_BAM {
    take:
    bam  // list of meta, bam, bai

    main:
    ch_versions                              = Channel.empty()

    // Remove *.bai from input channel
    ch_bam_for_angsd_dohaplocall = bam.map {
        meta, bam, bai -> [ meta, bam ] }

    ANGSD_DOHAPLOCALL ( ch_bam_for_angsd_dohaplocall )
    ch_versions                              = ch_versions.mix(ANGSD_DOHAPLOCALL.out.versions)

    ANGSD_HAPLOTOPLINK ( ANGSD_DOHAPLOCALL.out.haplo )
    ch_versions                              = ch_versions.mix(ANGSD_HAPLOTOPLINK.out.versions)

    emit:
    haplo                                    = ANGSD_DOHAPLOCALL.out.haplo                  // channel: [ val(meta), haplofile ]
    tfam                                     = ANGSD_HAPLOTOPLINK.out.tfam                  // channel: [ val(meta), tfam ]
    tped                                     = ANGSD_HAPLOTOPLINK.out.tped                  // channel: [ val(meta), tped ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}