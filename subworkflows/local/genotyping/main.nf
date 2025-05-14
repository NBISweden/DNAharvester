#! /usr/bin/env nextflow

include { ANGSD_GENOTYPING  } from '../../../modules/local/angsd/genotyping/main'

workflow GENOTYPING {
    take:
    bam  // list of meta, bam, bai
    reference
    fai

    main:
    ch_versions                              = Channel.empty()

    ANGSD_GENOTYPING ( bam, reference, fai )
    ch_versions                              = ch_versions.mix(ANGSD_GENOTYPING.out.versions)

    emit:
    angsd_log                                = ANGSD_GENOTYPING.out.angsd_log                   // channel: [ val(meta), log ]
    bamlist                                  = ANGSD_GENOTYPING.out.bamlist               // channel: [ val(meta), bamlist ]
    geno                                     = ANGSD_GENOTYPING.out.geno                   // channel: [ val(meta), geno ]
    mafs                                     = ANGSD_GENOTYPING.out.mafs                   // channel: [ val(meta), mafs ]
    beagle                                   = ANGSD_GENOTYPING.out.beagle                 // channel: [ val(meta), beagle ]
    bcf                                      = ANGSD_GENOTYPING.out.bcf                    // channel: [ val(meta), bcf ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}