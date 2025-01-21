#! /usr/bin/env nextflow

include { ANGSD_DOHAPLOCALL  } from '../../../modules/local/angsd/dohaplocall/main'
include { ANGSD_HAPLOTOPLINK } from '../../../modules/local/angsd/haplotoplink/main'
include { PLINK_RECODE       } from '../../../modules/local/plink/recode/main'
include { HAPLOTOFASTA       } from '../../../modules/local/haplotofasta/main'

workflow RANDOM_SAMPLING_BAM {
    take:
    bam  // list of meta, bam, bai
    reference
    fai

    main:
    ch_versions                              = Channel.empty()

    ANGSD_DOHAPLOCALL ( bam )
    ch_versions                              = ch_versions.mix(ANGSD_DOHAPLOCALL.out.versions)

    ANGSD_HAPLOTOPLINK ( ANGSD_DOHAPLOCALL.out.haplo )
    ch_versions                              = ch_versions.mix(ANGSD_HAPLOTOPLINK.out.versions)

    PLINK_RECODE ( ANGSD_HAPLOTOPLINK.out.tfam, ANGSD_HAPLOTOPLINK.out.tped, reference )
    ch_versions                              = ch_versions.mix(PLINK_RECODE.out.versions)

    HAPLOTOFASTA ( ANGSD_DOHAPLOCALL.out.haplo, fai )
    ch_versions                              = ch_versions.mix(HAPLOTOFASTA.out.versions)

    emit:
    haplo                                    = ANGSD_DOHAPLOCALL.out.haplo                  // channel: [ val(meta), haplofile ]
    tfam                                     = ANGSD_HAPLOTOPLINK.out.tfam                  // channel: [ val(meta), tfam ]
    tped                                     = ANGSD_HAPLOTOPLINK.out.tped                  // channel: [ val(meta), tped ]
    vcf                                      = PLINK_RECODE.out.vcf                         // channel: [ val(meta), vcf ]
    fasta                                    = HAPLOTOFASTA.out.fasta                       // channel: [ val(meta), fasta ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}