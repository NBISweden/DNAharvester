#! /usr/bin/env nextflow

include { GATK_REALIGNERTARGETCREATOR } from '../modules/nf-core/gatk/realignertargetcreator/main'  
include { GATK_INDELREALIGNER         } from '../modules/nf-core/gatk/indelrealigner/main'

workflow REALIGN_INDELS {
    take:
    bam, bai
    reference
    fai
    dict

    main:
    ch_versions = Channel.empty()

    GATK_REALIGNERTARGETCREATOR ( bam, bai, reference, fai, dict )
    ch_versions = ch_versions.mix(GATK_REALIGNERTARGETCREATOR.out.versions)

    GATK_INDELREALIGNER ( bam, bai, GATK_REALIGNERTARGETCREATOR.out.intervals, reference, fai, dict )
    ch_versions = ch_versions.mix(GATK_INDELREALIGNER.out.versions)

    emit:
    realigned         = GATK_INDELREALIGNER.out.bam                    // channel: [ val(meta), [ bam ] ]
    versions          = ch_versions                                    // channel: [ versions.yml ]
}