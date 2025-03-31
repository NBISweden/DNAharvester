#! /usr/bin/env nextflow

include { REPEATMODELER } from '../../../modules/local/repeat_cpg_masking/repeatmodeler'


workflow REPEAT_CPG_MASKING {
    take:
    reference

    main:
    ch_versions = Channel.empty()

    REPEATMODELER ( reference )
    ch_versions = ch_versions.mix(REPEATMODELER.out.versions)


    emit:
    consensi          = REPEATMODELER.out.consensi               // channel: path(consensi)
    families          = REPEATMODELER.out.families               // channel: path(families)
    versions          = ch_versions                              // channel: [ versions.yml ]
}