#! /usr/bin/env nextflow

include { BWA_INDEX } from '../../../modules/local/bwa/index.nf'
include { BWA_ALN } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE } from '../../../modules/local/bwa/samse.nf'


workflow MAPPING {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads

    main:
    ch_versions = Channel.empty()

    BWA_INDEX ( reference )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    BWA_ALN ( reads, BWA_INDEX.out.index )
    ch_versions = ch_versions.mix(BWA_ALN.out.versions)
    BWA_SAMSE ( BWA_ALN.out.reads, BWA_ALN.out.sai, BWA_INDEX.out.index )
    ch_versions = ch_versions.mix(BWA_SAMSE.out.versions)

    emit:
    index          = BWA_INDEX.out.index                         // channel: path(index)
    sai            = BWA_ALN.out.sai                             // channel: [ val(meta), [ sai ] ]
    bam            = BWA_SAMSE.out.bam                           // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions                                 // channel: [ versions.yml ]
}