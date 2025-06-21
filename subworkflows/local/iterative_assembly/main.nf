#! /usr/bin/env nextflow

include { MAPPING_ITERATIVE_ASSEMBLER } from '../../../modules/local/mapping_iterative_assembler/main.nf'

workflow ITERATIVE_ASSEMBLY {
    take:
    reads // merged paired-end reads or trimmed single-end reads
    mt_reference // mitochondrial reference genome

    main:
    ch_versions = Channel.empty()

    // Run MIA - Mapping Iterative Assembler


    mt_reference.view()


    MAPPING_ITERATIVE_ASSEMBLER (reads, mt_reference)
    ch_versions = ch_versions.mix(MAPPING_ITERATIVE_ASSEMBLER.out.versions)


    emit:
    mia_maln       = MAPPING_ITERATIVE_ASSEMBLER.out.mia_maln    // channel: [ val(meta), [ maln ] ]
    versions       = ch_versions                                 // channel: [ versions.yml ]
}