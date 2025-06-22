#! /usr/bin/env nextflow

include { MAPPING_ITERATIVE_ASSEMBLER } from '../../../modules/local/mapping_iterative_assembler/main.nf'
include { CONSENSUS_CALL_MIA } from '../../../modules/local/mapping_iterative_assembler/consensus_call_mia.nf'

workflow ITERATIVE_ASSEMBLY {
    take:
    reads // merged paired-end reads or trimmed single-end reads
    mt_reference // mitochondrial reference genome

    main:
    ch_versions = Channel.empty()

    // Run MIA - Mapping Iterative Assembler
    MAPPING_ITERATIVE_ASSEMBLER (reads, mt_reference)
    ch_versions = ch_versions.mix(MAPPING_ITERATIVE_ASSEMBLER.out.versions)

    // Consensus call
    CONSENSUS_CALL_MIA(MAPPING_ITERATIVE_ASSEMBLER.out.mia_maln_41)
    ch_versions = ch_versions.mix(CONSENSUS_CALL_MIA.out.versions)


    emit:
    mia_maln_41       = MAPPING_ITERATIVE_ASSEMBLER.out.mia_maln_41     // channel: [ val(meta), [ maln ] ]
    mia_fasta         = CONSENSUS_CALL_MIA.out.mia_fasta                // channel: [ val(meta), [ fasta ] ]
    versions          = ch_versions                                     // channel: [ versions.yml ]
}