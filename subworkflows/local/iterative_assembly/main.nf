#! /usr/bin/env nextflow

include { MAPPING_ITERATIVE_ASSEMBLER   as IA_MAPPING_ITERATIVE_ASSEMBLER   } from '../../../modules/local/mapping_iterative_assembler/mapping_iterative_assembler.nf'
include { CONSENSUS_CALL_MIA            as IA_CONSENSUS_CALL_MIA            } from '../../../modules/local/mapping_iterative_assembler/consensus_call_mia.nf'

workflow ITERATIVE_ASSEMBLY {
    take:
    reads // merged paired-end reads or trimmed single-end reads
    mt_reference // mitochondrial reference genome

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Merging all the fastq files per sample
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Merge all libraries/lanes per sample
    ch_reads_per_sample = reads.map { meta, fastq ->
        def new_meta = meta.clone()
            new_meta.id = meta.id.split("_")[0] // remove lane and library_id from id for merging all bams per sample_id
            new_meta.remove('single_end') // remove single_end info from meta as the same sample can have both single-end and paired-end data
            new_meta.remove('read_group') // remove read_group info from meta as the same sample can have multiple read groups
            new_meta.remove('library_type') // remove library_type info from meta as the same sample can have both double-stranded and single-stranded libraries
            [ new_meta, fastq ]
        }.groupTuple()

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Iterative Assembly using MIA and Consensus calling
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Run MIA - Mapping Iterative Assembler
    IA_MAPPING_ITERATIVE_ASSEMBLER ( ch_reads_per_sample, mt_reference )
    ch_versions         = ch_versions.mix( IA_MAPPING_ITERATIVE_ASSEMBLER.out.versions )

    // Consensus call
    IA_CONSENSUS_CALL_MIA( IA_MAPPING_ITERATIVE_ASSEMBLER.out.mia_maln_41 )
    ch_versions         = ch_versions.mix( IA_CONSENSUS_CALL_MIA.out.versions )
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    mia_maln_41         = IA_MAPPING_ITERATIVE_ASSEMBLER.out.mia_maln_41        // channel: [ val(meta), [ maln ] ]
    mia_fasta           = IA_CONSENSUS_CALL_MIA.out.mia_fasta                   // channel: [ val(meta), [ fasta ] ]
    versions            = ch_versions                                           // channel: [ versions.yml ]
}