#! /usr/bin/env nextflow

include { CONCATENATE_TARGET_DECOY_REFS                 } from '../../../modules/local/concatenate_target_decoy_refs/main'
include { BWA_INDEX as BWA_INDEX_CONCATENATED           } from '../../../modules/local/bwa/index.nf'
include { BWA_ALN as BWA_ALN_CONCATENATED               } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE as BWA_SAMSE_CONCATENATED           } from '../../../modules/local/bwa/samse.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CONCATENATED } from '../../../modules/nf-core/samtools/index/main'


workflow MAPPING_CONCATENATED_REFS {
    take:
    reference
    decoy
    reads // merged paired-end reads or trimmed single-end reads

    main:
    ch_versions = Channel.empty()

    // Concatenate the two references
    ch_concatenate_target_decoy_refs = reference.mix(decoy)
    CONCATENATE_TARGET_DECOY_REFS ( ch_concatenate_target_decoy_refs )

    // Index the concatented fasta file
    BWA_INDEX_CONCATENATED ( CONCATENATE_TARGET_DECOY_REFS.out.fasta )
    ch_versions = ch_versions.mix(BWA_INDEX_CONCATENATED.out.versions)

    // Map the reads to the concatenated fasta file
    BWA_ALN_CONCATENATED ( reads, BWA_INDEX_CONCATENATED.out.index )
    ch_versions = ch_versions.mix(BWA_ALN_CONCATENATED.out.versions)

    BWA_SAMSE_CONCATENATED ( BWA_ALN_CONCATENATED.out.reads, BWA_ALN_CONCATENATED.out.sai, BWA_INDEX_CONCATENATED.out.index )
    ch_versions = ch_versions.mix(BWA_SAMSE_CONCATENATED.out.versions)

    // Index the BAM file
    SAMTOOLS_INDEX_CONCATENATED ( BWA_SAMSE_CONCATENATED.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_CONCATENATED.out.versions)


    emit:
    index          = BWA_INDEX_CONCATENATED.out.index            // channel: path(index)
    sai            = BWA_ALN_CONCATENATED.out.sai                // channel: [ val(meta), [ sai ] ]
    bam            = BWA_SAMSE_CONCATENATED.out.bam              // channel: [ val(meta), [ bam ] ]
    bai            = SAMTOOLS_INDEX_CONCATENATED.out.bai         // channel: [ val(meta), [ bai ] ]
    versions       = ch_versions                                 // channel: [ versions.yml ]
}