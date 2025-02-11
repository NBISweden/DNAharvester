#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools/faidx/main'
include { BWA_INDEX      } from '../../../modules/local/bwa/index.nf'
include { BWA_ALN        } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE      } from '../../../modules/local/bwa/samse.nf'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'


workflow MAPPING {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // This will run only if the index is not already present
    BWA_INDEX ( reference )

    ch_reference_index = BWA_INDEX.out.index
        .map { id, files ->
            // Extracting the parent directory from the first file in the list
            def parentDir = files[0].getParent()
            return [id, parentDir] // Return necessary data
        }
        .collect()

    BWA_ALN ( reads, ch_reference_index )
    ch_versions = ch_versions.mix(BWA_ALN.out.versions)

    ch_bwa_samse         = reads.join(BWA_ALN.out.sai)
    BWA_SAMSE ( ch_bwa_samse, ch_reference_index )
    ch_versions = ch_versions.mix(BWA_SAMSE.out.versions)

    SAMTOOLS_INDEX ( BWA_SAMSE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    fai            = SAMTOOLS_FAIDX.out.fai                      // channel: path(index)
    index          = ch_reference_index                          // channel: path(index)
    bam            = BWA_SAMSE.out.bam                           // channel: [ val(meta), [ bam ] ]
    bai            = SAMTOOLS_INDEX.out.bai                      // channel: [ val(meta), [ bai ] ]
    versions       = ch_versions                                 // channel: [ versions.yml ]
}