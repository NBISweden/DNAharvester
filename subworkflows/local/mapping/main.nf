#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools/faidx/main'
include { BWA_INDEX      } from '../../../modules/local/bwa/index.nf'
include { BWA_ALN        } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE      } from '../../../modules/local/bwa/samse.nf'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'


workflow MAPPING {
    take:
    reference
    bwa_index_reference
    reads // merged paired-end reads or trimmed single-end reads

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_FAIDX ( reference )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.fai)

    //BWA_INDEX ( reference )
    //ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    BWA_ALN ( reads, bwa_index_reference )
    ch_versions = ch_versions.mix(BWA_ALN.out.versions)

    ch_bwa_samse         = reads.join(BWA_ALN.out.sai)
    BWA_SAMSE ( ch_bwa_samse, bwa_index_reference )
    ch_versions = ch_versions.mix(BWA_SAMSE.out.versions)

    SAMTOOLS_INDEX ( BWA_SAMSE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)


    emit:
    fai            = SAMTOOLS_FAIDX.out.fai                      // channel: path(index)
    //index          = BWA_INDEX.out.index                         // channel: path(index)
    bam            = BWA_SAMSE.out.bam                           // channel: [ val(meta), [ bam ] ]
    bai            = SAMTOOLS_INDEX.out.bai                      // channel: [ val(meta), [ bai ] ]
    versions       = ch_versions                                 // channel: [ versions.yml ]
}