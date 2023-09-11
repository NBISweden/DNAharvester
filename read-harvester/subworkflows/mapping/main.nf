#! /usr/bin/env nextflow

include { BWA_INDEX } from '../../modules/local/bwa/index.nf'
include { BWA_ALN } from '../../modules/local/bwa/aln.nf'
//include { BWA_SAMSE } from '../../modules/local/bwa/samse.nf'


workflow MAPPING {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads

    main:
    BWA_INDEX ( reference )
    BWA_ALN ( reads, BWA_INDEX.out.index )
    //BWA_SAMSE ( BWA_ALN.out.sai )

    emit:
    index          = BWA_INDEX.out.index                         // channel: path(index)
    sai            = BWA_ALN.out.sai                             // channel: [ val(meta), [ sai ] ]
    //bam            = BWA_SAMSE.out.bam                           // channel: [ val(meta), [ bam ] ]
    versions       = BWA_INDEX.out.versions                      // channel: [ versions.yml ]
}