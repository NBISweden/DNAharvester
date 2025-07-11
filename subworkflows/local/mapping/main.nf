#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools/samtools_faidx.nf'
include { BWA_INDEX      } from '../../../modules/local/bwa/index.nf'
include { BWA_ALN        } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE      } from '../../../modules/local/bwa/samse.nf'
include { BWA_ALN_MEM } from '../../../modules/local/bwa/bwa_aln_mem.nf'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'


workflow MAPPING {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads

    main:
    ch_versions         = Channel.empty()

    SAMTOOLS_FAIDX ( reference )
    ch_versions         = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Index the reference genome if it is not already indexed
    BWA_INDEX(reference, file(params.reference).getParent())
    ch_versions         = ch_versions.mix(BWA_INDEX.out.versions)
    ch_reference_index  = BWA_INDEX.out.index_dir

    // Map the reads to the reference genome using BWA
    if (params.mapping_tool == 'bwa-aln') {
        BWA_ALN ( reads, ch_reference_index )
        ch_versions         = ch_versions.mix(BWA_ALN.out.versions)
        ch_bwa_samse        = reads.join(BWA_ALN.out.sai)

        BWA_SAMSE ( ch_bwa_samse, ch_reference_index )
        ch_versions         = ch_versions.mix(BWA_SAMSE.out.versions)
        ch_bam              = BWA_SAMSE.out.bam

    } else if (params.mapping_tool == 'bwa-aln-mem') {
        BWA_ALN_MEM ( reads, ch_reference_index )
        ch_versions         = ch_versions.mix(BWA_ALN_MEM.out.versions)
        ch_bam              = BWA_ALN_MEM.out.bam
    }
    else {
        error "Invalid mapping tool specified: ${params.mapping_tool}. Use 'bwa-aln' or 'bwa-aln-mem'."
    }

    // Index the BAM file
    SAMTOOLS_INDEX ( ch_bam )
    ch_versions         = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    fai                 = SAMTOOLS_FAIDX.out.fai             // channel: path(index)
    index               = ch_reference_index                 // channel: path(index)
    bam                 = ch_bam                              // channel: [ val(meta), [ bam ] ]
    bai                 = SAMTOOLS_INDEX.out.bai             // channel: [ val(meta), [ bai ] ]
    versions            = ch_versions                        // channel: [ versions.yml ]
}