#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX } from '../../../modules/local/samtools/samtools_faidx.nf'
include { BWA_INDEX      } from '../../../modules/local/bwa/index.nf'
include { BOWTIE2_BUILD  } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BWA_ALN               } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE      } from '../../../modules/local/bwa/samse.nf'
include { BWA_ALN_MEM    } from '../../../modules/local/bwa/bwa_aln_mem.nf'
include { BOWTIE2        } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

include { BWA_ALN as BWA_ALN_R1 } from '../../../modules/local/bwa/aln.nf'
include { BWA_ALN as BWA_ALN_R2 } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMPE      } from '../../../modules/local/bwa/sampe.nf'
include { SAMTOOLS_MERGE } from '../../../modules/local/samtools/samtools_merge.nf'


workflow MAPPING {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads
    unmerged_reads // unmerged paired-end reads

    main:
    ch_versions         = Channel.empty()

    SAMTOOLS_FAIDX ( reference )
    ch_versions         = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)


    ////////////////////////////////////////////////////////////////////////////
    // Index the reference genome if it is not already indexed
    ////////////////////////////////////////////////////////////////////////////

    if (params.mapping_tool == 'bwa-aln' || params.mapping_tool == 'bwa-aln-mem') {
        // Build the BWA index
        BWA_INDEX (reference, file(params.reference).getParent())
        ch_versions         = ch_versions.mix(BWA_INDEX.out.versions)
        ch_reference_index  = BWA_INDEX.out.index_dir
    } else if (params.mapping_tool == 'bowtie2') {
        // Build the Bowtie2 index
        BOWTIE2_BUILD (reference, file(params.reference).getParent())
        ch_versions         = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        ch_reference_index  = BOWTIE2_BUILD.out.index_dir
    } else {
        error "Invalid mapping tool specified: ${params.mapping_tool}. Use 'bwa-aln', 'bwa-aln-mem', or 'bowtie2'."
    }

    ////////////////////////////////////////////////////////////////////////////
    // Mapping merged reads or single-end reads
    ////////////////////////////////////////////////////////////////////////////

    if (params.mapping_tool == 'bwa-aln') {
        BWA_ALN ( reads, ch_reference_index)
        ch_versions         = ch_versions.mix(BWA_ALN.out.versions)
        ch_bwa_samse        = reads.join(BWA_ALN.out.sai)

        BWA_SAMSE ( ch_bwa_samse, ch_reference_index )
        ch_versions         = ch_versions.mix(BWA_SAMSE.out.versions)
        ch_bam              = BWA_SAMSE.out.bam

    } else if (params.mapping_tool == 'bwa-aln-mem') {
        BWA_ALN_MEM ( reads, ch_reference_index )
        ch_versions         = ch_versions.mix(BWA_ALN_MEM.out.versions)
        ch_bam              = BWA_ALN_MEM.out.bam
    } else if (params.mapping_tool == 'bowtie2') {
        BOWTIE2 ( reads, ch_reference_index )
        ch_versions         = ch_versions.mix(BOWTIE2.out.versions)
        ch_bam              = BOWTIE2.out.bam
    } else {
        error "Invalid mapping tool specified: ${params.mapping_tool}. Use 'bwa-aln' or 'bwa-aln-mem' or 'bowtie2'."
    }

    ////////////////////////////////////////////////////////////////////////////
    // Mapping unmerged reads if params.keep_unmerged is true
    ////////////////////////////////////////////////////////////////////////////

    // processed unmerged reads if provided
    if (params.keep_unmerged.toBoolean()) {

        ch_paired_end_reads = ch_bam.filter { meta, file -> !meta.single_end }
        ch_single_end_reads = ch_bam.filter { meta, file -> meta.single_end }

        if (ch_paired_end_reads) {
            // Prepare channels for unmerged reads. Add R1 and R2 labels to the meta.id. This is needed to avoid file name clashes in the BWA_ALN process output
            ch_R1 = unmerged_reads.map { meta, file1, file2 ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}-R1"
                [new_meta, file1]
            }
            ch_R2 = unmerged_reads.map { meta, file1, file2 ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}-R2"
                [new_meta, file2]
            }

            // Align unmerged reads
            if (params.mapping_tool == 'bwa-aln') {
                BWA_ALN_R1 ( ch_R1, ch_reference_index)
                ch_versions             = ch_versions.mix(BWA_ALN_R1.out.versions)
                BWA_ALN_R2 ( ch_R2, ch_reference_index)
                ch_versions             = ch_versions.mix(BWA_ALN_R2.out.versions)
                // Fix meta.id to original id without -R1 or -R2 suffix after alignment
                ch_sai_R1 = BWA_ALN_R1.out.sai.map { meta, file ->
                    meta.id = meta.id.replaceAll(/-R1$/, '')
                    [meta, file]
                }
                ch_sai_R2 = BWA_ALN_R2.out.sai.map { meta, file ->
                    meta.id = meta.id.replaceAll(/-R2$/, '')
                    [meta, file]
                }
                // Join the aligned unmerged reads with their respective SAI files
                ch_bwa_sampe_unmerged = unmerged_reads.join(ch_sai_R1).join(ch_sai_R2)
                // Run BWA SAMPE
                BWA_SAMPE ( ch_bwa_sampe_unmerged, ch_reference_index )
                ch_versions             = ch_versions.mix(BWA_SAMPE.out.versions)
                ch_bam_unmerged         = BWA_SAMPE.out.bam

                // merge merged and unmerged bam files
                ch_bam_merged_unmerged = ch_bam.join(ch_bam_unmerged)
                        .map { meta, file1, file2 -> [meta, [file1, file2]] }

                SAMTOOLS_MERGE ( ch_bam_merged_unmerged, reference )
                ch_bam = SAMTOOLS_MERGE.out.bam.mix(ch_single_end_reads)
                ch_versions             = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    // Index the BAM file
    SAMTOOLS_INDEX ( ch_bam )
    ch_versions             = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    fai                     = SAMTOOLS_FAIDX.out.fai             // channel: path(index)
    index                   = ch_reference_index                 // channel: path(index)
    bam                     = ch_bam                             // channel: [ val(meta), [ bam ] ]
    bai                     = SAMTOOLS_INDEX.out.bai             // channel: [ val(meta), [ bai ] ]
    versions                = ch_versions                        // channel: [ versions.yml ]
}