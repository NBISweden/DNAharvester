#! /usr/bin/env nextflow

// BWA mapping modules
include { BWA_INDEX as TC_BWA_INDEX                         } from '../../../modules/local/bwa/bwa_index.nf'
include { BWA_ALN as TC_BWA_ALN                             } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_SAMSE as TC_BWA_SAMSE                         } from '../../../modules/local/bwa/bwa_samse.nf'
include { BWA_SAMPE as TC_BWA_SAMPE                         } from '../../../modules/local/bwa/bwa_sampe.nf'
include { BWA_MEM as TC_BWA_MEM                             } from '../../../modules/local/bwa/bwa_mem.nf'

// Bowtie2 mapping modules
include { BOWTIE2_BUILD as TC_BOWTIE2_BUILD                 } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BOWTIE2 as TC_BOWTIE2                             } from '../../../modules/local/bowtie2/bowtie2.nf'

// SAMTOOLS indexing module for raw BAM files
include { SAMTOOLS_INDEX as TC_RAW_BAM_INDEX                } from '../../../modules/nf-core/samtools/index/main'

// Mapping quality filter
include { SAMTOOLS_VIEW_MQ as TC_SAMTOOLS_VIEW_MQ           } from '../../../modules/local/samtools/samtools_view_mq.nf'

// Merge BAM files per library/PCR
include { SAMTOOLS_MERGE as TC_SAMTOOLS_MERGE_LIB           } from '../../../modules/local/samtools/samtools_merge.nf'

// Remove duplicates from BAM files merged per library/PCR
include { SAMREMOVEDUP as TC_SAMREMOVEDUP_LIB               } from '../../../modules/local/samremovedup/main'

// Merge BAM files per sample
include { SAMTOOLS_MERGE as TC_SAMTOOLS_MERGE_SAMPLE        } from '../../../modules/local/samtools/samtools_merge.nf'

// Remove duplicates from BAM files merged per sample
include { SAMREMOVEDUP as TC_SAMREMOVEDUP_SAMPLE            } from '../../../modules/local/samremovedup/main'
include { SAMTOOLS_INDEX as TC_SAMREMOVEDUP_SAMPLE_INDEX    } from '../../../modules/nf-core/samtools/index/main'

// SAMTOOLS_IDXSTATS
include { SAMTOOLS_IDXSTATS as TC_SAMTOOLS_IDXSTATS         } from '../../../modules/local/samtools/samtools_idxstats.nf'
include { MERGE_IDXSTATS as TC_MERGE_IDXSTATS               } from '../../../modules/local/merge_idxstats/merge_idxstats.nf'

workflow TAXONOMIC_CLASSIFICATION {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads
    workflow_name

    main:
    ch_versions         = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////
    // Index the reference database if it is not already indexed
    ////////////////////////////////////////////////////////////////////////////

    // Build the BWA index - only if BWA is selected as mapping tool
    if (params.tc_mapping_tool == 'bwa-aln' || params.tc_mapping_tool == 'bwa-mem') {
        TC_BWA_INDEX (reference, file(params.tc_reference_database).getParent())
        ch_versions         = ch_versions.mix(TC_BWA_INDEX.out.versions)
        ch_bwa_index        = TC_BWA_INDEX.out.index_dir
    }

    // Build the Bowtie2 index - only if Bowtie2 is selected as mapping tool
    if (params.tc_mapping_tool == 'bowtie2') {
        TC_BOWTIE2_BUILD (reference, file(params.tc_reference_database).getParent())
        ch_versions         = ch_versions.mix(TC_BOWTIE2_BUILD.out.versions)
        ch_bowtie2_index    = TC_BOWTIE2_BUILD.out.index_dir
    }

    ////////////////////////////////////////////////////////////////////////////
    // Mapping the reads to the reference database
    ////////////////////////////////////////////////////////////////////////////

    // BWA ALN
    if (params.tc_mapping_tool == 'bwa-aln' ) {
        TC_BWA_ALN ( reads, ch_bwa_index )
        ch_versions         = ch_versions.mix(TC_BWA_ALN.out.versions)
        ch_bwa_samse        = reads.join(TC_BWA_ALN.out.sai)
        TC_BWA_SAMSE ( ch_bwa_samse, ch_bwa_index )
        ch_versions         = ch_versions.mix(TC_BWA_SAMSE.out.versions)
        ch_bam              = TC_BWA_SAMSE.out.bam
    }
    // BWA MEM
    if (params.tc_mapping_tool == 'bwa-mem') {
        TC_BWA_MEM ( reads, ch_bwa_index )
        ch_versions         = ch_versions.mix(TC_BWA_MEM.out.versions)
        ch_bam              = TC_BWA_MEM.out.bam
    }
    // BOWTIE2
    if (params.tc_mapping_tool == 'bowtie2') {
        TC_BOWTIE2 ( reads, ch_bowtie2_index )
        ch_versions         = ch_versions.mix(TC_BOWTIE2.out.versions)
        ch_bam              = TC_BOWTIE2.out.bam
    }

    ////////////////////////////////////////////////////////////////////////////
    // BAM processing
    ////////////////////////////////////////////////////////////////////////////

    // Index the BAM file
    TC_RAW_BAM_INDEX ( ch_bam )
    ch_versions             = ch_versions.mix ( TC_RAW_BAM_INDEX.out.versions )

    // Filter the BAM file by mapping quality
    TC_SAMTOOLS_VIEW_MQ ( ch_bam )
    ch_versions         = ch_versions.mix(TC_SAMTOOLS_VIEW_MQ.out.versions)

    // Prepare the BAM files for merging by updating the metadata
    ch_bam_lib_to_merge = TC_SAMTOOLS_VIEW_MQ.out.bam.map { meta, bam ->
        // Update ID to "<sample>_<library>"
        def new_id = meta.id.split('_')
        meta.id = "${new_id[0]}_${new_id[1]}"
        [meta, bam]
    }
    .groupTuple()

    // Merge BAM files per library/PCR
    TC_SAMTOOLS_MERGE_LIB ( ch_bam_lib_to_merge, reference )
    ch_versions         = ch_versions.mix(TC_SAMTOOLS_MERGE_LIB.out.versions)

    // Remove duplicates from the merged BAM files per library/PCR
    TC_SAMREMOVEDUP_LIB ( TC_SAMTOOLS_MERGE_LIB.out.bam, reference )
    ch_versions         = ch_versions.mix(TC_SAMREMOVEDUP_LIB.out.versions)

    // Prepare the BAM files for merging by updating the metadata
    ch_bam_sample_to_merge = TC_SAMREMOVEDUP_LIB.out.dedup.map { meta, bam ->
        def new_id = meta.id.split('_')
        // Update ID to "<sample>"
        meta.id = "${new_id[0]}"
        [meta, bam]
    }.groupTuple()


    // Merge BAM files per sample
    TC_SAMTOOLS_MERGE_SAMPLE ( ch_bam_sample_to_merge, reference )
    ch_versions         = ch_versions.mix(TC_SAMTOOLS_MERGE_SAMPLE.out.versions)

    // Remove duplicates from the merged BAM files per sample
    TC_SAMREMOVEDUP_SAMPLE ( TC_SAMTOOLS_MERGE_SAMPLE.out.bam, reference )
    ch_versions         = ch_versions.mix(TC_SAMREMOVEDUP_SAMPLE.out.versions)
    // Index the BAM file after removing duplicates
    TC_SAMREMOVEDUP_SAMPLE_INDEX ( TC_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions         = ch_versions.mix(TC_SAMREMOVEDUP_SAMPLE_INDEX.out.versions)

    ////////////////////////////////////////////////////////////////////////////
    // Generate idxstats
    ////////////////////////////////////////////////////////////////////////////

    // Generate idxstats for the final BAM file
    TC_SAMTOOLS_IDXSTATS ( TC_SAMREMOVEDUP_SAMPLE.out.dedup )
    ch_versions         = ch_versions.mix(TC_SAMTOOLS_IDXSTATS.out.versions)

    // Merge idxstats from all samples
    ch_merged_input = TC_SAMTOOLS_IDXSTATS.out.idxstats
        .map { it -> [['id': "$workflow_name"], it[1]] }
        .groupTuple()

    // Merge idxstats
    TC_MERGE_IDXSTATS ( ch_merged_input )
    ch_versions         = ch_versions.mix(TC_MERGE_IDXSTATS.out.versions)

    ///////////////////////////////////////////////////////////////////////////

    emit:
    raw_bam                 = ch_bam                                      // channel: [ val(meta), [ bam ] ]
    raw_bam_bai             = TC_RAW_BAM_INDEX.out.bai                    // channel: [ val(meta), [ raw_bam_bai ] ]
    mq_filtered_bam         = TC_SAMTOOLS_VIEW_MQ.out.bam                 // channel: [ val(meta), [ mq_filtered_bam ] ]
    merged_bam_lib          = TC_SAMTOOLS_MERGE_LIB.out.bam               // channel: [ val(meta), [ merged_bam_lib ] ]
    dedup_bam_lib           = TC_SAMREMOVEDUP_LIB.out.dedup               // channel: [ val(meta), [ dedup_bam_lib ] ]
    merged_bam_sample       = TC_SAMTOOLS_MERGE_SAMPLE.out.bam            // channel: [ val(meta), [ merged_bam_sample ] ]
    dedup_bam_sample        = TC_SAMREMOVEDUP_SAMPLE.out.dedup            // channel: [ val(meta), [ dedup_bam_sample ] ]
    dedup_bam_sample_bai    = TC_SAMREMOVEDUP_SAMPLE_INDEX.out.bai        // channel: [ val(meta), [ dedup_bam_sample_bai ] ]
    idxstats                = TC_SAMTOOLS_IDXSTATS.out.idxstats           // channel: [ val(meta), [ idxstats ] ]
    merged_idxstats         = TC_MERGE_IDXSTATS.out.merged_idxstats       // channel: [ val(meta), [ merged_idxstats ] ]
    versions                = ch_versions                                 // channel: [ versions.yml ]
}