#! /usr/bin/env nextflow

include { SEQ_STATS as SEQ_STATS_LIB        } from '../../../modules/local/stats_output/seq_stats.nf'
include { SORT_STATS as SORT_STATS_LIB      } from '../../../modules/local/stats_output/sort_stats.nf'
include { SEQ_STATS as SEQ_STATS_SAMPLE     } from '../../../modules/local/stats_output/seq_stats.nf'
include { SORT_STATS as SORT_STATS_SAMPLE   } from '../../../modules/local/stats_output/sort_stats.nf'

workflow STATS_OUTPUT {
    take:
    reads
    fastp_log
    raw_bam_flagstat
    mq_filtered_bam_flagstat
    dedup_lib_flagstat
    dedup_lib
    dedup_sample_flagstat
    dedup_sample


    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////
    // Generating seq stats per library
    ////////////////////////////////////////

    // Processing channel for merging
    ch_reads = reads.map { meta, data ->
        [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1], 'library_type': meta.library_type, 'single_end': meta.single_end], data[0]]
    }.groupTuple()

    ch_fastp_log = fastp_log.map { meta, data ->
        [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1], 'library_type': meta.library_type, 'single_end': meta.single_end], data ]
    }.groupTuple()

    ch_raw_bam_flagstat = raw_bam_flagstat.map { meta, data ->
        [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1], 'library_type': meta.library_type, 'single_end': meta.single_end], data ]
    }.groupTuple()

    ch_mq_filtered_bam_flagstat = mq_filtered_bam_flagstat.map { meta, data ->
        [['id': meta.id.split("_")[0] + "_" + meta.id.split("_")[1], 'library_type': meta.library_type, 'single_end': meta.single_end], data ]
    }.groupTuple()

    // merging different channels
    ch_seq_stats = ch_reads
        .join(ch_fastp_log)
        .join(ch_raw_bam_flagstat)
        .join(ch_mq_filtered_bam_flagstat)
        .join(dedup_lib_flagstat)
        .join(dedup_lib)

    // run the SEQ_STATS process
    SEQ_STATS_LIB ( ch_seq_stats )
    ch_versions = ch_versions.mix(SEQ_STATS_LIB.out.versions)
    // Concatenate all output files
    def seq_stats_lib = SEQ_STATS_LIB.out.stats_txt
        .map { it[1] }
        .collectFile(name: "${params.workflowname}_lib_stats", keepHeader: true, skip: 1, sort: true)

    // sort the stats output file
    SORT_STATS_LIB ( seq_stats_lib )

    /////////////////////////////////////////
    // Generating seq stats per sample
    /////////////////////////////////////////

    // Prepare channels for sample statistics
    ch_reads_sample = reads.map { meta, data ->
        [['id': meta.id.split("_")[0]], data[0]]
    }.groupTuple()

    ch_fastp_log_sample = fastp_log.map { meta, data ->
        [['id': meta.id.split("_")[0]], data ]
    }.groupTuple()

    ch_raw_bam_flagstat_sample = raw_bam_flagstat.map { meta, data ->
        [['id': meta.id.split("_")[0]], data ]
    }.groupTuple()

    ch_mq_filtered_bam_flagstat_sample = mq_filtered_bam_flagstat.map { meta, data ->
        [['id': meta.id.split("_")[0]], data ]
    }.groupTuple()

    // merging different channels
    ch_seq_stats_sample = ch_reads_sample
        .join(ch_fastp_log_sample)
        .join(ch_raw_bam_flagstat_sample)
        .join(ch_mq_filtered_bam_flagstat_sample)
        .join(dedup_sample_flagstat)
        .join(dedup_sample)


    // run the SEQ_STATS process
    SEQ_STATS_SAMPLE ( ch_seq_stats_sample )
    ch_versions = ch_versions.mix(SEQ_STATS_SAMPLE.out.versions)
    // Concatenate all output files
    def seq_stats_sample = SEQ_STATS_SAMPLE.out.stats_txt
        .map { it[1] }
        .collectFile(name: "${params.workflowname}_sample_stats", keepHeader: true, skip: 1, sort: true)

    // sort the stats output file
    SORT_STATS_SAMPLE ( seq_stats_sample )

    // Emit channels
    emit:
    seq_stats_lib       = SORT_STATS_LIB.out.sorted_file        // channel: [ ${params.workflowname}_lib_stats.sorted.txt ]
    seq_stats_sample    = SORT_STATS_SAMPLE.out.sorted_file     // channel: [ ${params.workflowname}_sample_stats.sorted.txt ]
    versions            = ch_versions                           // channel: [ versions.yml ]
}