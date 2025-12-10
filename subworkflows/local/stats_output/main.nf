#! /usr/bin/env nextflow

include { SEQ_STATS_LIB                     } from '../../../modules/local/stats_output/seq_stats_lib.nf'
include { SORT_STATS as SORT_STATS_LIB      } from '../../../modules/local/stats_output/sort_stats.nf'
include { SEQ_STATS_SAMPLE                  } from '../../../modules/local/stats_output/seq_stats_sample.nf'
include { SORT_STATS as SORT_STATS_SAMPLE   } from '../../../modules/local/stats_output/sort_stats.nf'

workflow STATS_OUTPUT {
    take:
    workflow_name
    reads
    fastp_log
    raw_bam_flagstat
    mq_filtered_bam_flagstat
    dedup_lib_flagstat
    dedup_lib
    dedup_lib_index
    dedup_sample_flagstat
    dedup_sample
    dedup_sample_index
    decoy_flagstat
    dpstats

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Generating seq stats per library
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Processing channel for merging
    ch_reads = reads.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0] + "_" + meta.id.split("_")[1]]
        new_meta.remove('single_end')
        [ new_meta, data[0] ]
    }.groupTuple()

    ch_fastp_log = fastp_log.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0] + "_" + meta.id.split("_")[1]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_raw_bam_flagstat = raw_bam_flagstat.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0] + "_" + meta.id.split("_")[1]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_mq_filtered_bam_flagstat = mq_filtered_bam_flagstat.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0] + "_" + meta.id.split("_")[1]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_dedup_lib_flagstat = dedup_lib_flagstat.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0] + "_" + meta.id.split("_")[1]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_dedup_lib = dedup_lib.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0] + "_" + meta.id.split("_")[1]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    // merging different channels
    ch_seq_stats = ch_reads
        .join(ch_fastp_log)
        .join(ch_raw_bam_flagstat)
        .join(ch_mq_filtered_bam_flagstat)
        .join(ch_dedup_lib_flagstat)
        .join(ch_dedup_lib)

    // If competitive reference is used, include decoy flagstat
    if (params.competitive_reference) {
        ch_decoy_flagstat = decoy_flagstat.map { meta, data ->
            def new_meta = meta + [id: meta.id.split("_")[0] + "_" + meta.id.split("_")[1]]
            new_meta.remove('single_end')
            [ new_meta, data ]
        }.groupTuple()

        ch_seq_stats = ch_seq_stats.join(ch_decoy_flagstat)
    } else {
        // If no decoy flagstat is provided, append a dummy path directly
        ch_seq_stats = ch_seq_stats.map {
            it + [ [file('/dev/null')] ]
        }
    }

    // run the SEQ_STATS process
    SEQ_STATS_LIB ( ch_seq_stats )
    ch_versions = ch_versions.mix(SEQ_STATS_LIB.out.versions)
    // Concatenate all output files
    def seq_stats_lib = SEQ_STATS_LIB.out.stats_txt
        .map { it[1] }
        .collectFile(name: "${workflow_name}_lib_stats", keepHeader: true, skip: 1, sort: true)

    // sort the stats output file
    SORT_STATS_LIB ( seq_stats_lib )

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Generating seq stats per sample
    ////////////////////////////////////////////////////////////////////////////////////////////////

    // Prepare channels for sample statistics
    ch_reads_sample = reads.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0]]
        new_meta.remove('single_end')
        [ new_meta, data[0] ]
    }.groupTuple()

    ch_fastp_log_sample = fastp_log.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_raw_bam_flagstat_sample = raw_bam_flagstat.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_mq_filtered_bam_flagstat_sample = mq_filtered_bam_flagstat.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_dedup_sample_flagstat = dedup_sample_flagstat.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_dedup_sample = dedup_sample.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }.groupTuple()

    ch_dpstats = dpstats.map { meta, data ->
        def new_meta = meta + [id: meta.id.split("_")[0]]
        new_meta.remove('single_end')
        [ new_meta, data ]
    }

    // merging different channels
    ch_seq_stats_sample = ch_reads_sample
        .join(ch_fastp_log_sample)
        .join(ch_raw_bam_flagstat_sample)
        .join(ch_mq_filtered_bam_flagstat_sample)
        .join(ch_dedup_sample_flagstat)
        .join(ch_dedup_sample)
        .join(ch_dpstats)

    // If competitive reference is used, include decoy flagstat
    if (params.competitive_reference) {
        ch_decoy_flagstat_sample = decoy_flagstat.map { meta, data ->
            def new_meta = meta + [id: meta.id.split("_")[0]]
            new_meta.remove('single_end')
            [ new_meta, data ]
        }.groupTuple()

        ch_seq_stats_sample = ch_seq_stats_sample.join(ch_decoy_flagstat_sample)
    } else {
        // If no decoy flagstat is provided, append a dummy path directly
        ch_seq_stats_sample = ch_seq_stats_sample.map {
            it + [ [file('/dev/null')] ]
        }
    }

    // run the SEQ_STATS process
    SEQ_STATS_SAMPLE ( ch_seq_stats_sample )
    ch_versions = ch_versions.mix(SEQ_STATS_SAMPLE.out.versions)
    // Concatenate all output files
    def seq_stats_sample = SEQ_STATS_SAMPLE.out.stats_txt
        .map { it[1] }
        .collectFile(name: "${workflow_name}_sample_stats", keepHeader: true, skip: 1, sort: true)

    // sort the stats output file
    SORT_STATS_SAMPLE ( seq_stats_sample )

    // Emit channels
    emit:
    versions            = ch_versions                           // channel: [ versions.yml ]
}