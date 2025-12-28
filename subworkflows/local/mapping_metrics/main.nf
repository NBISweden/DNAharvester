#! /usr/bin/env nextflow

include { SEQ_STATS_LIB         as MM_SEQ_STATS_LIB       } from '../../../modules/local/mapping_metrics/seq_stats_lib.nf'
include { SORT_STATS            as MM_SORT_STATS_LIB      } from '../../../modules/local/mapping_metrics/sort_stats.nf'
include { SEQ_STATS_SAMPLE      as MM_SEQ_STATS_SAMPLE    } from '../../../modules/local/mapping_metrics/seq_stats_sample.nf'
include { SORT_STATS            as MM_SORT_STATS_SAMPLE   } from '../../../modules/local/mapping_metrics/sort_stats.nf'

workflow MAPPING_METRICS {
    take:
    workflow_name
    fastp_json
    raw_bam_flagstat
    mq_filtered_bam_flagstat
    dedup_lib_flagstat
    dedup_lib
    dedup_sample_flagstat
    dedup_sample
    decoy_flagstat
    dpstats

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Generating seq stats per library
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Prepare channels for library statistics
    ch_fastp_json_lib = fastp_json.map { meta, data ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] + "_" + meta.id.split("_")[1] // remove lane info
        new_meta.remove('single_end') // remove single_end info as same library can have both SE and PE
        new_meta.remove('read_group') // remove read_group info as same library can have multiple read groups
        [ new_meta, data ]
    }.groupTuple()

    ch_raw_bam_flagstat_lib = raw_bam_flagstat.map { meta, data ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] + "_" + meta.id.split("_")[1] // remove lane info
        new_meta.remove('single_end') // remove single_end info as same library can have both SE and PE
        new_meta.remove('read_group') // remove read_group info as same library can have multiple read groups
        [ new_meta, data ]
    }.groupTuple()

    ch_mq_filtered_bam_flagstat_lib = mq_filtered_bam_flagstat.map { meta, data ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] + "_" + meta.id.split("_")[1]
        new_meta.remove('single_end')
        new_meta.remove('read_group')
        [ new_meta, data ]
    }.groupTuple()

    // merging channels for library statistics
    ch_seq_stats_lib = ch_fastp_json_lib
        .join(ch_raw_bam_flagstat_lib)
        .join(ch_mq_filtered_bam_flagstat_lib)
        .join(dedup_lib_flagstat) // dedup_lib_flagstat is already in the corect format
        .join(dedup_lib) // dedup_lib is already in the corect format

    // If competitive reference is used, include decoy flagstat
    if (params.competitive_reference) {
        ch_decoy_flagstat_lib = decoy_flagstat.map { meta, data ->
            def new_meta = meta.clone()
            new_meta.id = meta.id.split("_")[0] + "_" + meta.id.split("_")[1]
            new_meta.remove('single_end')
            new_meta.remove('read_group')
            [ new_meta, data ]
        }.groupTuple()
        // join decoy flagstat channel to seq stats lib channel
        ch_seq_stats_lib = ch_seq_stats_lib.join(ch_decoy_flagstat_lib)
    } else {
        // If no decoy flagstat is provided, append a dummy path until nextflow supports optional inputs :(
        ch_seq_stats_lib = ch_seq_stats_lib.map {
            it + [ [file('/dev/null')] ]
        }
    }

    // run the SEQ_STATS process
    MM_SEQ_STATS_LIB ( ch_seq_stats_lib )
    ch_versions = ch_versions.mix(MM_SEQ_STATS_LIB.out.versions)
    // Concatenate all output files
    def ch_seq_stats_lib_all = MM_SEQ_STATS_LIB.out.stats_tsv
        .map { it[1] }
        .collectFile(name: "${workflow_name}_lib_stats", keepHeader: true, skip: 1, sort: true)

    // sort the stats output file
    MM_SORT_STATS_LIB ( ch_seq_stats_lib_all )

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Generating seq stats per sample
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Prepare channels for sample statistics
    ch_fastp_json_sample = ch_fastp_json_lib.map { meta, data ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] // remove library_id
        new_meta.remove('library_type') // remove library_type same sample can have multiple library_type
        [ new_meta, data ]
    }.groupTuple().map { meta, data -> [ meta, data.flatten() ] }

    ch_raw_bam_flagstat_sample = ch_raw_bam_flagstat_lib.map { meta, data ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] // remove library_id
        new_meta.remove('library_type') // remove library_type same sample can have multiple library_type
        [ new_meta, data ]
    }.groupTuple().map { meta, data -> [ meta, data.flatten() ] }

    ch_mq_filtered_bam_flagstat_sample = ch_mq_filtered_bam_flagstat_lib.map { meta, data ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] // remove library_id
        new_meta.remove('library_type') // remove library_type same sample can have multiple library_type
        [ new_meta, data ]
    }.groupTuple().map { meta, data -> [ meta, data.flatten() ] }

    // merging§ channels for sample statistics
    ch_seq_stats_sample = ch_fastp_json_sample
        .join(ch_raw_bam_flagstat_sample)
        .join(ch_mq_filtered_bam_flagstat_sample)
        .join(dedup_sample_flagstat) // dedup_sample_flagstat is already in the corect format
        .join(dedup_sample) // dedup_sample is already in the corect format
        .join(dpstats) // dpstats is already in the corect format

    // If competitive reference is used, include decoy flagstat
    if (params.competitive_reference) {
        ch_decoy_flagstat_sample = ch_decoy_flagstat_lib.map { meta, data ->
            def new_meta = meta + [id: meta.id.split("_")[0]]
            new_meta.remove('library_type')
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
    MM_SEQ_STATS_SAMPLE ( ch_seq_stats_sample )
    ch_versions = ch_versions.mix(MM_SEQ_STATS_SAMPLE.out.versions)
    // Concatenate all output files
    def ch_seq_stats_sample_all = MM_SEQ_STATS_SAMPLE.out.stats_tsv
        .map { it[1] }
        .collectFile(name: "${workflow_name}_sample_stats", keepHeader: true, skip: 1, sort: true)

    // sort the stats output file
    MM_SORT_STATS_SAMPLE ( ch_seq_stats_sample_all )

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    versions            = ch_versions                           // channel: [ versions.yml ]
}