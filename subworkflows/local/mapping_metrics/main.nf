#! /usr/bin/env nextflow

include { MAPPING_METRICS_LIB       as MM_MAPPING_METRICS_LIB       } from '../../../modules/local/mapping_metrics/mapping_metrics_lib.nf'
include { SORT_METRICS              as MM_SORT_METRICS_LIB          } from '../../../modules/local/mapping_metrics/sort_metrics.nf'
include { MAPPING_METRICS_SAMPLE    as MM_MAPPING_METRICS_SAMPLE    } from '../../../modules/local/mapping_metrics/mapping_metrics_sample.nf'
include { SORT_METRICS              as MM_SORT_METRICS_SAMPLE       } from '../../../modules/local/mapping_metrics/sort_metrics.nf'

workflow MAPPING_METRICS {
    take:
    workflow_name
    fastp_json
    raw_bam_flagstat
    filtered_bam_flagstat
    dedup_lib_flagstat
    dedup_lib
    dedup_sample_flagstat
    dedup_sample
    decoy_flagstat
    dpstats

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Generating mapping metrics per library
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

    ch_filtered_bam_flagstat_lib = filtered_bam_flagstat.map { meta, data ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] + "_" + meta.id.split("_")[1]
        new_meta.remove('single_end')
        new_meta.remove('read_group')
        [ new_meta, data ]
    }.groupTuple()

    // merging channels for library statistics
    ch_mapping_metrics_lib = ch_fastp_json_lib
        .join(ch_raw_bam_flagstat_lib)
        .join(ch_filtered_bam_flagstat_lib)
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
        // join decoy flagstat channel to mapping metrics lib channel
        ch_mapping_metrics_lib = ch_mapping_metrics_lib.join(ch_decoy_flagstat_lib)
    } else {
        // If no decoy flagstat is provided, append a dummy path until nextflow supports optional inputs :(
        ch_mapping_metrics_lib = ch_mapping_metrics_lib.map {
            it + [ [file('/dev/null')] ]
        }
    }

    // run the MAPPING_METRICS process
    MM_MAPPING_METRICS_LIB ( ch_mapping_metrics_lib )
    ch_versions = ch_versions.mix(MM_MAPPING_METRICS_LIB.out.versions)
    // Concatenate all output files
    def ch_mapping_metrics_lib_all = MM_MAPPING_METRICS_LIB.out.stats_tsv
        .map { it[1] }
        .collectFile(name: "${workflow_name}_library_metrics", keepHeader: true, skip: 1, sort: true)

    // sort the library metrics output file
    MM_SORT_METRICS_LIB ( ch_mapping_metrics_lib_all )

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Generating mapping metrics per sample
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

    ch_filtered_bam_flagstat_sample = ch_filtered_bam_flagstat_lib.map { meta, data ->
        def new_meta = meta.clone()
        new_meta.id = meta.id.split("_")[0] // remove library_id
        new_meta.remove('library_type') // remove library_type same sample can have multiple library_type
        [ new_meta, data ]
    }.groupTuple().map { meta, data -> [ meta, data.flatten() ] }

    // merging§ channels for sample statistics
    ch_mapping_metrics_sample = ch_fastp_json_sample
        .join(ch_raw_bam_flagstat_sample)
        .join(ch_filtered_bam_flagstat_sample)
        .join(dedup_sample_flagstat) // dedup_sample_flagstat is already in the corect format
        .join(dedup_sample) // dedup_sample is already in the corect format
        .join(dpstats) // dpstats is already in the corect format

    // If competitive reference is used, include decoy flagstat
    if (params.competitive_reference) {
        ch_decoy_flagstat_sample = ch_decoy_flagstat_lib.map { meta, data ->
            def new_meta = meta + [id: meta.id.split("_")[0]]
            new_meta.remove('library_type')
            [ new_meta, data ]
        }.groupTuple().map { meta, data -> [ meta, data.flatten() ] }

        ch_mapping_metrics_sample = ch_mapping_metrics_sample.join(ch_decoy_flagstat_sample)
    } else {
        // If no decoy flagstat is provided, append a dummy path directly
        ch_mapping_metrics_sample = ch_mapping_metrics_sample.map {
            it + [ [file('/dev/null')] ]
        }
    }

    // run the MAPPING_METRICS process
    MM_MAPPING_METRICS_SAMPLE ( ch_mapping_metrics_sample )
    ch_versions = ch_versions.mix(MM_MAPPING_METRICS_SAMPLE.out.versions)
    // Concatenate all output files
    def ch_mapping_metrics_sample_all = MM_MAPPING_METRICS_SAMPLE.out.stats_tsv
        .map { it[1] }
        .collectFile(name: "${workflow_name}_sample_metrics", keepHeader: true, skip: 1, sort: true)

    // sort the sample metrics output file
    MM_SORT_METRICS_SAMPLE ( ch_mapping_metrics_sample_all )

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    versions            = ch_versions                           // channel: [ versions.yml ]
}