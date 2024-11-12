#! /usr/bin/env nextflow

include { SEQ_STATS      } from '../../../modules/local/stats_output/seq_stats.nf'

workflow STATS_OUTPUT {
    take:
    json
    raw_bam_flagstat
    dedup_lib_flagstat


    main:
    SEQ_STATS ( json, raw_bam_flagstat, dedup_lib_flagstat)

    emit:
    stats_txt    = SEQ_STATS.out.stats_txt                            // channel: [ val(meta), [ txt ] ]

}



