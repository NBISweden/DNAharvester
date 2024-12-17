#! /usr/bin/env nextflow

include { SEQ_STATS      } from '../../../modules/local/stats_output/seq_stats.nf'

workflow STATS_OUTPUT {
    take:
    reads
    fastp_log
    raw_bam_flagstat
    dedup_lib_flagstat
    //realigned
    dedup_sample

    main:
    ch_versions = Channel.empty()

    // Processing channel for merging
    ch_reads = reads.map { meta, data ->
        [['id': meta.id.split("_")[0]], data[0] ]
    }
    ch_fastp_log = fastp_log.map { meta, data ->
        [['id': meta.id.split("_")[0]], data ]
    }
    ch_raw_bam_flagstat = raw_bam_flagstat.map { meta, data ->
        [['id': meta.id.split("_")[0]], data ]
    }
    ch_dedup_lib_flagstat = dedup_lib_flagstat.map { meta, data ->
        [['id': meta.id.split("_")[0]], data ]
    }
    // merging different channels
    ch_seq_stats = ch_reads
        .join(ch_fastp_log)
        .join(ch_raw_bam_flagstat)
        .join(ch_dedup_lib_flagstat)
        //.join(realigned)
        .join(dedup_sample)

    // run the SEQ_STATS process
    SEQ_STATS ( ch_seq_stats )

    // Concatenate all output files
    def all_stats_txt = SEQ_STATS.out.stats_txt
        .map { it[1] }
        .collectFile(name: "${workflow.runName}_sequencing_stats.txt", keepHeader: true, skip: 1, storeDir: "${params.outdir}/stats")

    // Emit channels
    emit:
    versions        = ch_versions                // channel: [ versions.yml ]
}