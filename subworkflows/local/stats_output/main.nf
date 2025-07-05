#! /usr/bin/env nextflow

include { SEQ_STATS      } from '../../../modules/local/stats_output/seq_stats.nf'

workflow STATS_OUTPUT {
    take:
    reads
    fastp_log
    raw_bam_flagstat
    mq_filtered_bam_flagstat
    dedup_lib_flagstat
    dedup_lib

    main:
    ch_versions = Channel.empty()

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

    //ch_dedup_lib_flagstat = dedup_lib_flagstat.map { meta, data ->
    //    [ meta + ['id': meta.id], data ]
    //}

    // merging different channels
    ch_seq_stats = ch_reads
        .join(ch_fastp_log)
        .join(ch_raw_bam_flagstat)
        .join(ch_mq_filtered_bam_flagstat)
        .join(dedup_lib_flagstat)
        .join(dedup_lib)

    // run the SEQ_STATS process
    SEQ_STATS ( ch_seq_stats )
    ch_versions = ch_versions.mix(SEQ_STATS.out.versions)
    // Concatenate all output files
    def all_stats_txt = SEQ_STATS.out.stats_txt
        .map { it[1] }
        .collectFile(name: "${params.workflowname}_sequencing_stats.txt", keepHeader: true, skip: 1, storeDir: "${params.outdir}/stats")

    // Emit channels
    emit:
    versions        = ch_versions                // channel: [ versions.yml ]
}