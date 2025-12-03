#! /usr/bin/env nextflow

include { FASTP         }       from '../../../modules/local/fastp/main.nf'
include { KRAKEN2       }       from '../../../modules/local/kraken2/kraken2.nf'
include { FILTER_FASTQ  }       from '../../../modules/local/kraken2/filter_fastq.nf'


workflow FASTQ_PROCESSING {
    take:
    kraken2_db
    reads

    main:
    ch_versions = Channel.empty()

    // Run FASTP for adapter trimming and read merging
    FASTP ( reads )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // Update meta.single_end to true for the merged reads output
    ch_reads = FASTP.out.reads.map { meta, reads ->
        meta.single_end = true
        [meta, reads]
    }

    // If Kraken2 classification is enabled, run Kraken2 and filter out classified reads
    if ( params.kraken2.toBoolean() ) {
        KRAKEN2 ( ch_reads, kraken2_db )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions)

        ch_kraken2_filtering = ch_reads.join(KRAKEN2.out.kraken2_output)

        FILTER_FASTQ ( ch_kraken2_filtering )
        ch_versions = ch_versions.mix(FILTER_FASTQ.out.versions)
    }

    // Join unmerged reads to a single channel for output
    ch_unmerged_reads = FASTP.out.reads_unmerged_R1.join(FASTP.out.reads_unmerged_R2)

    emit:
    reads               = params.kraken2.toBoolean() ? FILTER_FASTQ.out.filtered_reads : ch_reads               // Output filtered reads if Kraken is enabled, otherwise pass FASTP reads.
    json                = FASTP.out.json                                                                        // channel: [ val(meta), [ reads ] ]
    fastp_log           = FASTP.out.log                                                                         // channel: [ val(meta), [ reads ] ]
    unmerged_reads      = ch_unmerged_reads                                                                     // channel: [ val(meta), [ reads ] ]
    kraken2_output      = params.kraken2.toBoolean() ? KRAKEN2.out.kraken2_output : Channel.empty()             // channel: [ val(meta), [ reads ] ]
    kraken2_report      = params.kraken2.toBoolean() ? KRAKEN2.out.kraken2_report : Channel.empty()             // channel: [ val(meta), [ reads ] ]
    versions            = ch_versions                                                                           // channel: [ versions.yml ]
}