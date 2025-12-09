#! /usr/bin/env nextflow

include { FASTP         }       from '../../../modules/local/fastp/main.nf'
include { KRAKEN2       }       from '../../../modules/local/kraken2/kraken2.nf'
include { FILTER_FASTQ  }       from '../../../modules/local/kraken2/filter_fastq.nf'
include { ADAPTCLEAN    }       from '../../../modules/local/adaptclean/adaptclean.nf'


workflow FASTQ_PROCESSING {
    take:
    raw_reads
    kraken2_database

    main:
    ch_versions = Channel.empty()

    //////////////////////////////////////////////////////////////////////////////////////
    // 1. FASTQ Processing: Adapter trimming and read merging.
    //////////////////////////////////////////////////////////////////////////////////////

    ch_adapter_removed_reads = Channel.empty()

    // Run FASTP for adapter trimming and read merging
    FASTP ( raw_reads )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // If merge_reads set to 'true', update meta.single_end to true for the merged reads output
    if (params.merge_reads.toBoolean()) {
        ch_merged_reads = FASTP.out.reads_merged.map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.single_end = true
            [new_meta, reads]
        }
        ch_adapter_removed_reads = ch_adapter_removed_reads.mix( ch_merged_reads )
    }

    // Run AdaptClean to remove leftover adapter sequences from SE reads and PE reads and PE unmerged reads
    // Note: Merge reads dont need AdaptClean as the adapters are completely removed during merging
    if (params.adaptclean.toBoolean()) {
        // Prepare input channel for AdaptClean
        def ch_reads_se = FASTP.out.reads_se.join(FASTP.out.json)
        def ch_reads_pe = FASTP.out.reads_pe.join(FASTP.out.json)

        // Channel for unmerged reads if keep_unmerged_reads is set to 'true'
        if (params.keep_unmerged_reads.toBoolean()) {
            // Join unmerged reads with their respective JSON files
            def ch_reads_unmerged = FASTP.out.reads_unmerged
                .join(FASTP.out.json)
                .map { meta, reads, json ->
                    def new_meta = meta.clone()
                    // Update meta for unmerged reads to distinguish them
                    new_meta.id = meta.id + "-unmerged"
                    [new_meta, reads, json]
            }

            ch_adaptclean_input = ch_reads_se.mix(ch_reads_pe).mix(ch_reads_unmerged)
        } else {
            ch_adaptclean_input = ch_reads_se.mix(ch_reads_pe)
        }
        // Run AdaptClean
        ADAPTCLEAN ( ch_adaptclean_input )
        ch_versions = ch_versions.mix(ADAPTCLEAN.out.versions)

        // Add adaptclean reads to the ch_adapter_removed_reads channel
        ch_adapter_removed_reads = ch_adapter_removed_reads.mix(ADAPTCLEAN.out.adaptclean_reads)

    } else {
        if (params.keep_unmerged_reads.toBoolean()) {
            // If AdaptClean is not run, and keep_unmerged_reads is true, updated meta for unmerged reads
            def ch_unmerged_reads = FASTP.out.reads_unmerged.map { meta, reads ->
                def new_meta = meta.clone()
                new_meta.id = meta.id + "-unmerged"
                [new_meta, reads]
            }
            ch_adapter_removed_reads = ch_adapter_removed_reads
                .mix(FASTP.out.reads_se)
                .mix(FASTP.out.reads_pe)
                .mix(ch_unmerged_reads)
        }
        else {
            ch_adapter_removed_reads = ch_adapter_removed_reads
                .mix(FASTP.out.reads_se)
                .mix(FASTP.out.reads_pe)
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // 2. Kraken2 Classification and Filtering
    //////////////////////////////////////////////////////////////////////////////////////

    // If Kraken2 classification is enabled, run Kraken2 and filter out classified reads
    if ( params.kraken2_filtering.toBoolean() ) {
        KRAKEN2 ( ch_adapter_removed_reads, kraken2_database )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions)

        ch_kraken2_filtering = ch_adapter_removed_reads.join(KRAKEN2.out.kraken2_output)

        FILTER_FASTQ ( ch_kraken2_filtering )
        ch_versions = ch_versions.mix(FILTER_FASTQ.out.versions)
    }
    //////////////////////////////////////////////////////////////////////////////////////

    // Determine final processed reads based on Kraken2 filtering set or not
    def ch_processed_reads = params.kraken2_filtering.toBoolean() ? FILTER_FASTQ.out.filtered_reads : ch_adapter_removed_reads

    emit:
    processed_reads     = ch_processed_reads                                                                    // channel: [ val(meta), [ reads ] ]
    json                = FASTP.out.json                                                                        // channel: [ val(meta), [ json ] ]
    fastp_log           = FASTP.out.log                                                                         // channel: [ val(meta), [ log ] ]
    kraken2_output      = params.kraken2_filtering.toBoolean() ? KRAKEN2.out.kraken2_output : Channel.empty()   // channel: [ val(meta), [ kraken2_output ] ]
    kraken2_report      = params.kraken2_filtering.toBoolean() ? KRAKEN2.out.kraken2_report : Channel.empty()   // channel: [ val(meta), [ kraken2_report ] ]
    versions            = ch_versions                                                                           // channel: [ versions.yml ]
}