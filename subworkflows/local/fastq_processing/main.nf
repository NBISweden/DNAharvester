#! /usr/bin/env nextflow

include { FASTP         as FP_FASTP         }       from '../../../modules/local/fastp/main.nf'
include { LEEHOM        as FP_LEEHOM        }       from '../../../modules/local/leehom/leehom.nf'
include { KRAKEN2       as FP_KRAKEN2       }       from '../../../modules/local/kraken2/kraken2.nf'
include { ADAPTCLEAN    as FP_ADAPTCLEAN    }       from '../../../modules/local/adaptclean/adaptclean.nf'


workflow FASTQ_PROCESSING {
    take:
    raw_reads
    kraken2_database

    main:
    ch_versions  = Channel.empty()
    ch_trim_json = Channel.empty()
    ch_trim_log  = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. FASTQ Processing: Adapter trimming and read merging.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ch_adapter_removed_reads = Channel.empty()

    if (params.adapter_trimming_tool == 'leehom') {

        // Run leeHom for adapter trimming and read merging (ancient-DNA-focused alternative to fastp)
        FP_LEEHOM ( raw_reads )
        ch_versions  = ch_versions.mix(FP_LEEHOM.out.versions)
        ch_trim_json = FP_LEEHOM.out.json
        ch_trim_log  = FP_LEEHOM.out.log

        // leeHom always attempts to merge overlapping read pairs. Update meta.single_end
        // to true for the merged reads output, same as for fastp's merged reads
        ch_merged_reads = FP_LEEHOM.out.reads_merged.map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.single_end = true
            [new_meta, reads]
        }
        ch_adapter_removed_reads = ch_adapter_removed_reads
            .mix(FP_LEEHOM.out.reads_se)
            .mix(ch_merged_reads)

        // Note: leeHom already fully trims adapters from unmerged reads, so AdaptClean
        // is not applicable here (unlike for fastp's unmerged/non-merged reads)
        if (params.keep_unmerged_reads.toBoolean()) {
            def ch_unmerged_reads = FP_LEEHOM.out.reads_unmerged.map { meta, reads ->
                def new_meta = meta.clone()
                new_meta.id = meta.id + "-unmerged"
                [new_meta, reads]
            }
            ch_adapter_removed_reads = ch_adapter_removed_reads.mix(ch_unmerged_reads)
        }

    } else {

        // Run FASTP for adapter trimming and read merging
        FP_FASTP ( raw_reads )
        ch_versions  = ch_versions.mix(FP_FASTP.out.versions)
        ch_trim_json = FP_FASTP.out.json
        ch_trim_log  = FP_FASTP.out.log

        // If merge_reads set to 'true', update meta.single_end to true for the merged reads output
        if (params.merge_reads.toBoolean()) {
            ch_merged_reads = FP_FASTP.out.reads_merged.map { meta, reads ->
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
            def ch_reads_se = FP_FASTP.out.reads_se.join(FP_FASTP.out.json)
            def ch_reads_pe = FP_FASTP.out.reads_pe.join(FP_FASTP.out.json)

            // Channel for unmerged reads if keep_unmerged_reads is set to 'true'
            if (params.keep_unmerged_reads.toBoolean()) {
                // Join unmerged reads with their respective JSON files
                def ch_reads_unmerged = FP_FASTP.out.reads_unmerged
                    .join(FP_FASTP.out.json)
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
            FP_ADAPTCLEAN ( ch_adaptclean_input )
            ch_versions = ch_versions.mix(FP_ADAPTCLEAN.out.versions)

            // Add adaptclean reads to the ch_adapter_removed_reads channel
            ch_adapter_removed_reads = ch_adapter_removed_reads.mix(FP_ADAPTCLEAN.out.adaptclean_reads)

        } else {
            if (params.keep_unmerged_reads.toBoolean()) {
                // If AdaptClean is not run, and keep_unmerged_reads is true, updated meta for unmerged reads
                def ch_unmerged_reads = FP_FASTP.out.reads_unmerged.map { meta, reads ->
                    def new_meta = meta.clone()
                    new_meta.id = meta.id + "-unmerged"
                    [new_meta, reads]
                }
                ch_adapter_removed_reads = ch_adapter_removed_reads
                    .mix(FP_FASTP.out.reads_se)
                    .mix(FP_FASTP.out.reads_pe)
                    .mix(ch_unmerged_reads)
            }
            else {
                ch_adapter_removed_reads = ch_adapter_removed_reads
                    .mix(FP_FASTP.out.reads_se)
                    .mix(FP_FASTP.out.reads_pe)
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Kraken2 Classification and Filtering
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    def ch_processed_reads
    // If Kraken2 classification is enabled, run Kraken2 and filter out classified reads
    if ( params.kraken2_filtering.toBoolean() ) {
        FP_KRAKEN2 ( ch_adapter_removed_reads, kraken2_database )
        ch_versions = ch_versions.mix(FP_KRAKEN2.out.versions)

        // Keep only unclassified reads for downstream processing
        ch_processed_reads = FP_KRAKEN2.out.kraken2_unclassified
    } else {
        ch_processed_reads = ch_adapter_removed_reads
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    processed_reads     = ch_processed_reads                                                                        // channel: [ val(meta), [ reads ] ]
    fastp_json          = ch_trim_json                                                                              // channel: [ val(meta), [ json ] ] (fastp or leehom, depending on adapter_trimming_tool)
    fastp_log           = ch_trim_log                                                                               // channel: [ val(meta), [ log ] ]  (fastp or leehom, depending on adapter_trimming_tool)
    kraken2_output      = params.kraken2_filtering.toBoolean() ? FP_KRAKEN2.out.kraken2_output : Channel.empty()    // channel: [ val(meta), [ kraken2_output ] ]
    kraken2_report      = params.kraken2_filtering.toBoolean() ? FP_KRAKEN2.out.kraken2_report : Channel.empty()    // channel: [ val(meta), [ kraken2_report ] ]
    versions            = ch_versions                                                                               // channel: [ versions.yml ]
}