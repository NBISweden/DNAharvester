#! /usr/bin/env nextflow

include { FASTP         }       from '../../../modules/local/fastp/main.nf'
include { KRAKEN2       }       from '../../../modules/local/kraken2/kraken2.nf'
include { FILTER_FASTQ  }       from '../../../modules/local/kraken2/filter_fastq.nf'


workflow FASTQ_PROCESSING {
    take:
    raw_reads
    kraken2_db

    main:
    ch_versions = Channel.empty()

    //////////////////////////////////////////////////////////////////////////////////////
    // 1. FASTQ Processing: Adapter trimming and read merging.
    //////////////////////////////////////////////////////////////////////////////////////

    // Run FASTP for adapter trimming and read merging
    FASTP ( raw_reads )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    // If merge_reads set to 'true', update meta.single_end to true for the merged reads output
    if (params.merge_reads.toBoolean()) {
        ch_processed_reads = FASTP.out.processed_reads.map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.single_end = true
            [new_meta, reads]
        }
    } else {
        ch_processed_reads = FASTP.out.processed_reads
    }

    // If keep_unmerged_reads is set to 'true', add "-unmerged" suffix to meta.id for unmerged reads
    if (params.keep_unmerged_reads.toBoolean()) {
        ch_unmerged_reads = FASTP.out.unmerged_reads.map { meta, reads ->
            def new_meta = meta.clone()
            new_meta.id = meta.id + '-unmerged'
            [ new_meta, reads ]
        }
        ch_processed_reads = ch_processed_reads.mix(ch_unmerged_reads)
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // 2. Kraken2 Classification and Filtering
    //////////////////////////////////////////////////////////////////////////////////////

    // If Kraken2 classification is enabled, run Kraken2 and filter out classified reads
    if ( params.kraken2_filtering.toBoolean() ) {
        KRAKEN2 ( ch_processed_reads, kraken2_db )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions)

        ch_kraken2_filtering = ch_processed_reads.join(KRAKEN2.out.kraken2_output)

        FILTER_FASTQ ( ch_kraken2_filtering )
        ch_versions = ch_versions.mix(FILTER_FASTQ.out.versions)
    }
    //////////////////////////////////////////////////////////////////////////////////////

    // Determine final processed reads based on Kraken2 filtering set or not
    def ch_processed_reads_final = params.kraken2_filtering.toBoolean() ? FILTER_FASTQ.out.filtered_reads : ch_processed_reads

    emit:
    processed_reads     = ch_processed_reads_final                                                              // channel: [ val(meta), [ reads ] ]
    json                = FASTP.out.json                                                                        // channel: [ val(meta), [ json ] ]
    fastp_log           = FASTP.out.log                                                                         // channel: [ val(meta), [ log ] ]
    kraken2_output      = params.kraken2_filtering.toBoolean() ? KRAKEN2.out.kraken2_output : Channel.empty()   // channel: [ val(meta), [ kraken2_output ] ]
    kraken2_report      = params.kraken2_filtering.toBoolean() ? KRAKEN2.out.kraken2_report : Channel.empty()   // channel: [ val(meta), [ kraken2_report ] ]
    versions            = ch_versions                                                                           // channel: [ versions.yml ]
}