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

    FASTP ( reads )
    ch_versions = ch_versions.mix(FASTP.out.versions)

    if ( params.kraken2 == 'true' ) {
        KRAKEN2 ( FASTP.out.reads, kraken2_db )
        ch_versions = ch_versions.mix(KRAKEN2.out.versions)

        ch_kraken2_filtering = FASTP.out.reads.join(KRAKEN2.out.kraken2_output)

        FILTER_FASTQ ( ch_kraken2_filtering )
        ch_versions = ch_versions.mix(FILTER_FASTQ.out.versions)
    }


    emit:
    reads          = params.kraken2 == 'true' ? FILTER_FASTQ.out.filtered_reads : FASTP.out.reads   // Output filtered reads if Kraken is enabled, otherwise pass FASTP reads.
    json           = FASTP.out.json                                                                 // channel: [ val(meta), [ reads ] ]
    fastp_log      = FASTP.out.log                                                                  // channel: [ val(meta), [ reads ] ]
    reads_unmerged_R1 = FASTP.out.reads_unmerged_R1                                                 // channel: [ val(meta), [ reads ] ]
    reads_unmerged_R2 = FASTP.out.reads_unmerged_R2                                                 // channel: [ val(meta), [ reads ] ]
    kraken2_output  = params.kraken2 == 'true' ? KRAKEN2.out.kraken2_output : Channel.empty()       // channel: [ val(meta), [ reads ] ]
    kraken2_report  = params.kraken2 == 'true' ? KRAKEN2.out.kraken2_report : Channel.empty()       // channel: [ val(meta), [ reads ] ]
    versions       = ch_versions                                                                    // channel: [ versions.yml ]
}