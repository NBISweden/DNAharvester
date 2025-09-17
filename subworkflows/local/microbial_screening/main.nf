#! /usr/bin/env nextflow

include { SAMTOOLS_UNMAPPED_READS                   } from '../../../modules/local/samtools/samtools_unmapped_reads.nf'
include { BWA_INDEX as PS_BWA_INDEX                 } from '../../../modules/local/bwa/index.nf'
include { BOWTIE2_BUILD as PS_BOWTIE2_BUILD         } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BWA_ALN as PS_BWA_ALN                     } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE as PS_BWA_SAMSE                 } from '../../../modules/local/bwa/samse.nf'
include { BWA_ALN_MEM as PS_BWA_ALN_MEM             } from '../../../modules/local/bwa/bwa_aln_mem.nf'
include { BOWTIE2 as PS_BOWTIE2                     } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_INDEX                            } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW_MQ as PS_SAMTOOLS_VIEW_MQ   } from '../../../modules/local/samtools/samtools_view_mq.nf'
include { FILTERBAM as PS_FILTERBAM                 } from '../../../modules/local/stats_output/filterbam.nf'
include { FILTERBAM_PLOT as PS_FILTERBAM_PLOT       } from '../../../modules/local/pathogen_screening/filterbam_plot.nf'


workflow MICROBIAL_SCREENING {
    take:
    raw_bam
    ms_reference_database

    main:
    ch_versions             = Channel.empty()

    // get unmapped reads from the raw BAM file
    SAMTOOLS_UNMAPPED_READS ( raw_bam )
    ch_versions             = ch_versions.mix(SAMTOOLS_UNMAPPED_READS.out.versions)
    ch_unmapped_reads       = SAMTOOLS_UNMAPPED_READS.out.unmapped_fastq

    ////////////////////////////////////////////////////////////////////////////
    // Index the reference database if the index files are not present in the dir
    ////////////////////////////////////////////////////////////////////////////

    if (params.ms_mapping_tool == 'bwa-aln' || params.ms_mapping_tool == 'bwa-aln-mem') {
        // Build the BWA index
        PS_BWA_INDEX (ms_reference_database, file(params.ms_reference_database).getParent())
        ch_versions         = ch_versions.mix(PS_BWA_INDEX.out.versions)
        ch_ms_reference_database_index  = PS_BWA_INDEX.out.index_dir
    } else if (params.ms_mapping_tool == 'bowtie2') {
        // Build the Bowtie2 index
        PS_BOWTIE2_BUILD (ms_reference_database, file(params.ms_reference_database).getParent())
        ch_versions         = ch_versions.mix(PS_BOWTIE2_BUILD.out.versions)
        ch_ms_reference_database_index  = PS_BOWTIE2_BUILD.out.index_dir
    } else {
        error "Invalid mapping tool specified: ${params.ms_mapping_tool}. Use 'bwa-aln', 'bwa-aln-mem', or 'bowtie2'."
    }

    ////////////////////////////////////////////////////////////////////////////
    // Map the unmapped reads to the microbial reference database
    ////////////////////////////////////////////////////////////////////////////

    if (params.ms_mapping_tool == 'bwa-aln') {
        PS_BWA_ALN ( ch_unmapped_reads, ch_ms_reference_database_index )
        ch_versions         = ch_versions.mix(PS_BWA_ALN.out.versions)
        ch_bwa_samse        = ch_unmapped_reads.join(PS_BWA_ALN.out.sai)

        PS_BWA_SAMSE ( ch_bwa_samse, ch_ms_reference_database_index )
        ch_versions         = ch_versions.mix(PS_BWA_SAMSE.out.versions)
        ch_raw_bam              = PS_BWA_SAMSE.out.bam

    } else if (params.ms_mapping_tool == 'bwa-aln-mem') {
        PS_BWA_ALN_MEM ( ch_unmapped_reads, ch_ms_reference_database_index )
        ch_versions         = ch_versions.mix(PS_BWA_ALN_MEM.out.versions)
        ch_raw_bam              = PS_BWA_ALN_MEM.out.bam
    } else if (params.ms_mapping_tool == 'bowtie2') {
        PS_BOWTIE2 ( ch_unmapped_reads, ch_ms_reference_database_index )
        ch_versions         = ch_versions.mix(PS_BOWTIE2.out.versions)
        ch_raw_bam              = PS_BOWTIE2.out.bam
    } else {
        error "Invalid mapping tool specified: ${params.ms_mapping_tool}. Use 'bwa-aln' or 'bwa-aln-mem' or 'bowtie2'."
    }

    ////////////////////////////////////////////////////////////////////////////
    // Filter and Merge the BAM files
    ////////////////////////////////////////////////////////////////////////////

    // mapping quality filter
    PS_SAMTOOLS_VIEW_MQ ( ch_raw_bam )
    ch_versions             = ch_versions.mix ( PS_SAMTOOLS_VIEW_MQ.out.versions )


    // Index the BAM file
    SAMTOOLS_INDEX ( PS_SAMTOOLS_VIEW_MQ.out.bam )
    ch_bam_bai              = PS_SAMTOOLS_VIEW_MQ.out.bam.join(SAMTOOLS_INDEX.out.bai)
    ch_versions             = ch_versions.mix ( SAMTOOLS_INDEX.out.versions )

    // run filterBAM to generate stats
    PS_FILTERBAM ( ch_bam_bai )
    ch_versions             = ch_versions.mix ( PS_FILTERBAM.out.versions )


    // generate plot from filterBAM output
    PS_FILTERBAM_PLOT ( ch_bam_bai, PS_FILTERBAM.out.filterBAM_stats )
    ch_versions             = ch_versions.mix ( PS_FILTERBAM_PLOT.out.versions )


    emit:
    unmapped_fastq                = SAMTOOLS_UNMAPPED_READS.out.unmapped_fastq  // channel: [ val(meta), [ fastq ] ]
    ms_reference_database_index   = ch_ms_reference_database_index              // channel: path(index)
    bam                           = ch_raw_bam                                  // channel: [ val(meta), [ bam ] ]
    ms_filterBAM_stats            = PS_FILTERBAM.out.filterBAM_stats            // channel: [ val(meta), [ filterBAM.csv ] ]
    ms_screening_plot             = PS_FILTERBAM_PLOT.out.plot                  // channel: [ val(meta), [ pathogen_screening_plot.pdf ] ]
    versions                      = ch_versions                                 // channel: [ versions.yml ]
}