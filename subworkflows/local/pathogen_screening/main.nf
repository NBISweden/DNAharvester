#! /usr/bin/env nextflow

include { SAMTOOLS_UNMAPPED_READS   } from '../../../modules/local/samtools/samtools_unmapped_reads.nf'
include { BOWTIE2_BUILD             } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BOWTIE2 as PS_BOWTIE2     } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_INDEX            } from '../../../modules/nf-core/samtools/index/main'
include { FILTERBAM                 } from '../../../modules/local/stats_output/filterbam.nf'
include { FILTERBAM_PLOT            } from '../../../modules/local/pathogen_screening/filterbam_plot.nf'


workflow PATHOGEN_SCREENING {
    take:
    raw_bam
    reference_database

    main:
    ch_versions             = Channel.empty()
    // define workflow name to be used in summary stats output
    def workflow_name       = params.workflow_run_name ?: workflow.runName


    // get unmapped reads from the raw BAM file
    SAMTOOLS_UNMAPPED_READS ( raw_bam )
    ch_versions             = ch_versions.mix(SAMTOOLS_UNMAPPED_READS.out.versions)
    ch_unmapped_reads       = SAMTOOLS_UNMAPPED_READS.out.unmapped_fastq


    // Index the pathogen reference database if it is not already indexed
    BOWTIE2_BUILD (reference_database, file(params.pathogen_reference_database).getParent())
    ch_versions             = ch_versions.mix(BOWTIE2_BUILD.out.versions)


    // Map the reads to the reference genome
    PS_BOWTIE2 ( ch_unmapped_reads, BOWTIE2_BUILD.out.index_dir )
    ch_versions             = ch_versions.mix(PS_BOWTIE2.out.versions)
    ch_bam                  = PS_BOWTIE2.out.bam
    // Index the BAM file
    SAMTOOLS_INDEX ( ch_bam )
    ch_versions             = ch_versions.mix ( SAMTOOLS_INDEX.out.versions )


    // run filterBAM to generate stats
    FILTERBAM ( ch_bam )
    ch_versions             = ch_versions.mix ( FILTERBAM.out.versions )


    // generate plot from filterBAM output
    FILTERBAM_PLOT ( ch_bam, FILTERBAM.out.filterBAM_stats )
    ch_versions             = ch_versions.mix ( FILTERBAM_PLOT.out.versions )


    emit:
    unmapped_fastq          = SAMTOOLS_UNMAPPED_READS.out.unmapped_fastq  // channel: [ val(meta), [ fastq ] ]
    reference_index         = BOWTIE2_BUILD.out.index_dir                 // channel: path(index)
    raw_bam                 = PS_BOWTIE2.out.bam                          // channel: [ val(meta), [ bam ] ]
    raw_bam_bai             = SAMTOOLS_INDEX.out.bai                      // channel: [ val(meta), [ raw_bam_bai ] ]
    filterBAM_stats         = FILTERBAM.out.filterBAM_stats               // channel: [ val(meta), [ filterBAM.csv ] ]
    pathogen_screening_plot = FILTERBAM_PLOT.out.pathogen_screening_plot  // channel: [ val(meta), [ pathogen_screening_plot.pdf ] ]
    versions                = ch_versions                                 // channel: [ versions.yml ]
}