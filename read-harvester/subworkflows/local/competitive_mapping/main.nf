#! /usr/bin/env nextflow

include { BWA_ALN as BWA_ALN_COMPETITIVE                 } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE BWA_SAMSE_COMPETITIVE                } from '../../../modules/local/bwa/samse.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_COMPETITIVE   } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_COMPETITIVE   } from '../../../modules/local/samtools/faidx/main'
include { FAI_TO_BED as FAI_TO_BED_COMPETITIVE           } from '../../../modules/local/fai2bed/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_TARGET         } from '../../../modules/local/samtools/faidx/main'
include { FAI_TO_BED as FAI_TO_BED_TARGET                 } from '../../../modules/local/fai2bed/main'
include { BEDTOOLS_SUBTRACT as BEDTOOLS_SUBTRACT_TARGET   } from '../../../modules/local/bedtools/intersect/main'
include { SAMTOOLS_VIEW_REGIONS as SAMTOOLS_VIEW_DECOY    } from '../../../modules/local/samtools/view_regions/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DECOY          } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_DECOY             } from '../../../modules/nf-core/samtools/flagstat/main'
include { MULTIQC as MULTIQC_DECOY                        } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_VIEW_REGIONS as SAMTOOLS_VIEW_TARGET   } from '../../../modules/local/samtools/view_regions/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_TARGET         } from '../../../modules/nf-core/samtools/index/main'

workflow COMPETITIVE_MAPPING {
    take:
    competitive_reference
    reference
    reads

    main:
    ch_versions = Channel.empty()

    // Create channel for bwa index files
    ch_bwa_index_competitive         = 

    // Map the reads to the concatenated fasta file
    BWA_ALN_COMPETITIVE ( reads, ch_bwa_index_competitive )
    ch_versions                      = ch_versions.mix(BWA_ALN_COMPETITIVE.out.versions)

    BWA_SAMSE_COMPETITIVE ( BWA_ALN_COMPETITIVE.out.reads, BWA_ALN_COMPETITIVE.out.sai, ch_bwa_index_competitive )
    ch_versions                      = ch_versions.mix(BWA_SAMSE_COMPETITIVE.out.versions)

    // Index the BAM file
    SAMTOOLS_INDEX_COMPETITIVE ( BWA_SAMSE_COMPETITIVE.out.bam )
    ch_versions                      = ch_versions.mix(SAMTOOLS_INDEX_COMPETITIVE.out.versions)

    // Split the BAM file into target genome and decoy genome
    ch_concatenated_bam_index        = BWA_SAMSE_COMPETITIVE.out.bam.mix( SAMTOOLS_INDEX_COMPETITIVE.out.bai )

    // Generate *.fai index for the concatenated reference
    SAMTOOLS_FAIDX_COMPETITIVE ( competitive_reference )
    ch_versions                      = ch_versions.mix(SAMTOOLS_FAIDX_COMPETITIVE.out.versions)

    // Convert to *.bed format
    FAI_TO_BED_COMPETITIVE ( SAMTOOLS_FAIDX_COMPETITIVE.out.fai )

    // Generate *.fai index for the target reference genome
    SAMTOOLS_FAIDX_TARGET ( reference )
    ch_versions                      = ch_versions.mix(SAMTOOLS_FAIDX_TARGET.out.versions)

    // Convert to *.bed format
    FAI_TO_BED_TARGET ( SAMTOOLS_FAIDX_TARGET.out.fai )
   
    // Decoy genome
    // Extract the decoy genome chromosomes from the concatenated genome BED file
    ch_bedtools_subtract_target_intervals  = FAI_TO_BED_COMPETITIVE.out.bed.combine( FAI_TO_BED_TARGET.out.bed )
    BEDTOOLS_SUBTRACT_TARGET ( ch_bedtools_subtract_target_intervals )

    // Extract the region from the BAM file
    SAMTOOLS_VIEW_DECOY ( ch_concatenated_bam_index, competitive_reference, BEDTOOLS_SUBTRACT_TARGET.out.bed )
    ch_versions                      = ch_versions.mix(SAMTOOLS_VIEW_DECOY.out.versions)
     // Index the BAM file
    SAMTOOLS_INDEX_DECOY ( SAMTOOLS_VIEW_DECOY.out.bam )
    ch_versions                      = ch_versions.mix(SAMTOOLS_INDEX_DECOY.out.versions)   

    // Run samtools flagstat and MultiQC
    ch_flagstat_decoy                = SAMTOOLS_VIEW_DECOY.out.bam.join(SAMTOOLS_INDEX_DECOY.out.bai)
    FLAGSTAT_DECOY ( ch_flagstat_decoy )
    ch_versions                      = ch_versions.mix(FLAGSTAT_DECOY.out.versions)

    ch_multiqc_decoy_files           = FLAGSTAT_DECOY.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
    ch_multiqc_config                = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config          = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                  = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    MULTIQC_DECOY (
        ch_multiqc_decoy_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                      = ch_versions.mix(MULTIQC_DECOY.out.versions)

    // Target genome
    // Extract the region from the BAM file
    SAMTOOLS_VIEW_TARGET ( ch_concatenated_bam_index, competitive_reference, FAI_TO_BED_TARGET.out.bed )
    ch_versions                      = ch_versions.mix(SAMTOOLS_VIEW_TARGET.out.versions)
    // Index the BAM file containing only the target genome
    SAMTOOLS_INDEX_TARGET ( SAMTOOLS_VIEW_TARGET.out.bam )
    ch_versions                      = ch_versions.mix(SAMTOOLS_INDEX_TARGET.out.versions)

    emit:
    fai                              = SAMTOOLS_FAIDX_TARGET.out.fai         // channel: path(index)
    multiqc_decoy_report             = MULTIQC_DECOY.out.report.toList()     // channel: [ val(meta), path(report) ]
    bam                              = SAMTOOLS_VIEW_TARGET.out.bam          // channel: [ val(meta), [ bam ] ]
    bai                              = SAMTOOLS_INDEX_TARGET.out.bai         // channel: [ val(meta), [ bai ] ]
    versions                         = ch_versions                           // channel: [ versions.yml ]
}