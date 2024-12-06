#! /usr/bin/env nextflow

include { CONCATENATE_TARGET_DECOY_REFS                 } from '../../../modules/local/concatenate_target_decoy_refs/main'
include { BWA_INDEX as BWA_INDEX_CONCATENATED           } from '../../../modules/local/bwa/index.nf'
include { BWA_ALN as BWA_ALN_CONCATENATED               } from '../../../modules/local/bwa/aln.nf'
include { BWA_SAMSE as BWA_SAMSE_CONCATENATED           } from '../../../modules/local/bwa/samse.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_CONCATENATED } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_DECOY        } from '../../../modules/local/samtools/faidx/main'
include { FAI_TO_BED as FAI_TO_BED_DECOY                } from '../../../modules/local/fai2bed/main'
include { SAMTOOLS_VIEW_REGIONS as SAMTOOLS_VIEW_DECOY  } from '../../../modules/local/samtools/view_regions/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DECOY        } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_TARGET       } from '../../../modules/local/samtools/faidx/main'
include { FAI_TO_BED as FAI_TO_BED_TARGET               } from '../../../modules/local/fai2bed/main'
include { SAMTOOLS_VIEW_REGIONS as SAMTOOLS_VIEW_TARGET } from '../../../modules/local/samtools/view_regions/main'
include { SAMTOOLS_FLAGSTAT as FLAGSTAT_DECOY           } from '../../../modules/nf-core/samtools/flagstat/main'
include { MULTIQC as MULTIQC_DECOY                      } from '../../../modules/nf-core/multiqc/main'

workflow MAPPING_CONCATENATED_REFS {
    take:
    reference
    decoy
    reads // merged paired-end reads or trimmed single-end reads

    main:
    ch_versions = Channel.empty()

    // Concatenate the two references
    CONCATENATE_TARGET_DECOY_REFS ( reference, decoy )

    // Index the concatented fasta file
    BWA_INDEX_CONCATENATED ( CONCATENATE_TARGET_DECOY_REFS.out.fasta )
    ch_versions                      = ch_versions.mix(BWA_INDEX_CONCATENATED.out.versions)

    // Map the reads to the concatenated fasta file
    BWA_ALN_CONCATENATED ( reads, BWA_INDEX_CONCATENATED.out.index )
    ch_versions                      = ch_versions.mix(BWA_ALN_CONCATENATED.out.versions)

    BWA_SAMSE_CONCATENATED ( BWA_ALN_CONCATENATED.out.reads, BWA_ALN_CONCATENATED.out.sai, BWA_INDEX_CONCATENATED.out.index )
    ch_versions                      = ch_versions.mix(BWA_SAMSE_CONCATENATED.out.versions)

    // Index the BAM file
    SAMTOOLS_INDEX_CONCATENATED ( BWA_SAMSE_CONCATENATED.out.bam )
    ch_versions                      = ch_versions.mix(SAMTOOLS_INDEX_CONCATENATED.out.versions)

    // Split the BAM file into target genome and decoy genome
    ch_concatenated_bam_index        = BWA_SAMSE_CONCATENATED.out.bam.mix( SAMTOOLS_INDEX_CONCATENATED.out.bai )

    // Decoy genome
    // Generate *.fai index
    SAMTOOLS_FAIDX_DECOY ( decoy )
    ch_versions                      = ch_versions.mix(SAMTOOLS_FAIDX_DECOY.out.versions)
    // Convert to *.bed format
    FAI_TO_BED_DECOY ( SAMTOOLS_FAIDX_DECOY.out.fai )
    // Extract the region from the BAM file
    SAMTOOLS_VIEW_DECOY ( ch_concatenated_bam_index, decoy, FAI_TO_BED_DECOY.out.bed )
    ch_versions                      = ch_versions.mix(SAMTOOLS_VIEW_DECOY.out.versions)
     // Index the BAM file
    SAMTOOLS_INDEX_DECOY ( SAMTOOLS_VIEW_DECOY.out.bam )
    ch_versions                      = ch_versions.mix(SAMTOOLS_INDEX_DECOY.out.versions)   

    // Run samtools flagstat and MultiQC
    ch_flagstat_decoy                = SAMTOOLS_VIEW_DECOY.out.bam.join(SAMTOOLS_INDEX_DECOY.out.bai)
    FLAGSTAT_DECOY ( ch_flagstat_decoy )
    ch_versions                      = ch_versions.mix(FLAGSTAT_DECOY.out.versions)

    ch_multiqc_decoy_files           = FLAGSTAT_DECOY.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
    MULTIQC_DECOY (
        ch_multiqc_decoy_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                      = ch_versions.mix(MULTIQC_DECOY.out.versions)

    // Target genome
    // Generate *.fai index
    SAMTOOLS_FAIDX_TARGET ( reference )
    ch_versions                      = ch_versions.mix(SAMTOOLS_FAIDX_TARGET.out.versions)
    // Convert to *.bed format
    FAI_TO_BED_TARGET ( SAMTOOLS_FAIDX_TARGET.out.fai )

    // Extract the region from the BAM file
    SAMTOOLS_VIEW_TARGET ( ch_concatenated_bam_index, reference, FAI_TO_BED_TARGET.out.bed )
    ch_versions                      = ch_versions.mix(SAMTOOLS_VIEW_TARGET.out.versions)
    // Index the BAM file containing only the target genome
    SAMTOOLS_INDEX_TARGET ( SAMTOOLS_VIEW_TARGET.out.bam )
    ch_versions                      = ch_versions.mix(SAMTOOLS_INDEX_CONCATENATED.out.versions)

    emit:
    multiqc_decoy_report             = MULTIQC_DECOY.out.report.toList()     // channel: [ val(meta), path(report) ]
    bam                              = SAMTOOLS_VIEW_TARGET.out.bam          // channel: [ val(meta), [ bam ] ]
    bai                              = SAMTOOLS_INDEX_TARGET.out.bai         // channel: [ val(meta), [ bai ] ]
    versions                         = ch_versions                           // channel: [ versions.yml ]
}