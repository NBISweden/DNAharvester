#! /usr/bin/env nextflow

// Analysis script for DNAharvester pipeline

nextflow.enable.dsl = 2

// Import subworkflows
include { GUNZIP                     } from "$projectDir/modules/local/gunzip/main"
include { INPUT_CHECK                } from "$projectDir/subworkflows/local/input_check/main"
include { RAW_FASTQ_QC               } from "$projectDir/subworkflows/local/raw_fastq_qc/main"
include { FASTQ_PROCESSING           } from "$projectDir/subworkflows/local/fastq_processing/main"
include { PROCESSED_FASTQ_QC         } from "$projectDir/subworkflows/local/processed_fastq_qc/main"
include { MAPPING                    } from "$projectDir/subworkflows/local/mapping/main"
include { COMPETITIVE_MAPPING        } from "$projectDir/subworkflows/local/competitive_mapping/main"
include { MICROBIAL_SCREENING        } from "$projectDir/subworkflows/local/microbial_screening/main"
include { ITERATIVE_ASSEMBLY         } from "$projectDir/subworkflows/local/iterative_assembly/main"
include { TAXONOMIC_CLASSIFICATION   } from "$projectDir/subworkflows/local/taxonomic_classification/main"
include { REPEAT_CPG_IDENTIFICATION  } from "$projectDir/subworkflows/local/repeat_cpg_identification/main"
include { RAW_BAM_QC                 } from "$projectDir/subworkflows/local/raw_bam_qc/main"
include { BAM_PROCESSING             } from "$projectDir/subworkflows/local/bam_processing/main"
include { PROCESSED_BAM_QC           } from "$projectDir/subworkflows/local/processed_bam_qc/main"
include { RANDOM_SAMPLING_BAM        } from "$projectDir/subworkflows/local/random_sampling_bam/main"
include { VARIANT_CALLING_BCFTOOLS   } from "$projectDir/subworkflows/local/variant_calling/variant_calling_bcftools.nf"
include { VARIANT_CALLING_ANGSD      } from "$projectDir/subworkflows/local/variant_calling/variant_calling_angsd.nf"
include { STATS_OUTPUT               } from "$projectDir/subworkflows/local/stats_output/main"


workflow {

    // Set the workflow name
    def workflow_name = params.workflow_run_name ?: workflow.runName

    // The primary workflow for the DNAharvester pipeline
    log.info("""
    Running DNAharvester. Workflow run name: $workflow_name
    """)

    ////////////////////////////////////////////////////////////////////////////
    // Prepare input channels
    ////////////////////////////////////////////////////////////////////////////

    // Channel to collect all software versions
    ch_all_versions = Channel.empty()

    // Input channels for reference genome
    ch_reference_raw = params.reference ? Channel.fromPath(params.reference, checkIfExists: true) : Channel.empty()
    // Unzip the gzipped reference genome if it is gzipped
    if (params.reference.endsWith('.gz')) {
        ch_reference = ch_reference_raw
            .map { file -> tuple([id: file.name.replaceAll(/\.gz$/, '')], file) }
            .collect()
        GUNZIP(ch_reference)
        ch_reference = GUNZIP.out.unzip_fasta
    } else {
        ch_reference = ch_reference_raw
            .map { file -> tuple([id: file.name], file) }
            .collect()
    }

    // Input channel for competitive reference genome
    ch_competitive_reference = params.competitive_reference ? Channel.fromPath( params.competitive_reference, checkIfExists: true )
        .map { it -> [[id:it.Name], it] }.collect() : Channel.empty()

    // Warn if the reference genome or competitive reference genome is larger than 20GB
    def warnIfLarge = { Path file, String label ->
        if (file.size() > 20L * 1024 * 1024 * 1024) {
            log.warn """
            ${label} '${file.name}' is larger than 20GB. This might take a long time to process.
            Consider increasing resources for indexing and mapping in the `config/<cluster_name>.config` file. However, Pipeline will continue with the current settings.
            """
        }
    }
    ch_reference.subscribe { tuple -> warnIfLarge(tuple[1], "Reference genome")}
    ch_competitive_reference.subscribe { tuple -> warnIfLarge(tuple[1], "Competitive reference genome")}

    ////////////////////////////////////////////////////////////////////////////
    // Input check and Fastq processing
    ////////////////////////////////////////////////////////////////////////////

    INPUT_CHECK ( params.samplesheet )
    ch_all_versions = ch_all_versions.mix(INPUT_CHECK.out.versions)

    if ( params.raw_fastq_qc.toBoolean() ) {
        RAW_FASTQ_QC ( INPUT_CHECK.out.reads )
        ch_all_versions = ch_all_versions.mix(RAW_FASTQ_QC.out.versions)
    }

    if ( params.fastq_processing.toBoolean() ) {
        FASTQ_PROCESSING (
            INPUT_CHECK.out.reads,
            params.kraken2_db ? file(params.kraken2_db, checkIfExists: true ) : [],
        )
        ch_all_versions = ch_all_versions.mix(FASTQ_PROCESSING.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Processed Fastq QC - Run FastQC, MultiQC and read statistics on processed reads
    ////////////////////////////////////////////////////////////////////////////

    if ( params.processed_fastq_qc.toBoolean() ) {
        PROCESSED_FASTQ_QC (
            FASTQ_PROCESSING.out.processed_reads,
            FASTQ_PROCESSING.out.json
        )
        ch_all_versions = ch_all_versions.mix(PROCESSED_FASTQ_QC.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Mapping/Competitive mapping
    ////////////////////////////////////////////////////////////////////////////

    if ( params.mapping.toBoolean() ) {
        // Competitive mapping to a concatenated reference (target plus decoy)
        if (params.competitive_reference && file( params.competitive_reference ).exists()) {
            COMPETITIVE_MAPPING (
                    ch_competitive_reference,
                    ch_reference,
                    FASTQ_PROCESSING.out.processed_reads
            )
            ch_all_versions = ch_all_versions.mix(COMPETITIVE_MAPPING.out.versions)
        // Map to the reference genome assembly
        } else {
            MAPPING (
                ch_reference,
                FASTQ_PROCESSING.out.processed_reads
            )
            ch_all_versions = ch_all_versions.mix(MAPPING.out.versions)
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Microbial screening
    ////////////////////////////////////////////////////////////////////////////

    if ( params.microbial_screening.toBoolean() ) {
        ch_ms_reference = Channel.fromPath ( params.ms_reference, checkIfExists: true )
            .map { it -> [[id:it.Name], it] }.collect()

        MICROBIAL_SCREENING (
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bam : MAPPING.out.raw_bam,
            ch_ms_reference
        )
        ch_all_versions = ch_all_versions.mix(MICROBIAL_SCREENING.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // MIA - Mapping Iterative Assembler
    ////////////////////////////////////////////////////////////////////////////

    if ( params.iterative_assembly.toBoolean() ) {
        ch_mt_reference = Channel.fromPath( params.mtDNA_reference, checkIfExists: true )
                .map { it -> [[id:it.Name], it] }.collect()

        ITERATIVE_ASSEMBLY (FASTQ_PROCESSING.out.processed_reads, ch_mt_reference)
        ch_all_versions = ch_all_versions.mix(ITERATIVE_ASSEMBLY.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Taxonomic Classification
    ////////////////////////////////////////////////////////////////////////////

    if ( params.taxonomic_classification.toBoolean() ) {
        ch_reference_database = Channel.fromPath( params.tc_reference_database, checkIfExists: true )
            .map { it -> [[id:it.Name], it] }.collect()

        TAXONOMIC_CLASSIFICATION (ch_reference_database, FASTQ_PROCESSING.out.processed_reads, workflow_name)
        ch_all_versions = ch_all_versions.mix(TAXONOMIC_CLASSIFICATION.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Repeat and CpG site identification
    ////////////////////////////////////////////////////////////////////////////

    if ( params.repeat_cpg_identification.toBoolean() ) {
        REPEAT_CPG_IDENTIFICATION ( ch_reference )
        ch_all_versions = ch_all_versions.mix(REPEAT_CPG_IDENTIFICATION.out.versions)
    }

    // Create a channel from repeat masked bed file
    ch_intervals = params.intervals ?
        Channel.fromPath(params.intervals, checkIfExists: true)
            .map { it -> [[id: it.name], it] }.collect() :
        (params.repeat_cpg_identification.toBoolean() ?
            REPEAT_CPG_IDENTIFICATION.out.repma_bed :
            Channel.value([[id: 'null'], file('null')]) // Provide null file
        )

    ////////////////////////////////////////////////////////////////////////////
    // Raw BAM QC - Run samtools flagstat, MapDamage2, AMBER and MultiQC on raw bam files
    ////////////////////////////////////////////////////////////////////////////

    if ( params.raw_bam_qc.toBoolean() ) {
        RAW_BAM_QC (
            params.competitive_reference ? ch_competitive_reference : ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bam : MAPPING.out.raw_bam,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bai : MAPPING.out.raw_bai
        )
        ch_all_versions = ch_all_versions.mix(RAW_BAM_QC.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // BAM processing - Merge bam files and remove duplicates
    ////////////////////////////////////////////////////////////////////////////

    if ( params.bam_processing.toBoolean() ) {
        BAM_PROCESSING (
            params.competitive_reference ? ch_competitive_reference : ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bam : MAPPING.out.raw_bam,
            RAW_BAM_QC.out.amber_txt,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.target_fai : MAPPING.out.fai
        )
        ch_all_versions = ch_all_versions.mix(BAM_PROCESSING.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Processed BAM QC - Run flagstat and MultiQC on processed bam files
    ////////////////////////////////////////////////////////////////////////////

    if ( params.processed_bam_qc.toBoolean() ) {
        PROCESSED_BAM_QC (
            params.competitive_reference ? ch_competitive_reference : ch_reference,
            BAM_PROCESSING.out.mq_filtered_bam,
            BAM_PROCESSING.out.mq_filtered_index,
            BAM_PROCESSING.out.rm_short_reads_bam,
            BAM_PROCESSING.out.rm_short_reads_index,
            BAM_PROCESSING.out.merged_bam_lib,
            BAM_PROCESSING.out.merged_bam_lib_index,
            BAM_PROCESSING.out.dedup_lib,
            BAM_PROCESSING.out.dedup_lib_index,
            BAM_PROCESSING.out.merged_bam_sample,
            BAM_PROCESSING.out.merged_bam_sample_index,
            BAM_PROCESSING.out.dedup_sample,
            BAM_PROCESSING.out.dedup_sample_index,
            ch_intervals
        )
        ch_all_versions = ch_all_versions.mix(PROCESSED_BAM_QC.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Random sampling of BAM files with ANGSD -doHaploCall 1
    ////////////////////////////////////////////////////////////////////////////

    if ( params.random_sampling_bam.toBoolean() ) {
        RANDOM_SAMPLING_BAM (
            BAM_PROCESSING.out.dedup_sample,
            ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.target_fai : MAPPING.out.fai
        )
        ch_all_versions = ch_all_versions.mix(RANDOM_SAMPLING_BAM.out.versions)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Variant calling with ANGSD and bcftools
    ////////////////////////////////////////////////////////////////////////////

    if ( params.variant_calling.toBoolean() ) {
        // Collecting processed BAM files for all samples
        ch_all_dedup_samples = BAM_PROCESSING.out.dedup_sample
            .map { meta, bam -> tuple([id: workflow_name], bam)}
            .groupTuple()

        // Variant calling with BCFTOOLS
        if ( params.variant_calling_bcftools.toBoolean() ) {
            VARIANT_CALLING_BCFTOOLS (
                ch_all_dedup_samples,
                ch_reference,
                params.competitive_reference ? COMPETITIVE_MAPPING.out.target_fai : MAPPING.out.fai
            )
            ch_all_versions = ch_all_versions.mix(VARIANT_CALLING_BCFTOOLS.out.versions)
        }
        // Variant calling with ANGSD
        if ( params.variant_calling_angsd.toBoolean() ) {
            VARIANT_CALLING_ANGSD (
                ch_all_dedup_samples,
                ch_reference,
                params.competitive_reference ? COMPETITIVE_MAPPING.out.target_fai : MAPPING.out.fai
            )
            ch_all_versions = ch_all_versions.mix(VARIANT_CALLING_ANGSD.out.versions)
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Output stats
    ////////////////////////////////////////////////////////////////////////////

    if ( params.stats_output.toBoolean() ) {
        STATS_OUTPUT (
            INPUT_CHECK.out.reads,
            FASTQ_PROCESSING.out.fastp_log,
            RAW_BAM_QC.out.flagstat,
            PROCESSED_BAM_QC.out.mq_filtered_bam_flagstat,
            PROCESSED_BAM_QC.out.dedup_lib_flagstat,
            BAM_PROCESSING.out.dedup_lib,
            BAM_PROCESSING.out.dedup_lib_index,
            PROCESSED_BAM_QC.out.dedup_sample_flagstat,
            BAM_PROCESSING.out.dedup_sample,
            BAM_PROCESSING.out.dedup_sample_index,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.decoy_flagstat : Channel.empty(),
            PROCESSED_BAM_QC.out.dpstats
        )
        ch_all_versions = ch_all_versions.mix(STATS_OUTPUT.out.versions)
    }

    //////////////////////////////////////////////////////////////////////////////

    // output software versions
    ch_all_versions.collectFile(name: "versions.yml", storeDir: "${params.outdir}")

}

workflow.onComplete {
    if( workflow.success ){
        log.info("""
        Thank you for using DNAharvester.

        Results are located in the folder: $params.outdir
        """)
    } else {
        log.info("""
        The pipeline completed unsuccessfully.

        Please read the error message. If you need help to solve your issue,
        feel free to reach out via slack or by opening an issue at
        https://github.com/NBISweden/DNAharvester/issues.
        """)
    }
}