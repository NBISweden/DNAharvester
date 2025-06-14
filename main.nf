#! /usr/bin/env nextflow

// Analysis script for DNAharvester pipeline

nextflow.enable.dsl = 2

// Import subworkflows
include { GUNZIP                     } from "$projectDir/modules/local/gunzip/main"
include { INPUT_CHECK                } from "$projectDir/subworkflows/local/input_check/main"
include { FASTQ_PROCESSING           } from "$projectDir/subworkflows/local/fastq_processing/main"
include { PROCESSED_FASTQ_QC         } from "$projectDir/subworkflows/local/processed_fastq_qc/main"
include { MAPPING                    } from "$projectDir/subworkflows/local/mapping/main"
include { COMPETITIVE_MAPPING        } from "$projectDir/subworkflows/local/competitive_mapping/main"
include { REPEAT_CPG_IDENTIFICATION  } from "$projectDir/subworkflows/local/repeat_cpg_identification/main"
include { RAW_BAM_QC                 } from "$projectDir/subworkflows/local/raw_bam_qc/main"
include { BAM_PROCESSING             } from "$projectDir/subworkflows/local/bam_processing/main"
include { PROCESSED_BAM_QC           } from "$projectDir/subworkflows/local/processed_bam_qc/main"
include { RANDOM_SAMPLING_BAM        } from "$projectDir/subworkflows/local/random_sampling_bam/main"
include { STATS_OUTPUT               } from "$projectDir/subworkflows/local/stats_output/main"


workflow {

    // Define workflow stages
    def recognized_workflow_stages = ['fastq_processing', 'mapping', 'repeat_cpg_identification', 'processed_fastq_qc', 'raw_bam_qc', 'bam_processing', 'processed_bam_qc', 'random_sampling_bam', 'stats_output']

    // Check input
    def workflow_steps = params.steps.tokenize(",")
    if ( ! workflow_steps.every { it in recognized_workflow_stages } ) {
        error "Unrecognised workflow step in $params.steps ( $recognized_workflow_stages )"
    }

    // The primary workflow for the DNAharvester pipeline
    log.info("""
    Running DNAharvester.
    """)

    ch_all_versions = Channel.empty()

    // Input channels for reference genome
    ch_reference_raw = Channel.fromPath(params.reference, checkIfExists: true)
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
            Consider increasing the resources or pre-indexing it with BWA index. However, Pipeline will continue with the current settings.
            """
        }
    }
    ch_reference.subscribe { tuple -> warnIfLarge(tuple[1], "Reference genome")}
    ch_competitive_reference.subscribe { tuple -> warnIfLarge(tuple[1], "Competitive reference genome")}

    // Input check, Merge paired-end reads, trim adapters and filter for minimum read length
    if ( 'fastq_processing' in workflow_steps ) {
        INPUT_CHECK ( params.samplesheet )
        ch_all_versions = ch_all_versions.mix(INPUT_CHECK.out.versions)

        FASTQ_PROCESSING (
            params.kraken2_db ? file(params.kraken2_db, checkIfExists: true ) : [],
            INPUT_CHECK.out.reads
        )
        ch_all_versions = ch_all_versions.mix(FASTQ_PROCESSING.out.versions)
    }

    // Run FastQC, MultiQC and read statistics on processed reads
    if ( 'processed_fastq_qc' in workflow_steps ) {
        PROCESSED_FASTQ_QC (
            FASTQ_PROCESSING.out.reads,
            FASTQ_PROCESSING.out.json
        )
        ch_all_versions = ch_all_versions.mix(PROCESSED_FASTQ_QC.out.versions)
    }

    // Map with bwa-aln (aDNA parameters) and convert to bam
    if ( 'mapping' in workflow_steps ) {
        // Competitive mapping to a concatenated reference (target plus decoy)
        if (params.competitive_reference && file( params.competitive_reference ).exists()) {
            COMPETITIVE_MAPPING (
                    ch_competitive_reference,
                    ch_reference,
                    FASTQ_PROCESSING.out.reads
            )
            ch_all_versions = ch_all_versions.mix(COMPETITIVE_MAPPING.out.versions)
        // Map to the reference genome assembly
        } else {
            MAPPING (
                ch_reference,
                FASTQ_PROCESSING.out.reads
            )
            ch_all_versions = ch_all_versions.mix(MAPPING.out.versions)
        }
    }

    // Run RepeatModeler and RepeatMasker to identify repeats and a custom script to identify CpG sites
    if ( 'repeat_cpg_identification' in workflow_steps ) {
        REPEAT_CPG_IDENTIFICATION ( ch_reference )
        ch_all_versions = ch_all_versions.mix(REPEAT_CPG_IDENTIFICATION.out.versions)
    }

    // Create a channel from repeat masked bed file
    ch_intervals = params.intervals ? Channel.fromPath(params.intervals, checkIfExists: true)
        .map { it -> [[id: it.name], it] }.collect()
        : ('repeat_cpg_identification' in workflow_steps ? REPEAT_CPG_IDENTIFICATION.out.repma_bed : Channel.empty())

    // Run samtools flagstat, MapDamage2, AMBER and MultiQC on raw bam files
    if ( 'raw_bam_qc' in workflow_steps ) {
        RAW_BAM_QC (
            params.competitive_reference ? ch_competitive_reference : ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bam : MAPPING.out.bam,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bai : MAPPING.out.bai,
        )
        ch_all_versions = ch_all_versions.mix(RAW_BAM_QC.out.versions)
    }
    // Merge bam files per index, remove duplicates, merge bam files per sample, remove duplicates
    if ( 'bam_processing' in workflow_steps ) {
        BAM_PROCESSING (
            ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.target_fai : MAPPING.out.fai,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bam : MAPPING.out.bam,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.bai : MAPPING.out.bai,
            RAW_BAM_QC.out.amber_txt
        )
        ch_all_versions = ch_all_versions.mix(BAM_PROCESSING.out.versions)
    }

    // Run flagstat and MultiQC on processed bam files
    if ( 'processed_bam_qc' in workflow_steps ) {
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

    // Run ANGSD -doHaploCall 1 to sample a random base at each site from bam files
    if ( 'random_sampling_bam' in workflow_steps ) {
        RANDOM_SAMPLING_BAM (
            BAM_PROCESSING.out.dedup_sample,
            ch_reference,
            params.competitive_reference ? COMPETITIVE_MAPPING.out.target_fai : MAPPING.out.fai
        )
        ch_all_versions = ch_all_versions.mix(RANDOM_SAMPLING_BAM.out.versions)
    }

    // Output stats
    if ( 'stats_output' in workflow_steps ) {
        STATS_OUTPUT (
            INPUT_CHECK.out.reads,
            FASTQ_PROCESSING.out.fastp_log,
            RAW_BAM_QC.out.flagstat,
            PROCESSED_BAM_QC.out.mq_filtered_bam_flagstat,
            PROCESSED_BAM_QC.out.dedup_lib_flagstat,
            BAM_PROCESSING.out.dedup_lib
        )
        ch_all_versions = ch_all_versions.mix(STATS_OUTPUT.out.versions)
    }

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