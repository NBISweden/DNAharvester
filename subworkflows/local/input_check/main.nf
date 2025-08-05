#! /usr/bin/env nextflow

// Check input samplesheet and get read channels
// Modified from nf-core pipeline template
// https://github.com/nf-core/tools/blob/e5ce6ce20304835bd40f102f038b7e1aadc888b2/nf_core/pipeline-template/subworkflows/local/input_check.nf

include { SAMPLESHEET_CHECK      } from '../../../modules/local/samplesheet_check'
include { FASTQC as FASTQC_RAW   } from '../../../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_RAW } from '../../../modules/nf-core/multiqc/main'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    ch_versions                              = Channel.empty()

    // Check if `variant_calling` is set to true but neither `variant_calling_bcftools` nor `variant_calling_angsd` is enabled
    if (params.variant_calling.toBoolean()) {
        if (!(params.variant_calling_bcftools.toBoolean() || params.variant_calling_angsd.toBoolean())) {
            log.error """`variant_calling` is set to true, but neither `variant_calling_bcftools` nor `variant_calling_angsd` is enabled. Please set at least one of them to true.
            Exiting the pipeline......!
            """
            System.exit(1)
        }
    }

    // Check if both `mapdamage2_rescale` and `remove_transitions` are set to true
    if ( params.mapdamage2_rescale.toBoolean() && params.remove_transitions.toBoolean() ) {
        log.error """Both `mapdamage2_rescale` and `remove_transitions` are set to true. Please select only one option.
        Exiting the pipeline......!
        """
    System.exit(1)
    }

    // Dont allow readlength set to auto for mystery_sample analysis
     if (params.mystery_sample.toBoolean() && params.readlength == 'auto') {
        log.error """For mystery_sample analysis, `readlength` should not be set to 'auto' as the mapping is done against multiple reference genomes.
        Please set `readlength` to a specific value.
        Exiting the pipeline......!
        """
        System.exit(1)
    }

    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }
    ch_versions                              = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    FASTQC_RAW ( reads )
    ch_versions                              = ch_versions.mix(FASTQC_RAW.out.versions)

    // Run MultiQC on FastQC output
    ch_multiqc_processed_files                   = FASTQC_RAW.out.zip.map{ meta, qcfile -> qcfile }.collect()
    ch_multiqc_config                        = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config                  = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                          = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    MULTIQC_RAW (
        ch_multiqc_processed_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                              = ch_versions.mix(MULTIQC_RAW.out.versions)

    emit:
    reads                                                                                   // channel: [ val(meta), [ reads ] ]
    csv                                      = SAMPLESHEET_CHECK.out.csv                    // channel: [ samplesheet.valid.csv ]
    fastqc_html                              = FASTQC_RAW.out.html                          // channel: [ val(meta), path(html) ]
    fastqc_zip                               = FASTQC_RAW.out.zip                           // channel: [ val(meta), path(zip) ]
    multiqc_report                           = MULTIQC_RAW.out.report.toList()              // channel: [ val(meta), path(report) ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map: [id, single_end]
    // TO DO Update to take strandedness and damage treatment
    def meta = [:]
    meta.id                = row.sample + "_" + row.library_id + "_" + row.lane
    meta.library_type      = row.library_type
    meta.single_end        = row.single_end.toBoolean()
    // readgroup: ID = ID = readgroup id (flowcell-id.lane-nr.library-index-nr), SM = sample-id, PL = sequencing platform (e.g. Illumina, NovaSeq), LB = library-index-nr based on the number of unique libraries per sample in samplesheet
    meta.read_group        = "@RG\\tID:" + row.flowcell_id + "." + row.lane + "." + row.library_id + "\\tSM:" + row.sample + "\\tPL:" + row.seq_platform + "\\tLB:" + row.library_id

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return fastq_meta
}