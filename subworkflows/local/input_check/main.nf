#! /usr/bin/env nextflow

// Check input samplesheet and get read channels
// Modified from nf-core pipeline template
// https://github.com/nf-core/tools/blob/e5ce6ce20304835bd40f102f038b7e1aadc888b2/nf_core/pipeline-template/subworkflows/local/input_check.nf

include { SAMPLESHEET_CHECK      } from '../../../modules/local/samplesheet_check'

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

    // If species_identification is set to true, check if si_reference_database is provided
    if (params.species_identification.toBoolean() && !params.si_reference_database) {
        log.error """`species_identification` is set to true, but `si_reference_database` is not provided. Please provide the path to the species identification reference database.
        Exiting the pipeline......!
        """
        System.exit(1)
    }

    // Dont allow readlength set to auto for species_identification analysis
    if (params.species_identification.toBoolean() && params.readlength == 'auto') {
        log.error """For species_identification analysis, `readlength` should not be set to 'auto' as the mapping is done against multiple reference genomes.
        Please set `readlength` to a specific value.
        Exiting the pipeline......!
        """
        System.exit(1)
    }

    // if iterative_assembly is set to true, check if mtDNA_reference is provided
    if (params.iterative_assembly.toBoolean() && !params.mtDNA_reference) {
        log.error """`iterative_assembly` is set to true, but `mtDNA_reference` is not provided. Please provide the path to the mitochondrial reference genome from any closely related species.
        Exiting the pipeline......!
        """
        System.exit(1)
    }

    // Make sure either repeat_cpg_identification is set to true or intervals (repeats_masked.bed) is provided. both cannot be true
    if (params.repeat_cpg_identification.toBoolean() && params.intervals) {
        log.error """Both `repeat_cpg_identification` is set to true and `intervals` (file: ${params.intervals}) is provided. Please select only one option.
        Either provide path to BED file with reference genome positions to include in the downstream analysis OR
        set `repeat_cpg_identification` to identify repeats and use the generated repeats_masked.bed file for downstream analysis.
        Exiting the pipeline......!
        """
        System.exit(1)
    }

    // Check if microbial_screening is set to true and reads are being filtered with kraken2 before mapping
    if (params.microbial_screening.toBoolean() && params.kraken2.toBoolean()) {
        log.warn """Both `microbial_screening` and `kraken2` are set to true.
        Please note that if reads are filtered with kraken2 before mapping, the microbial screening step
        will includes unmapped reads and also reads filtered out by kraken2.
        """
    }

    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

    emit:
    reads                                                                                   // channel: [ val(meta), [ reads ] ]
    csv                                      = SAMPLESHEET_CHECK.out.csv                    // channel: [ samplesheet.valid.csv ]
    versions                                 = ch_versions                                  // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map: [id, single_end]
    // TO DO Update to take strandedness and damage treatment
    def meta = [:]
    meta.id                = row.sample_id + "_" + row.library_id + "_" + row.lane
    meta.library_type      = row.library_type
    meta.single_end        = row.single_end.toBoolean()
    // readgroup: ID = ID = readgroup id (flowcell-id.lane-nr.library-index-nr), SM = sample-id, PL = sequencing platform (e.g. Illumina, NovaSeq), LB = library-index-nr based on the number of unique libraries per sample in samplesheet
    meta.read_group        = "@RG\\tID:" + row.flowcell_id + "." + row.lane + "." + row.library_id + "\\tSM:" + row.sample_id + "\\tPL:" + row.seq_platform + "\\tLB:" + row.library_id

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