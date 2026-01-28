#! /usr/bin/env nextflow

include { BWA_INDEX                 as CM_BWA_INDEX                     } from '../../../modules/local/bwa/bwa_index.nf'
include { BOWTIE2_BUILD             as CM_BOWTIE2_BUILD                 } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BWA_ALN                   as CM_BWA_ALN                       } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_MEM                   as CM_BWA_MEM                       } from '../../../modules/local/bwa/bwa_mem.nf'
include { SPLIT_FASTQ               as CM_SPLIT_FASTQ                   } from '../../../modules/local/awk/split_fastq.nf'
include { BWA_ALN                   as CM_BWA_ALN_SHORT                 } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_MEM                   as CM_BWA_MEM_LONG                  } from '../../../modules/local/bwa/bwa_mem.nf'
include { SAMTOOLS_MERGE            as CM_BWA_ALN_MEM_MERGE             } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_MERGE            as CM_MERGED_UNMERGED_READS_BAM     } from '../../../modules/local/samtools/samtools_merge.nf'
include { BOWTIE2                   as CM_BOWTIE2                       } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_INDEX            as CM_SAMTOOLS_INDEX_RAW            } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX            as CM_SAMTOOLS_FAIDX                } from '../../../modules/local/samtools/samtools_faidx.nf'
include { FAI_TO_BED                as CM_FAI_TO_BED                    } from '../../../modules/local/fai2bed/main'
include { SAMTOOLS_FAIDX            as CM_SAMTOOLS_FAIDX_TARGET         } from '../../../modules/local/samtools/samtools_faidx.nf'
include { FAI_TO_BED                as CM_FAI_TO_BED_TARGET             } from '../../../modules/local/fai2bed/main'
include { BEDTOOLS_SUBTRACT         as CM_BEDTOOLS_SUBTRACT_TARGET      } from '../../../modules/local/bedtools/subtract/main'
include { SAMTOOLS_VIEW_REGIONS     as CM_SAMTOOLS_VIEW_DECOY           } from '../../../modules/local/samtools/samtools_view_regions.nf'
include { SAMTOOLS_INDEX            as CM_SAMTOOLS_INDEX_DECOY          } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT         as CM_FLAGSTAT_DECOY                } from '../../../modules/local/samtools/samtools_flagstat.nf'
include { MULTIQC                   as CM_MULTIQC_DECOY                 } from '../../../modules/nf-core/multiqc/main'
include { SAMTOOLS_VIEW_REGIONS     as CM_SAMTOOLS_VIEW_TARGET          } from '../../../modules/local/samtools/samtools_view_regions.nf'
include { SAMTOOLS_INDEX            as CM_SAMTOOLS_INDEX_TARGET         } from '../../../modules/nf-core/samtools/index/main'

workflow COMPETITIVE_MAPPING {
    take:
    competitive_reference
    reference
    reads

    main:
    ch_versions = Channel.empty()

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1. Index the competitive reference genome if it is not already indexed
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Build the BWA index - only if BWA is selected as mapping tool
    def bwa_tools = ['bwa-aln', 'bwa-mem', 'bwa-aln-mem']
    if (bwa_tools.contains( params.mapping_tool_ancient ) || bwa_tools.contains( params.mapping_tool_modern )) {
        CM_BWA_INDEX (competitive_reference, file(params.reference).getParent())
        ch_versions                 = ch_versions.mix(CM_BWA_INDEX.out.versions)
        ch_competitive_bwa_index    = CM_BWA_INDEX.out.index_dir
    }
    // Build the Bowtie2 index - only if Bowtie2 is selected as mapping tool
    def bowtie2_tools = ['bowtie2']
    if (bowtie2_tools.contains( params.mapping_tool_ancient ) || bowtie2_tools.contains( params.mapping_tool_modern )) {
        CM_BOWTIE2_BUILD (competitive_reference, file(params.reference).getParent())
        ch_versions                  = ch_versions.mix(CM_BOWTIE2_BUILD.out.versions)
        ch_competitive_bowtie2_index = CM_BOWTIE2_BUILD.out.index_dir
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 2. Mapping merged reads or single-end reads
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ch_raw_bam = Channel.empty()

    // Branch reads into ancient and modern based on sample_type metadata
    reads.branch {
        ancient: it[0].sample_type == 'ancient'
        modern : it[0].sample_type == 'modern'
    }.set { ch_reads_branched }

    // Route reads to specific tools
    ch_reads_bwa_aln        = Channel.empty()
    ch_reads_bwa_mem        = Channel.empty()
    ch_reads_bwa_aln_mem    = Channel.empty()
    ch_reads_bowtie2        = Channel.empty()

    // Configure routing for ancient samples
    if (params.mapping_tool_ancient == 'bwa-aln')           { ch_reads_bwa_aln = ch_reads_bwa_aln.mix(ch_reads_branched.ancient) }
    else if (params.mapping_tool_ancient == 'bwa-mem')      { ch_reads_bwa_mem = ch_reads_bwa_mem.mix(ch_reads_branched.ancient) }
    else if (params.mapping_tool_ancient == 'bwa-aln-mem')  { ch_reads_bwa_aln_mem = ch_reads_bwa_aln_mem.mix(ch_reads_branched.ancient) }
    else if (params.mapping_tool_ancient == 'bowtie2')      { ch_reads_bowtie2 = ch_reads_bowtie2.mix(ch_reads_branched.ancient) }

    // Configure routing for modern samples
    if (params.mapping_tool_modern == 'bwa-aln')            { ch_reads_bwa_aln = ch_reads_bwa_aln.mix(ch_reads_branched.modern) }
    else if (params.mapping_tool_modern == 'bwa-mem')       { ch_reads_bwa_mem = ch_reads_bwa_mem.mix(ch_reads_branched.modern) }
    else if (params.mapping_tool_modern == 'bwa-aln-mem')   { ch_reads_bwa_aln_mem = ch_reads_bwa_aln_mem.mix(ch_reads_branched.modern) }
    else if (params.mapping_tool_modern == 'bowtie2')       { ch_reads_bowtie2 = ch_reads_bowtie2.mix(ch_reads_branched.modern) }


    // Run Mapping Tools
    // BWA ALN
    if (params.mapping_tool_ancient == 'bwa-aln' || params.mapping_tool_modern == 'bwa-aln') {
        CM_BWA_ALN ( ch_reads_bwa_aln, ch_competitive_bwa_index )
        ch_versions         = ch_versions.mix(CM_BWA_ALN.out.versions)
        ch_raw_bam          = ch_raw_bam.mix(CM_BWA_ALN.out.bam)
    }
    // BWA MEM
    if (params.mapping_tool_ancient == 'bwa-mem' || params.mapping_tool_modern == 'bwa-mem') {
        CM_BWA_MEM ( ch_reads_bwa_mem, ch_competitive_bwa_index )
        ch_versions         = ch_versions.mix(CM_BWA_MEM.out.versions)
        ch_raw_bam          = ch_raw_bam.mix(CM_BWA_MEM.out.bam)
    }
    // BWA ALN-MEM
    if (params.mapping_tool_ancient == 'bwa-aln-mem' || params.mapping_tool_modern == 'bwa-aln-mem') {
        //split fastq
        CM_SPLIT_FASTQ ( ch_reads_bwa_aln_mem )
        ch_versions         = ch_versions.mix(CM_SPLIT_FASTQ.out.versions)
        //align short reads with BWA ALN
        CM_BWA_ALN_SHORT ( CM_SPLIT_FASTQ.out.short_reads, ch_competitive_bwa_index )
        ch_versions         = ch_versions.mix(CM_BWA_ALN_SHORT.out.versions)
        //align long reads with BWA MEM
        CM_BWA_MEM_LONG ( CM_SPLIT_FASTQ.out.long_reads,  ch_competitive_bwa_index )
        ch_versions         = ch_versions.mix(CM_BWA_MEM_LONG.out.versions)
        //merge BAMs from short and long reads
        ch_raw_bam_aln_mem      = CM_BWA_ALN_SHORT.out.bam.join(CM_BWA_MEM_LONG.out.bam)
                .map { meta, file1, file2 -> [meta, [file1, file2]] }
        CM_BWA_ALN_MEM_MERGE ( ch_raw_bam_aln_mem, reference )
        ch_versions         = ch_versions.mix(CM_BWA_ALN_MEM_MERGE.out.versions)
        ch_raw_bam          = ch_raw_bam.mix(CM_BWA_ALN_MEM_MERGE.out.bam)
    }
    // BOWTIE2
    if (params.mapping_tool_ancient == 'bowtie2' || params.mapping_tool_modern == 'bowtie2') {
        CM_BOWTIE2 ( ch_reads_bowtie2, ch_competitive_bowtie2_index )
        ch_versions         = ch_versions.mix(CM_BOWTIE2.out.versions)
        ch_raw_bam          = ch_raw_bam.mix(CM_BOWTIE2.out.bam)
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 3. Merge the mapped unmerged reads if provided
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    def ch_merged_raw_bam = null
    if (params.merge_reads.toBoolean() && params.keep_unmerged_reads.toBoolean()) {
        // Group the unmerged reads BAMs by sample ID (removing the '-unmerged' suffix)
        ch_raw_bam_grouped = ch_raw_bam.map { meta, bam ->
                def new_meta = meta.clone()
                new_meta.id = meta.id.replace("-unmerged", "")
                tuple(new_meta.id, new_meta, bam)
            }
            .groupTuple()
            .map { id, metas, bams ->
                // pick one meta; choose the one with single_end == false if present
                def final_meta = metas.find { !it.single_end } ?: metas[0]
                // remove grouping id, return meta + joined bam list
                tuple(final_meta, bams)
            }

        CM_MERGED_UNMERGED_READS_BAM ( ch_raw_bam_grouped, reference )
        ch_versions = ch_versions.mix(CM_MERGED_UNMERGED_READS_BAM.out.versions)
        ch_merged_raw_bam = CM_MERGED_UNMERGED_READS_BAM.out.bam
    }

    // Use merged raw bam if created, else use original raw bam
    ch_raw_bam_for_processing = ch_merged_raw_bam ?: ch_raw_bam

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 4. Processing BAM files
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Index the BAM file
    CM_SAMTOOLS_INDEX_RAW ( ch_raw_bam_for_processing )
    ch_versions                      = ch_versions.mix(CM_SAMTOOLS_INDEX_RAW.out.versions)

    // Split the BAM file into target genome and decoy genome
    ch_concatenated_bam_index        = ch_raw_bam_for_processing.join( CM_SAMTOOLS_INDEX_RAW.out.bai )

    // Generate *.fai index for the concatenated reference
    CM_SAMTOOLS_FAIDX ( competitive_reference )
    ch_versions                      = ch_versions.mix(CM_SAMTOOLS_FAIDX.out.versions)

    // Convert to *.bed format
    CM_FAI_TO_BED ( CM_SAMTOOLS_FAIDX.out.fai )

    // Generate *.fai index for the target reference genome
    CM_SAMTOOLS_FAIDX_TARGET ( reference )
    ch_versions                      = ch_versions.mix(CM_SAMTOOLS_FAIDX_TARGET.out.versions)

    // Convert to *.bed format
    CM_FAI_TO_BED_TARGET ( CM_SAMTOOLS_FAIDX_TARGET.out.fai )

    // Decoy genome
    // Extract the decoy genome chromosomes from the concatenated genome BED file
    CM_BEDTOOLS_SUBTRACT_TARGET ( CM_FAI_TO_BED.out.bed, CM_FAI_TO_BED_TARGET.out.bed )

    // Extract the region from the BAM file
    CM_SAMTOOLS_VIEW_DECOY ( ch_concatenated_bam_index, CM_BEDTOOLS_SUBTRACT_TARGET.out.bed )
    ch_versions                      = ch_versions.mix(CM_SAMTOOLS_VIEW_DECOY.out.versions)
     // Index the BAM file
    CM_SAMTOOLS_INDEX_DECOY ( CM_SAMTOOLS_VIEW_DECOY.out.bam )
    ch_versions                      = ch_versions.mix(CM_SAMTOOLS_INDEX_DECOY.out.versions)

    // Run samtools flagstat and MultiQC
    ch_flagstat_decoy                = CM_SAMTOOLS_VIEW_DECOY.out.bam.join(CM_SAMTOOLS_INDEX_DECOY.out.bai)
    CM_FLAGSTAT_DECOY ( ch_flagstat_decoy )
    ch_versions                      = ch_versions.mix(CM_FLAGSTAT_DECOY.out.versions)

    ch_multiqc_decoy_files           = CM_FLAGSTAT_DECOY.out.flagstat.map{ meta, flagstat -> flagstat }.collect()
    ch_multiqc_config                = params.multiqc_config       ? Channel.fromPath( params.multiqc_config,       checkIfExists: true ) : Channel.empty()
    ch_multiqc_extra_config          = params.multiqc_extra_config ? Channel.fromPath( params.multiqc_extra_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                  = params.multiqc_logo         ? Channel.fromPath( params.multiqc_logo,         checkIfExists: true ) : Channel.empty()

    CM_MULTIQC_DECOY (
        ch_multiqc_decoy_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_extra_config.toList(),
        ch_multiqc_logo.toList()
    )
    ch_versions                      = ch_versions.mix(CM_MULTIQC_DECOY.out.versions)

    // Target genome
    // Extract the region from the BAM file
    CM_SAMTOOLS_VIEW_TARGET ( ch_concatenated_bam_index, CM_FAI_TO_BED_TARGET.out.bed )
    ch_versions                      = ch_versions.mix(CM_SAMTOOLS_VIEW_TARGET.out.versions)
    // Index the BAM file containing only the target genome
    CM_SAMTOOLS_INDEX_TARGET ( CM_SAMTOOLS_VIEW_TARGET.out.bam )
    ch_versions                      = ch_versions.mix(CM_SAMTOOLS_INDEX_TARGET.out.versions)

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    emit:
    competitive_fai                  = CM_SAMTOOLS_FAIDX.out.fai                   // channel: path(index)
    target_fai                       = CM_SAMTOOLS_FAIDX_TARGET.out.fai            // channel: path(index)
    multiqc_decoy_report             = CM_MULTIQC_DECOY.out.report.toList()        // channel: [ val(meta), path(report) ]
    decoy_bam                        = CM_SAMTOOLS_VIEW_DECOY.out.bam              // channel: [ val(meta), [ bam ] ]
    decoy_bai                        = CM_SAMTOOLS_INDEX_DECOY.out.bai             // channel: [ val(meta), [ bai ] ]
    decoy_flagstat                   = CM_FLAGSTAT_DECOY.out.flagstat              // channel: [ val(meta), [ flagstat ] ]
    target_bam                       = CM_SAMTOOLS_VIEW_TARGET.out.bam             // channel: [ val(meta), [ bam ] ]
    target_bai                       = CM_SAMTOOLS_INDEX_TARGET.out.bai            // channel: [ val(meta), [ bai ] ]
    versions                         = ch_versions                                 // channel: [ versions.yml ]
}