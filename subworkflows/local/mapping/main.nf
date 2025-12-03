#! /usr/bin/env nextflow

include { SAMTOOLS_FAIDX            } from '../../../modules/local/samtools/samtools_faidx.nf'
include { BWA_INDEX                 } from '../../../modules/local/bwa/bwa_index.nf'
include { BOWTIE2_BUILD             } from '../../../modules/local/bowtie2/bowtie2_build.nf'
include { BWA_ALN                   } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_ALN as BWA_ALN_R1     } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_ALN as BWA_ALN_R2     } from '../../../modules/local/bwa/bwa_aln.nf'
include { BWA_SAMSE                 } from '../../../modules/local/bwa/bwa_samse.nf'
include { BWA_SAMPE                 } from '../../../modules/local/bwa/bwa_sampe.nf'
include { BWA_MEM                   } from '../../../modules/local/bwa/bwa_mem.nf'
include { SPLIT_FASTQ               } from '../../../modules/local/awk/split_fastq.nf'


include { BOWTIE2                   } from '../../../modules/local/bowtie2/bowtie2.nf'
include { SAMTOOLS_MERGE            } from '../../../modules/local/samtools/samtools_merge.nf'
include { SAMTOOLS_INDEX            } from '../../../modules/nf-core/samtools/index/main'


workflow MAPPING {
    take:
    reference
    reads // merged paired-end reads or trimmed single-end reads
    unmerged_reads // unmerged paired-end reads

    main:
    ch_versions         = Channel.empty()

    SAMTOOLS_FAIDX ( reference )
    ch_versions         = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ////////////////////////////////////////////////////////////////////////////
    // Index the reference genome if it is not already indexed
    ////////////////////////////////////////////////////////////////////////////

    // Build the BWA index - only if BWA is selected as mapping tool
    def bwa_tools = ['bwa-aln', 'bwa-mem', 'bwa-aln-mem']
    if (bwa_tools.contains( params.mapping_tool_ancient ) ||
        bwa_tools.contains( params.mapping_tool_modern )) {
        BWA_INDEX (reference, file(params.reference).getParent())
        ch_versions         = ch_versions.mix(BWA_INDEX.out.versions)
        ch_bwa_index        = BWA_INDEX.out.index_dir
    }

    // Build the Bowtie2 index - only if Bowtie2 is selected as mapping tool
    def bowtie2_tools = ['bowtie2']
    if (bowtie2_tools.contains( params.mapping_tool_ancient ) ||
        bowtie2_tools.contains( params.mapping_tool_modern )) {
        BOWTIE2_BUILD (reference, file(params.reference).getParent())
        ch_versions         = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        ch_bowtie2_index    = BOWTIE2_BUILD.out.index_dir
    }

    ////////////////////////////////////////////////////////////////////////////
    // Mapping merged reads or single-end reads
    ////////////////////////////////////////////////////////////////////////////

    ch_bam = Channel.empty()

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
    if (params.mapping_tool_ancient == 'bwa-aln') ch_reads_bwa_aln = ch_reads_bwa_aln.mix(ch_reads_branched.ancient)
    if (params.mapping_tool_ancient == 'bwa-mem') ch_reads_bwa_mem = ch_reads_bwa_mem.mix(ch_reads_branched.ancient)
    if (params.mapping_tool_ancient == 'bwa-aln-mem') ch_reads_bwa_aln_mem = ch_reads_bwa_aln_mem.mix(ch_reads_branched.ancient)
    if (params.mapping_tool_ancient == 'bowtie2') ch_reads_bowtie2 = ch_reads_bowtie2.mix(ch_reads_branched.ancient)

    // Configure routing for modern samples
    if (params.mapping_tool_modern == 'bwa-aln') ch_reads_bwa_aln = ch_reads_bwa_aln.mix(ch_reads_branched.modern)
    if (params.mapping_tool_modern == 'bwa-mem') ch_reads_bwa_mem = ch_reads_bwa_mem.mix(ch_reads_branched.modern)
    if (params.mapping_tool_modern == 'bwa-aln-mem') ch_reads_bwa_aln_mem = ch_reads_bwa_aln_mem.mix(ch_reads_branched.modern)
    if (params.mapping_tool_modern == 'bowtie2') ch_reads_bowtie2 = ch_reads_bowtie2.mix(ch_reads_branched.modern)


    // Run Mapping Tools
    // BWA ALN
    if (params.mapping_tool_ancient == 'bwa-aln' || params.mapping_tool_modern == 'bwa-aln') {
        BWA_ALN ( ch_reads_bwa_aln, ch_bwa_index )
        ch_versions         = ch_versions.mix(BWA_ALN.out.versions)
        ch_bwa_samse        = ch_reads_bwa_aln.join(BWA_ALN.out.sai)
        BWA_SAMSE ( ch_bwa_samse, ch_bwa_index )
        ch_versions         = ch_versions.mix(BWA_SAMSE.out.versions)
        ch_bam              = ch_bam.mix(BWA_SAMSE.out.bam)
    }
    // BWA MEM
    if (params.mapping_tool_ancient == 'bwa-mem' || params.mapping_tool_modern == 'bwa-mem') {
        BWA_MEM ( ch_reads_bwa_mem, ch_bwa_index )
        ch_versions         = ch_versions.mix(BWA_MEM.out.versions)
        ch_bam              = ch_bam.mix(BWA_MEM.out.bam)
    }
    // BWA ALN-MEM
    if (params.mapping_tool_ancient == 'bwa-aln-mem' || params.mapping_tool_modern == 'bwa-aln-mem') {
        //split fastq
        SPLIT_FASTQ ( ch_reads_bwa_aln_mem )
        ch_versions         = ch_versions.mix(SPLIT_FASTQ.out.versions)
        //align short reads with BWA ALN
        BWA_ALN ( SPLIT_FASTQ.out.short_reads, ch_bwa_index )
        ch_versions         = ch_versions.mix(BWA_ALN.out.versions)
        ch_bwa_samse_short  = SPLIT_FASTQ.out.short_reads.join(BWA_ALN.out.sai)
        BWA_SAMSE ( ch_bwa_samse_short, ch_bwa_index )
        ch_versions         = ch_versions.mix(BWA_SAMSE.out.versions)
        //align long reads with BWA MEM
        BWA_MEM ( SPLIT_FASTQ.out.long_reads,  ch_bwa_index )
        ch_versions         = ch_versions.mix(BWA_MEM.out.versions)
        //merge BAMs from short and long reads
        ch_bam_aln_mem      = BWA_SAMSE.out.bam.join(BWA_MEM.out.bam)
                .map { meta, file1, file2 -> [meta, [file1, file2]] }
        SAMTOOLS_MERGE ( ch_bam_aln_mem, reference )
        ch_versions         = ch_versions.mix(SAMTOOLS_MERGE.out.versions)
        ch_bam              = ch_bam.mix(SAMTOOLS_MERGE.out.bam)
    }
    // BOWTIE2
    if (params.mapping_tool_ancient == 'bowtie2' || params.mapping_tool_modern == 'bowtie2') {
        BOWTIE2 ( ch_reads_bowtie2, ch_bowtie2_index )
        ch_versions         = ch_versions.mix(BOWTIE2.out.versions)
        ch_bam              = ch_bam.mix(BOWTIE2.out.bam)
    }

    ////////////////////////////////////////////////////////////////////////////
    // Mapping unmerged reads if params.keep_unmerged is true
    ////////////////////////////////////////////////////////////////////////////

    // processed unmerged reads if provided
    if (params.keep_unmerged.toBoolean()) {
        // Split channels into paired-end and single-end reads
        ch_paired_end_reads = ch_bam.filter { meta, file -> !meta.single_end }
        ch_single_end_reads = ch_bam.filter { meta, file -> meta.single_end }

        if (ch_paired_end_reads) {
            // Prepare channels for unmerged reads. Add R1 and R2 labels to the meta.id. This is needed to avoid file name clashes in the BWA_ALN process output
            ch_R1 = unmerged_reads.map { meta, file1, file2 ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}-R1"
                [new_meta, file1]
            }
            ch_R2 = unmerged_reads.map { meta, file1, file2 ->
                def new_meta = meta.clone()
                new_meta.id = "${meta.id}-R2"
                [new_meta, file2]
            }

            // Align unmerged reads with BWA ALN and then run BWA SAMPE
            if (params.mapping_tool == 'bwa-aln') {
                BWA_ALN_R1 ( ch_R1, ch_reference_index)
                ch_versions             = ch_versions.mix(BWA_ALN_R1.out.versions)
                BWA_ALN_R2 ( ch_R2, ch_reference_index)
                ch_versions             = ch_versions.mix(BWA_ALN_R2.out.versions)
                // Fix meta.id to original id without -R1 or -R2 suffix after alignment
                ch_sai_R1 = BWA_ALN_R1.out.sai.map { meta, file ->
                    meta.id = meta.id.replaceAll(/-R1$/, '')
                    [meta, file]
                }
                ch_sai_R2 = BWA_ALN_R2.out.sai.map { meta, file ->
                    meta.id = meta.id.replaceAll(/-R2$/, '')
                    [meta, file]
                }
                // Join the aligned unmerged reads with their respective SAI files
                ch_bwa_sampe_unmerged = unmerged_reads.join(ch_sai_R1).join(ch_sai_R2)
                // Run BWA SAMPE
                BWA_SAMPE ( ch_bwa_sampe_unmerged, ch_reference_index )
                ch_versions             = ch_versions.mix(BWA_SAMPE.out.versions)
                ch_bam_unmerged         = BWA_SAMPE.out.bam
            } else if (params.mapping_tool == 'bowtie2') {
                // Align unmerged reads with Bowtie2
                BOWTIE2 ( unmerged_reads, ch_reference_index )
                ch_versions             = ch_versions.mix(BOWTIE2.out.versions)
                ch_bam_unmerged         = BOWTIE2.out.bam
            }

            // merge merged and unmerged bam files
            ch_bam_merged_unmerged = ch_bam.join(ch_bam_unmerged)
                    .map { meta, file1, file2 -> [meta, [file1, file2]] }
            SAMTOOLS_MERGE ( ch_bam_merged_unmerged, reference )
            ch_versions             = ch_versions.mix(SAMTOOLS_MERGE.out.versions)

            // Final BAM channel with both merged and unmerged reads
            ch_bam = SAMTOOLS_MERGE.out.bam.mix(ch_single_end_reads)

        }
    }

    ////////////////////////////////////////////////////////////////////////////

    // Index the BAM file
    SAMTOOLS_INDEX ( ch_bam )
    ch_versions             = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    emit:
    fai                     = SAMTOOLS_FAIDX.out.fai             // channel: path(index)
    // index                   = ch_reference_index                 // channel: path(index)
    bam                     = ch_bam                             // channel: [ val(meta), [ bam ] ]
    bai                     = SAMTOOLS_INDEX.out.bai             // channel: [ val(meta), [ bai ] ]
    versions                = ch_versions                        // channel: [ versions.yml ]
}