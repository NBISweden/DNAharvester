process MERGE_IDXSTATS {
    tag "${meta.id}"
    label 'process_medium'

    conda "bioconda::htslib=1.21 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.apptainer_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9e/9edc2564215d5cd137a8b25ca8a311600987186d406b092022444adf3c4447f7/data' :
        'community.wave.seqera.io/library/htslib_samtools:1.21--6cb89bfd40cbaabf' }"

    input:
    tuple val(meta), path(input_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${meta.id}.idxstats.txt"), emit: merged_idxstats

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_names = input_files.collect { it.getBaseName().replace(/\.idxstats$/, "") }
    def header = (['reference_name', 'size'] + sample_names).join('\t')
    def input_str = input_files.collect { it.toString() }.join(' ')
    def nfiles = input_files.size()

    """
    paste ${input_str} | awk -v OFS='\\t' -v N=${nfiles} '
    BEGIN { print "${header}" }
    {
        out = \$1 OFS \$2
        for (i = 0; i < N; i++) {
            col = (i * 4) + 3
            out = out OFS \$col
        }
        print out
    }' > ${prefix}.idxstats.txt

    """
}
