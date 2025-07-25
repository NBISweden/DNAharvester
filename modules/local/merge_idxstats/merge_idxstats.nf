process MERGE_IDXSTATS {
    tag "${meta.id}"
    label 'process_medium'

    conda "conda-forge::gawk=5.3.1"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.1' }"


    input:
    tuple val(meta), path(input_files, stageAs: "?/*")

    output:
    tuple val(meta), path("${meta.id}.idxstats.txt")    , emit: merged_idxstats
    path "versions.yml"                                 , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | head -n 1 | awk '{print \$1, \$2, \$3}')
    END_VERSIONS
    """
}
