process CREATE_REPEATS_BED {
    tag "$fasta"
    label 'process_low'

    conda "conda-forge::python=3.13.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/python:3.13.0--a8086dc1de1c4e39' :
        'community.wave.seqera.io/library/python:3.13.0--a025ad9838d75455' }"

    input:
    tuple val(meta2), path(fasta)
    tuple val(meta2), path(repeatmasker_out)

    output:
    tuple val(meta2), path("*_repeats.bed")     , emit: repeats_bed
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta2.id}"
    def args = task.ext.args ?: ''

    """
    python3 -c '
from sys import argv
input_rep_out = argv[1]
output_rep_bed = argv[2]
with open(input_rep_out, "r") as f, open(output_rep_bed, "w") as b:
    for line in f:
        edited = line.strip().split()
        if len(edited) > 6 and edited[0].isdigit():
            start = str(int(edited[5]) - 1)
            b.write(edited[4] + "\\t" + start + "\\t" + edited[6] + "\\n")
' ${repeatmasker_out} ${prefix}_repeats.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(echo \$(python --version 2>&1) | awk '{print \$2}')
    END_VERSIONS
    """
}
