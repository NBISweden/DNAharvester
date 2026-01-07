process ADAPTCLEAN {
    tag "$meta.id"
    label 'process_adaptclean'

    conda "bioconda::gcc"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/gcc_gxx_jq_wget:549b249d6ae2ebd5' :
        'community.wave.seqera.io/library/gcc_gxx_jq_wget:2458d207f4475602' }"

    input:
    tuple val(meta), path(processed_reads), path(json)

    output:
    tuple val(meta), path("*.adaptclean.fastq.gz"), emit: adaptclean_reads
    path "versions.yml"                           , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    adaptclean = "https://raw.githubusercontent.com/bilalbioinfo/AdaptClean/refs/heads/main/AdaptClean.cpp"


    if (meta.single_end) {
        """
        if [ ! -f ${projectDir}/bin/AdaptClean ]; then
            wget $adaptclean
            g++ -std=c++17 -O2 -Wall AdaptClean.cpp -o AdaptClean -lz
            mv AdaptClean ${projectDir}/bin/
        fi


        adapter_seq=\$(jq -r '.adapter_cutting.read1_adapter_sequence' ${json} | head -c4)
        sequencing_cycle=\$(jq -r '.read1_before_filtering.total_cycles' ${json})
        AdaptClean \\
            ${processed_reads} \\
            ${prefix}.se.adaptclean.fastq.gz \\
            \${sequencing_cycle} \\
            \${adapter_seq}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            AdaptClean: 0.1.0
        END_VERSIONS

        """
    }
    else {
        """
        if [ ! -f ${projectDir}/bin/AdaptClean ]; then
            wget $adaptclean
            g++ -std=c++17 -O2 -Wall AdaptClean.cpp -o AdaptClean -lz
            mv AdaptClean ${projectDir}/bin/
        fi

        adapter_seq_r1=\$(jq -r '.adapter_cutting.read1_adapter_sequence' ${json} | head -c4)
        adapter_seq_r2=\$(jq -r '.adapter_cutting.read2_adapter_sequence' ${json} | head -c4)
        sequencing_cycle_r1=\$(jq -r '.read1_before_filtering.total_cycles' ${json})
        sequencing_cycle_r2=\$(jq -r '.read2_before_filtering.total_cycles' ${json})

        AdaptClean \\
            ${processed_reads[0]} \\
            ${prefix}.R1.adaptclean.fastq.gz \\
            \${sequencing_cycle_r1} \\
            \${adapter_seq_r1}
        AdaptClean \\
            ${processed_reads[1]} \\
            ${prefix}.R2.adaptclean.fastq.gz \\
            \${sequencing_cycle_r2} \\
            \${adapter_seq_r2} \\

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            AdaptClean: 0.1.0
        END_VERSIONS

        """
    }
}

