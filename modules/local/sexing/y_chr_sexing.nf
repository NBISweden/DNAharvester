process Y_CHR_SEXING {
    tag "$meta.id"
    label 'process_y_chr_sexing'

    conda "conda-forge::pandas=2.3.3 conda-forge::matplotlib=3.10.1"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_pandas:8965e9dfe0245807' :
        'community.wave.seqera.io/library/matplotlib_pandas:76f3c63ec67531f0' }"

    input:
    tuple val(meta), path(idxstats)
    val(y_chromosome)
    val(x_chromosome)

    output:
    tuple val(meta), path("${meta.id}.y_chr_sexing.tsv")    , emit: sexing_report
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    # Parameters
    sample_id = "${meta.id}"
    y_chromosome = "${y_chromosome}"
    x_chromosome = "${x_chromosome}"

    # Read idxstats file
    # Format: chr_name, chr_length, mapped_reads, unmapped_reads
    df = pd.read_csv("${idxstats}", sep="\\t", header=None,
                     names=["chr", "length", "mapped", "unmapped"])

    # Get Y chromosome reads
    y_reads = df[df["chr"] == y_chromosome]["mapped"].sum()

    # Get X chromosome reads
    x_reads = df[df["chr"] == x_chromosome]["mapped"].sum()

    # Calculate Y / (X + Y) ratio
    total_sex_chr_reads = x_reads + y_reads
    y_ratio = y_reads / total_sex_chr_reads if total_sex_chr_reads > 0 else 0

    # Create output
    result = {
        "sample_id": sample_id,
        "y_over_xy": round(y_ratio, 6)
    }

    result_df = pd.DataFrame([result])
    result_df.to_csv("${meta.id}.y_chr_sexing.tsv", sep="\\t", index=False)

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
