process SEX_DETERMINATION {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.3' :
        'quay.io/biocontainers/pandas:2.0.3' }"

    input:
    tuple val(meta), path(idxstats)
    val(x_chromosomes)
    val(autosomes)
    val(y_chromosomes)

    output:
    tuple val(meta), path("${meta.id}.sex_determination.tsv"), emit: sex_report
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def x_chrs = x_chromosomes ? x_chromosomes.tokenize(' ').collect { "'${it}'" }.join(', ') : ""
    def auto_chrs = autosomes ? autosomes.tokenize(' ').collect { "'${it}'" }.join(', ') : ""
    def y_chrs = y_chromosomes ? y_chromosomes.tokenize(' ').collect { "'${it}'" }.join(', ') : ""

    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    # Parameters
    sample_id = "${meta.id}"
    x_chromosomes = [${x_chrs}] if "${x_chromosomes}" else []
    autosomes = [${auto_chrs}] if "${autosomes}" else []
    y_chromosomes = [${y_chrs}] if "${y_chromosomes}" else []

    # Determine which methods are enabled
    x_ratio_method = bool(x_chromosomes and autosomes)
    y_chr_method = bool(y_chromosomes)

    # Read idxstats file
    # Format: chr_name, chr_length, mapped_reads, unmapped_reads
    df = pd.read_csv("${idxstats}", sep="\\t", header=None,
                     names=["chr", "length", "mapped", "unmapped"])

    # Calculate read counts
    x_reads = df[df["chr"].isin(x_chromosomes)]["mapped"].sum() if x_chromosomes else 0
    auto_reads = df[df["chr"].isin(autosomes)]["mapped"].sum() if autosomes else 0
    y_reads = df[df["chr"].isin(y_chromosomes)]["mapped"].sum() if y_chromosomes else 0

    # Calculate lengths for normalization
    x_length = df[df["chr"].isin(x_chromosomes)]["length"].sum() if x_chromosomes else 0
    auto_length = df[df["chr"].isin(autosomes)]["length"].sum() if autosomes else 0
    y_length = df[df["chr"].isin(y_chromosomes)]["length"].sum() if y_chromosomes else 0

    # Total reads used in analysis
    total_reads_analyzed = x_reads + auto_reads + y_reads
    total_mapped_reads = df["mapped"].sum()

    # Initialize results
    results = {
        "sample_id": sample_id,
        "total_mapped_reads": int(total_mapped_reads),
        "reads_used_for_analysis": int(total_reads_analyzed),
        "x_chr_reads": int(x_reads) if x_ratio_method else "NA",
        "x_chr_length": int(x_length) if x_ratio_method else "NA",
        "autosome_reads": int(auto_reads) if x_ratio_method else "NA",
        "autosome_length": int(auto_length) if x_ratio_method else "NA",
        "y_chr_reads": int(y_reads) if y_chr_method else "NA",
        "y_chr_length": int(y_length) if y_chr_method else "NA",
        "x_to_auto_ratio_raw": "NA",
        "x_to_auto_ratio_normalized": "NA",
        "y_to_auto_ratio_normalized": "NA",
        "x_ratio_method_enabled": x_ratio_method,
        "y_chr_method_enabled": y_chr_method
    }

    # X-to-autosome ratio method
    if x_ratio_method and auto_reads > 0:
        # Raw ratio (not length normalized)
        raw_ratio = x_reads / auto_reads if auto_reads > 0 else 0
        results["x_to_auto_ratio_raw"] = round(raw_ratio, 6)

        # Normalized ratio (per base)
        if x_length > 0 and auto_length > 0:
            x_normalized = x_reads / x_length
            auto_normalized = auto_reads / auto_length
            norm_ratio = x_normalized / auto_normalized if auto_normalized > 0 else 0
            results["x_to_auto_ratio_normalized"] = round(norm_ratio, 6)

    # Y-chromosome method
    if y_chr_method and auto_length > 0 and y_length > 0:
        y_normalized = y_reads / y_length
        auto_normalized = auto_reads / auto_length if auto_length > 0 else 0
        y_ratio = y_normalized / auto_normalized if auto_normalized > 0 else 0
        results["y_to_auto_ratio_normalized"] = round(y_ratio, 6)

    # Write output
    output_df = pd.DataFrame([results])
    output_df.to_csv("${meta.id}.sex_determination.tsv", sep="\\t", index=False)

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
