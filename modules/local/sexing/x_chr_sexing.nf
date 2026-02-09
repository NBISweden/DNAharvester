process X_CHR_SEXING {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.3' :
        'quay.io/biocontainers/pandas:2.0.3' }"

    input:
    tuple val(meta), path(idxstats)
    val (x_chromosome)
    val (autosomes)

    output:
    tuple val(meta), path("${meta.id}.x_chr_sexing.tsv"), emit: sex_report
    tuple val(meta), path("${meta.id}.x_chr_sexing_summary.tsv"), emit: sex_summary
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def auto_chrs = autosomes ? autosomes.tokenize(' ').collect { "'${it}'" }.join(', ') : ""

    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    # Parameters
    sample_id = "${meta.id}"
    x_chromosome = "${x_chromosome}"
    autosomes = [${auto_chrs}] if "${autosomes}" else []

    # Read idxstats file
    # Format: chr_name, chr_length, mapped_reads, unmapped_reads
    df = pd.read_csv("${idxstats}", sep="\\t", header=None,
                     names=["chr", "length", "mapped", "unmapped"])

    # Filter for autosomes and X chromosome only
    chromosomes_of_interest = autosomes + [x_chromosome]
    df_filtered = df[df["chr"].isin(chromosomes_of_interest)].copy()

    # Calculate X chromosome statistics
    x_df = df[df["chr"] == x_chromosome]
    x_reads = x_df["mapped"].sum()
    x_length = x_df["length"].sum()
    x_normalized = x_reads / x_length if x_length > 0 else 0

    # Calculate per-chromosome statistics
    df_filtered["chr_size"] = df_filtered["length"]
    df_filtered["read_mapped"] = df_filtered["mapped"]
    df_filtered["normalized_mapped_reads"] = df_filtered["mapped"] / df_filtered["length"]

    # Calculate X/autosome ratio for each chromosome
    df_filtered["X/autosome"] = x_normalized / df_filtered["normalized_mapped_reads"]

    # Mark chromosome type
    df_filtered["chr_type"] = df_filtered["chr"].apply(
        lambda c: "X" if c == x_chromosome else "autosome"
    )

    # Select and order columns for detailed output
    output_df = df_filtered[["chr", "chr_size", "read_mapped", "normalized_mapped_reads", "X/autosome", "chr_type"]]

    # Write detailed per-chromosome output
    output_df.to_csv("${meta.id}.x_chr_sexing.tsv", sep="\\t", index=False)

    # Calculate summary statistics
    auto_df = df[df["chr"].isin(autosomes)]
    auto_reads = auto_df["mapped"].sum()
    auto_length = auto_df["length"].sum()
    auto_normalized = auto_reads / auto_length if auto_length > 0 else 0

    # X-to-autosome ratio (normalized)
    x_auto_ratio = x_normalized / auto_normalized if auto_normalized > 0 else 0

    # Sex determination based on ratio
    # Females (XX) typically have ratio ~1.0, Males (XY) typically have ratio ~0.5
    if x_auto_ratio >= 0.75:
        sex_call = "Female"
    elif x_auto_ratio <= 0.6:
        sex_call = "Male"
    else:
        sex_call = "Undetermined"

    # Create summary
    summary = {
        "sample_id": sample_id,
        "total_mapped_reads": int(df["mapped"].sum()),
        "x_chr_reads": int(x_reads),
        "x_chr_length": int(x_length),
        "x_normalized_reads": round(x_normalized, 9),
        "autosome_reads": int(auto_reads),
        "autosome_length": int(auto_length),
        "autosome_normalized_reads": round(auto_normalized, 9),
        "x_to_autosome_ratio": round(x_auto_ratio, 6),
        "sex_call": sex_call
    }

    summary_df = pd.DataFrame([summary])
    summary_df.to_csv("${meta.id}.x_chr_sexing_summary.tsv", sep="\\t", index=False)

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
