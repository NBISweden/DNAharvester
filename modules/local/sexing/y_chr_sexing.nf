process Y_CHR_SEXING {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.3' :
        'quay.io/biocontainers/pandas:2.0.3' }"

    input:
    tuple val(meta), path(idxstats)
    val(y_chromosomes)
    val(autosomes)

    output:
    tuple val(meta), path("${meta.id}.y_chr_sexing.tsv"), emit: sex_report
    tuple val(meta), path("${meta.id}.y_chr_sexing_summary.tsv"), emit: sex_summary
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def y_chrs = y_chromosomes ? y_chromosomes.tokenize(' ').collect { "'${it}'" }.join(', ') : ""
    def auto_chrs = autosomes ? autosomes.tokenize(' ').collect { "'${it}'" }.join(', ') : ""

    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    # Parameters
    sample_id = "${meta.id}"
    y_chromosomes = [${y_chrs}] if "${y_chromosomes}" else []
    autosomes = [${auto_chrs}] if "${autosomes}" else []

    # Read idxstats file
    # Format: chr_name, chr_length, mapped_reads, unmapped_reads
    df = pd.read_csv("${idxstats}", sep="\\t", header=None,
                     names=["chr", "length", "mapped", "unmapped"])

    # Filter for autosomes and Y chromosomes only
    chromosomes_of_interest = autosomes + y_chromosomes
    df_filtered = df[df["chr"].isin(chromosomes_of_interest)].copy()

    # Calculate Y chromosome statistics
    y_df = df[df["chr"].isin(y_chromosomes)]
    y_reads = y_df["mapped"].sum()
    y_length = y_df["length"].sum()
    y_normalized = y_reads / y_length if y_length > 0 else 0

    # Calculate per-chromosome statistics
    df_filtered["chr_size"] = df_filtered["length"]
    df_filtered["read_mapped"] = df_filtered["mapped"]
    df_filtered["normalized_mapped_reads"] = df_filtered["mapped"] / df_filtered["length"]

    # Calculate Y/autosome ratio for each chromosome
    df_filtered["Y/autosome"] = y_normalized / df_filtered["normalized_mapped_reads"]

    # Mark chromosome type
    df_filtered["chr_type"] = df_filtered["chr"].apply(
        lambda x: "Y" if x in y_chromosomes else "autosome"
    )

    # Select and order columns for detailed output
    output_df = df_filtered[["chr", "chr_size", "read_mapped", "normalized_mapped_reads", "Y/autosome", "chr_type"]]

    # Write detailed per-chromosome output
    output_df.to_csv("${meta.id}.y_chr_sexing.tsv", sep="\\t", index=False)

    # Calculate summary statistics
    auto_df = df[df["chr"].isin(autosomes)]
    auto_reads = auto_df["mapped"].sum()
    auto_length = auto_df["length"].sum()
    auto_normalized = auto_reads / auto_length if auto_length > 0 else 0

    # Y-to-autosome ratio (normalized)
    y_auto_ratio = y_normalized / auto_normalized if auto_normalized > 0 else 0

    # Sex determination based on Y chromosome presence
    # Males typically have Y reads, Females typically have very few/no Y reads
    if y_auto_ratio >= 0.3:
        sex_call = "Male"
    elif y_auto_ratio <= 0.1:
        sex_call = "Female"
    else:
        sex_call = "Undetermined"

    # Create summary
    summary = {
        "sample_id": sample_id,
        "total_mapped_reads": int(df["mapped"].sum()),
        "y_chr_reads": int(y_reads),
        "y_chr_length": int(y_length),
        "y_normalized_reads": round(y_normalized, 9),
        "autosome_reads": int(auto_reads),
        "autosome_length": int(auto_length),
        "autosome_normalized_reads": round(auto_normalized, 9),
        "y_to_autosome_ratio": round(y_auto_ratio, 6),
        "sex_call": sex_call
    }

    summary_df = pd.DataFrame([summary])
    summary_df.to_csv("${meta.id}.y_chr_sexing_summary.tsv", sep="\\t", index=False)

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
