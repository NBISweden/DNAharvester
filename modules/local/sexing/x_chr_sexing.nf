process X_CHR_SEXING {
    tag "$meta.id"
    label 'process_x_chr_sexing'

    conda "conda-forge::pandas=2.3.3 conda-forge::matplotlib=3.10.1"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_pandas:8965e9dfe0245807' :
        'community.wave.seqera.io/library/matplotlib_pandas:76f3c63ec67531f0' }"

    input:
    tuple val(meta), path(idxstats)
    val (x_chromosome)
    val (autosomes)

    output:
    tuple val(meta), path("${meta.id}.x_chr_sexing.tsv")            , emit: sexing_report
    tuple val(meta), path("${meta.id}.x_chr_sexing_summary.tsv")    , emit: sexing_summary
    tuple val(meta), path("${meta.id}.x_chr_ploidy.pdf")            , emit: sexing_ploidy_plot
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def auto_chrs = autosomes ? autosomes.tokenize(' ').collect { "'${it}'" }.join(', ') : ""

    """
    #!/usr/bin/env python3

    import pandas as pd
    import matplotlib.pyplot as plt
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

    # Create ploidy plot
    # Calculate expected ploidy: normalize each chr to autosomal average (diploid = 2x)
    auto_mean_normalized = auto_normalized

    # Add ploidy column to output_df
    plot_df = output_df.copy()
    plot_df["ploidy"] = (plot_df["normalized_mapped_reads"] / auto_mean_normalized) * 2 if auto_mean_normalized > 0 else 0

    # Separate autosomes and X chromosome
    auto_plot = plot_df[plot_df["chr_type"] == "autosome"]
    x_plot = plot_df[plot_df["chr_type"] == "X"]

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot autosomes as black dots
    ax.scatter(range(len(auto_plot)), auto_plot["ploidy"],
               c='black', s=80, marker='o', label='Autosome', zorder=3)

    # Plot X chromosome as red X marker
    x_position = len(auto_plot)
    if len(x_plot) > 0:
        ax.scatter([x_position], x_plot["ploidy"],
                   c='red', s=150, marker='x', linewidths=3, label='X chromosome', zorder=3)

    # Add horizontal lines at highest and lowest autosome ploidy values
    if len(auto_plot) > 0:
        auto_max = auto_plot["ploidy"].max()
        auto_min = auto_plot["ploidy"].min()
        ax.axhline(y=auto_max, color='gray', linestyle='--', alpha=0.5)
        ax.axhline(y=auto_min, color='gray', linestyle='--', alpha=0.5)

    # X-axis labels
    all_chrs = list(auto_plot["chr"]) + list(x_plot["chr"])
    ax.set_xticks(range(len(all_chrs)))
    ax.set_xticklabels(all_chrs, rotation=45, ha='right', fontsize=8)

    # Labels and title
    ax.set_xlabel("Chromosome", fontsize=10)
    ax.set_ylabel("Relative ploidy", fontsize=10)
    ax.set_title(f"Chromosome ploidy - {sample_id}", fontsize=12, fontweight='bold')

    # Set y-axis limits and integer ticks only
    y_max = max(3, int(plot_df["ploidy"].max()) + 1)
    ax.set_ylim(0, y_max)
    ax.set_yticks(range(0, y_max + 1))

    # Legend
    ax.legend(loc='upper right', fontsize=9)

    plt.tight_layout()
    plt.savefig("${meta.id}.x_chr_ploidy.pdf", dpi=300, bbox_inches='tight')
    plt.close()

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
        f.write(f"    matplotlib: {plt.matplotlib.__version__}\\n")
    """
}
