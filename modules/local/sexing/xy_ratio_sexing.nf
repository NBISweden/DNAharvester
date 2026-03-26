process XY_RATIO_SEXING {
    tag "${workflow_name}"
    label 'process_xy_ratio_sexing'

    conda "conda-forge::pandas=2.3.3 conda-forge::matplotlib=3.10.1"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'oras://community.wave.seqera.io/library/matplotlib_pandas:8965e9dfe0245807' :
        'community.wave.seqera.io/library/matplotlib_pandas:76f3c63ec67531f0' }"

    input:
    val(sample_ids)
    path(idxstats_files)
    val(x_chromosome)
    val(y_chromosome)
    val(autosomes)
    val(workflow_name)

    output:
    path("${workflow_name}_xy_ratio_sexing_plot.pdf")    , emit: sexing_plot
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Convert lists to comma-separated strings
    def sample_ids_str = sample_ids.join(',')
    def idxstats_files_str = idxstats_files.collect { it.name }.join(',')
    def auto_chrs = autosomes ? autosomes.tokenize(' ').collect { "'${it}'" }.join(', ') : ""

    """
    #!/usr/bin/env python3

    import pandas as pd
    import matplotlib.pyplot as plt
    import sys

    # Parameters
    sample_ids = "${sample_ids_str}".split(',')
    idxstats_files = "${idxstats_files_str}".split(',')
    x_chromosome = "${x_chromosome}"
    y_chromosome = "${y_chromosome}"
    autosomes = [${auto_chrs}] if "${autosomes}" else []

    results = []

    for sample_id, idxstats_file in zip(sample_ids, idxstats_files):
        # Read idxstats file
        # Format: chr_name, chr_length, mapped_reads, unmapped_reads
        df = pd.read_csv(idxstats_file, sep="\\t", header=None,
                         names=["chr", "length", "mapped", "unmapped"])

        # Calculate autosomal baseline (Na) - sum of reads aligning to autosomes
        auto_reads = df[df["chr"].isin(autosomes)]["mapped"].sum()

        # Get X chromosome reads
        x_reads = df[df["chr"] == x_chromosome]["mapped"].sum()

        # Get Y chromosome reads
        y_reads = df[df["chr"] == y_chromosome]["mapped"].sum()

        # Calculate Rx and Ry ratios
        # Rx = X reads / autosomal reads
        # Ry = Y reads / autosomal reads
        rx = x_reads / auto_reads if auto_reads > 0 else 0
        ry = y_reads / auto_reads if auto_reads > 0 else 0

        # Sex determination based on Rx and Ry
        # Females: high Rx (~0.04-0.05 depending on genome), low/zero Ry
        # Males: lower Rx (~0.02-0.025), measurable Ry
        if ry < 0.001 and rx > 0.03:
            sex_call = "Female"
        elif ry > 0.001:
            sex_call = "Male"
        else:
            sex_call = "Undetermined"

        results.append({
            "sample_id": sample_id,
            "x_chr_reads": int(x_reads),
            "y_chr_reads": int(y_reads),
            "autosome_reads": int(auto_reads),
            "Rx": round(rx, 6),
            "Ry": round(ry, 6),
            "sex_call": sex_call
        })

    # Create results dataframe
    results_df = pd.DataFrame(results)


    # Create the Rx vs Ry plot
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot each sample as a dot
    scatter = ax.scatter(results_df["Rx"], results_df["Ry"],
                         s=100, c='steelblue', edgecolors='black', linewidths=0.5,
                         alpha=0.8, zorder=3)

    # Add sample labels next to each dot
    for idx, row in results_df.iterrows():
        ax.annotate(row["sample_id"],
                    (row["Rx"], row["Ry"]),
                    xytext=(5, 5), textcoords='offset points',
                    fontsize=8, alpha=0.8)

    # Labels and title
    ax.set_xlabel("Rx (X chromosome reads / Autosomal reads)", fontsize=11)
    ax.set_ylabel("Ry (Y chromosome reads / Autosomal reads)", fontsize=11)
    ax.set_title("Sex Determination: Rx vs Ry Ratio", fontsize=12, fontweight='bold')

    # Set axis limits with padding (15% margin on each side)
    rx_range = results_df["Rx"].max() - results_df["Rx"].min()
    ry_range = results_df["Ry"].max() - results_df["Ry"].min()
    rx_padding = max(rx_range * 0.15, results_df["Rx"].max() * 0.15)
    ry_padding = max(ry_range * 0.15, results_df["Ry"].max() * 0.15)
    ax.set_xlim(left=-rx_padding, right=results_df["Rx"].max() + rx_padding)
    ax.set_ylim(bottom=-ry_padding, top=results_df["Ry"].max() + ry_padding)

    # Add grid for better readability
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)

    # Add legend with sex indicators
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='none', edgecolor='none', label='High Rx, Low Ry → Female'),
        Patch(facecolor='none', edgecolor='none', label='Lower Rx, Ry > 0 → Male')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(f"${workflow_name}_xy_ratio_sexing_plot.pdf", dpi=300, bbox_inches='tight')
    plt.close()

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
        f.write(f"    matplotlib: {plt.matplotlib.__version__}\\n")
    """
}
