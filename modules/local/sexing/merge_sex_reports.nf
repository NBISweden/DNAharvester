process MERGE_SEX_REPORTS {
    label 'process_low'

    conda "conda-forge::python=3.11 conda-forge::pandas=2.0"
    container "${ (workflow.containerEngine == 'apptainer' || workflow.containerEngine == 'singularity') && !task.ext.apptainer_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.0.3' :
        'quay.io/biocontainers/pandas:2.0.3' }"

    input:
    path(sex_reports)

    output:
    path("sex_determination_summary.tsv"), emit: summary
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import glob
    import sys

    # Read all individual reports
    reports = glob.glob("*.sex_determination.tsv")
    dfs = [pd.read_csv(f, sep="\\t") for f in reports]

    # Merge all reports
    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.sort_values("sample_id")

    # Write summary
    merged.to_csv("sex_determination_summary.tsv", sep="\\t", index=False)

    # Versions
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
