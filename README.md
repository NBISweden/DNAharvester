This is the Git repository for the Bioinformatics long-term support project L_Dalen_2302:
"A novel pipeline for processing and analyzing extremely degraded DNA"

- The pipeline code will be available in the folder "read-harvester"

- "main" branch: tested and reviewed code
- "dev" branch: code in development
- feature branches: for developing and testing code for specific analyses or features


The workflow uses https://github.com/NBISweden/Earth-Biogenome-Project-pilot nextflow workflow and other nf-core pipelines as example and structured it this way:

- main.nf - script that is executed and that connects all subworkflows, depending which workflow step is chosen in configs/analysis.config
- subworkflows/ - contains one subworkflow per workflow step that is executed in main.nf . Each subworkflow runs several modules. Current subworkflows: input_check  merge_filter_reads mapping data_qc (more below about each of them)
- modules/ - contains one module per process (a process is equivalent to a Snakemake rule) that are connected in the subworkflows.


Currently the pipeline does the following:
- input_check subworkflow:
- Check the samplesheet and read in the short-read data
- merge_filter_reads subworkflow:
- runs FastP, currently only implemented for paired-end reads, with read merging, --correction and automatic adapter detection (--detect_adapter_for_pe ) and with the parameter -l for minimum read length accessible for adjustment in the analysis.config file. Do we want a hard lower limit for -l ? Or is it enough that a default value is provided in the config file?

mapping subworkflow:
- bwa index for the reference genome
- bwa aln with -l 16500 -n 0.01 -o 2
- bwa samse, taking the readgroup info from the samplesheet

data_qc subworkflow:
- FastQC on raw reads
- FastQC on trimmed / merged reads
- MapDamage2 without rescaling base qualities
- MultiQC (currently doesn’t recognize MapDamage2 output, working on that)

Next is to implement AMBER here
