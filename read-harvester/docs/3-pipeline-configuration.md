# Pipeline configuration

## Sample sheet

Create a sample sheet, listing metadata information for 
each sample, including the location of raw Fastq files with 
the sequencing data on the same compute system or links to 
publicly available read data in a repository, and place it 
in `assets/`. Here are examples: `assets/samplesheet.csv` or 
`assets/samplesheet_testdata-eager-mammoth.csv`. Note that 
the pipeline currently only supports paired-end data. 

## Configuration file

Next, create a custom pipeline configuration file with paths 
to the reference genome (on the same compute system or a 
link to a publicly available sequence) and the sample sheet, 
pipeline steps to be run, a path to a results directory, and 
tool-specific parameters. `assets/custom.config` is a template, 
a filled-out example is available here: 
`read-harvester/assets/test-dardel-eager_verena.config`. 

If you are running readharvester on a cluster using slurm, 
such as Dardel from PDC/KTH, add your slurm compute project 
ID to the parameter `process.clusterOptions` in the custom 
config file (see `read-harvester/assets/test-dardel-eager_verena.config` 
for an example). 

## Compute resources

Default compute resources for processes run by readharvester 
are specified in `configs/nf-core-defaults.config`. Compute 
resources for the PDC/KTH cluster Dardel on which the pipeline 
was developed are specified in `configs/dardel.config`. Note 
that if a process fails due to limited memory or time, changes 
to labels in the `configs/*.config` files will be applied to 
all processes with this label. If this is not wanted, the 
process label of a specific process can be changed directly 
in its `*.nf` file under `modules/local/` or `modules/nf-core/` 
and therein. 