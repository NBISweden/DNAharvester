# Read harvester

A novel pipeline for processing and analyzing extremely degraded DNA

## How to run read harvester

- Clone this repository to a directory on Rackham

- The conda environment from `environment.yaml` has 
already been created in a directory accessible for 
everyone in the storage project using the following 
code:

```
module load conda
export CONDA_ENVS_PATH=/proj/sllstore2017093/b2016342/b2016342_nobackup/lts/conda_environments
mamba env create -p /proj/sllstore2017093/b2016342/b2016342_nobackup/lts/conda_environments/read-harvester -f environment.yaml
```

- Create a sample sheet in `configs/`, listing metadata 
information for each sample, including the location of raw 
fastq files with the sequencing data. Here is an example: 
`configs/samplesheet.csv`. Note that the pipeline currently 
only supports paired-end data.

- Create a custom pipeline configuration file with paths 
to input data, pipeline steps to be run, path to results 
directory, and tool-specific parameters. `configs/custom.config` 
is a template, a filled-out example is available here: 
`configs/test-rackham.config`.

- Add your UPPMAX compute project ID to the parameter 
`process.clusterOptions` in line 54 of the custom config 
file (see `configs/test-rackham.config`)

- Open a tmux or screen session on Rackham, e.g. 

```
tmux new-session -s rh
```

- Activate the conda environment in the tmux session

```
module load conda
export CONDA_ENVS_PATH=/proj/sllstore2017093/b2016342/b2016342_nobackup/lts/conda_environments
conda activate read-harvester
```

- Start the pipeline in the tmux session with the activated 
conda environment as follows (replace `config/test-rackham.config` 
with your custom config file)

```
nextflow run -c configs/test-rackham.config -profile uppmax main.nf
```