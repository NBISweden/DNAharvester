# Read harvester

A novel pipeline for processing and analyzing extremely degraded DNA

## How to run read harvester

- Clone this repository to a directory on Rackham
- Create the conda environment in `environment.yaml` in a 
different directory

```
module load conda
export CONDA_ENVS_PATH=/path/to/your/directory/with/conda/environments
mamba env create -p /path/to/your/directory/with/conda/environments/read-harvester -f environment.yaml
```

- Open a tmux or screen session on Rackham, e.g. 

```
tmux new-session -s rh
```

- Activate the conda environment in the tmux session

```
module load conda
export CONDA_ENVS_PATH=/path/to/your/directory/with/conda/environments
conda activate read-harvester
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
`configs/test.config`.

- Add your UPPMAX compute project ID to the file `nextflow.config` 
in line 115. (not tested yet)

- Start the pipeline in the tmux session with the activated 
conda environment as follows

