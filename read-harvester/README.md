# Read harvester

A novel pipeline for processing and analyzing extremely degraded DNA

## How to run read harvester

### Rackham (HPC cluster)

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
file (see `configs/test-rackham.config`).

- Read harvester runs AMBER, a tool that is not available 
as container image or conda package yet. To run the pipeline, 
clone the AMBER Github repository (`https://github.com/tvandervalk/AMBER/`) 
into a different location and copy the file `AMBER` to 
`read-harvester/bin`. 

- Open a tmux or screen session on Rackham, e.g.: 

```
tmux new-session -s rh
```

- Activate the conda environment in the tmux session: 

```
module load conda
export CONDA_ENVS_PATH=/proj/sllstore2017093/b2016342/b2016342_nobackup/lts/conda_environments
conda activate read-harvester
```

- Start the pipeline in the tmux session with the activated 
conda environment as follows (replace `config/test-rackham.config` 
with your custom config file): 

```
nextflow run -c configs/test-rackham.config -profile uppmax main.nf &> YYMMDD_rh.out
```

> The `uppmax` profile is set up to submit each process as 
a job to the slurm queue, and to use singularity (apptainer) 
to run each process in a container with the required software. 

### Locally (small test dataset)

- Clone this repository to a directory on your computer

- Create the conda environment from `environment.yaml`, 
if you haven't done so yet (with mamba or conda): 

```
mamba env create -f environment.yaml
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
`configs/test-local.config`.

- Read harvester runs AMBER, a tool that is not available 
as container image or conda package yet. To run the pipeline, 
clone the AMBER Github repository (`https://github.com/tvandervalk/AMBER/`) 
into a different location and copy the file `AMBER` to 
`read-harvester/bin`. 

- Open a tmux or screen session, e.g.: 

```
tmux new-session -s rh
```

- Activate the conda environment in the tmux session: 

```
conda activate read-harvester
```

- Start the pipeline in the tmux session with the activated 
conda environment as follows (replace `config/test-local.config` 
with your custom config file): 

```
nextflow run -c configs/test-local.config -profile conda main.nf &> YYMMDD_rh.out
```

> The `conda` profile is set up to to use conda to run each 
process in a conda environment with the required software. 
