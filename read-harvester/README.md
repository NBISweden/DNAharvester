# Read harvester

A novel pipeline for processing and analyzing extremely degraded DNA

## How to run read harvester

### Dardel (HPC cluster at PDC/KTH)

- Load the following modules:

```
module load PDC bioinfo-tools conda apptainer tmux
```

> Note that `tmux` is only available as a module on Dardel 
but the equivalent tool `screen` is pre-installed and does 
not need to be loaded. 

> Apptainer (former singularity) can use your `scratch` for 
caching, which is a temporary directory with unlimited space 
by adding this row to your `~/.bashrc`: 
`export NXF_SINGULARITY_CACHEDIR=$PDC_TMP`. 

- Clone this repository to a directory on Dardel

- Create the pipeline conda environment from `environment.yaml`. 
Since home directories on Dardel are limited in storage space, 
you need to create a directory in your storage project for the 
conda environment to be installed in, and run the following 
command: 

```
conda env create -f environment.yml -p /cfs/klemming/projects/supr/sllstore.../read-harvester
```

> Note that you can save storage space in your storage project 
on Dardel by creating a common pipeline conda environment for 
several people. 

- Create a sample sheet, listing metadata information for 
each sample, including the location of raw fastq files with 
the sequencing data. Here is an example: `assets/samplesheet.csv`. 
Note that the pipeline currently only supports paired-end data. 

- Create a custom pipeline configuration file with paths 
to input data, pipeline steps to be run, path to results 
directory, and tool-specific parameters. `assets/custom.config` 
is a template, a filled-out example is available here: 
`read-harvester/assets/test-dardel-eager_verena.config`. 

- Add your UPPMAX compute project ID to the parameter 
`process.clusterOptions` in line 68 of the custom config 
file (see `read-harvester/assets/test-dardel-eager_verena.config`). 

- Read harvester runs AMBER, a tool that is not available 
as container image or conda package yet. To run the pipeline, 
clone the AMBER Github repository (`https://github.com/tvandervalk/AMBER/`) 
into a different location and copy the file `AMBER` to 
`read-harvester/bin`. 

- Open a tmux or screen session on Dardel, e.g.: 

```
tmux new-session -s rh
```

- Activate the conda environment in the tmux session, replacing 
the path to the directory where you created the conda environment: 

```
export CONDA_ENVS_PATH=/cfs/klemming/projects/supr/sllstore.../conda_environments/
conda activate read-harvester
```

- Start the pipeline in the tmux session with the activated 
conda environment as follows (replace `assets/test-dardel.config` 
with your custom config file): 

```
nextflow run -c assets/test-dardel.config -profile dardel main.nf &> YYMMDD_rh.out
```

> The `dardel` profile is set up to submit each process as 
a job to the slurm queue, and to use Apptainer (former singularity) 
to run each process in a container with the required software. 

### Locally (small test dataset)

- Clone this repository to a directory on your computer

- Create the conda environment from `environment.yaml`, 
if you haven't done so yet (with mamba or conda): 

```
mamba env create -f environment.yaml
```

- Create a sample sheet in `assets/`, listing metadata 
information for each sample, including the location of raw 
fastq files with the sequencing data. Here is an example: 
`assets/samplesheet.csv`. Note that the pipeline currently 
only supports paired-end data.

- Create a custom pipeline configuration file with paths 
to input data, pipeline steps to be run, path to results 
directory, and tool-specific parameters. `assets/custom.config` 
is a template, a filled-out example is available here: 
`assets/test-local-eager_verena.config`.

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
conda environment as follows (replace `assets/test-local.config` 
with your custom config file): 

```
nextflow run -c assets/test-local.config -profile docker main.nf &> YYMMDD_rh.out
```

> The `docker` profile is set up to to use Docker to run each 
process in a container with the required software. It requires 
Docker to be installed and running on your local machine.
Currently, `-profile mamba` or `-profile conda` can only be 
used on Linux machines to run the pipeline because the tool 
ANGSD is only available as conda package for Linux. 
