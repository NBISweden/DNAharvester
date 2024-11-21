# How to run readharvester

## HPC cluster Dardel (PDC/KTH)

- Load the following modules:

```
module load PDC bioinfo-tools apptainer tmux
```

> Note that `tmux` is only available as a module on Dardel 
but the equivalent tool `screen` is pre-installed and does 
not need to be loaded. 

> Apptainer (former singularity) can use your `scratch` for 
caching, which is a temporary directory with unlimited space 
by adding this row to your `~/.bashrc`: 
`export NXF_SINGULARITY_CACHEDIR=$PDC_TMP` 

- Open a tmux or screen session, e.g.: 

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
a job to the slurm queue with compute resources specified in 
`configs/dardel.config`, and to use apptainer (former singularity) 
to run each process in a container with the required software. 

> Other available profiles: `conda`, `mamba`, `singularity`, 
`docker`. These profiles will use compute resources specified 
in `configs/nf-core-defaults.config`. Currently, `conda` and 
`mamba` can only be used on Linux machines to run the pipeline 
because the tool ANGSD is only available as conda package for 
Linux. 

## Locally (small test dataset)

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
with your custom config file and `docker` with your profile of 
choice): 

```
nextflow run -c assets/test-local.config -profile docker main.nf &> YYMMDD_rh.out
```

> The `docker` profile is set up to to use Docker to run each 
process in a container with the required software. It requires 
Docker to be installed and running on your local machine.
Currently, `conda` and `mamba` can only be used on Linux 
machines to run the pipeline because the tool ANGSD is only 
available as conda package for Linux. 