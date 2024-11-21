# Pipeline requirements

This pipeline analyzes whole genome re-sequencing data and 
requires a compute system with at least 64 GB RAM to analyze 
data for a species with a genome size of around 3 Gbp. 

Nextflow has to be installed, e.g. in the conda environment 
specified in `environment.yml`, to start a pipeline run. 
The pipeline processes can be run with Docker, apptainer 
(only on HPC clusters), or conda/mamba (currently only on 
Linux systems). Readharvester has been tested on Dardel, an 
HPC system with apptainer. 

A terminal multiplexer like tmux or screen is useful to send 
the Nextflow process to the background since it can take a 
while for the pipeline to finish.

Make sure there is storage space available as the pipeline 
produces a lot of large files. 

Clone this repository to the directory where you want to run 
the pipeline (hereafter called `path/to/your/readharvester_run`). 
All output files will be automatically written to subdirectories 
in `results/`. Intermediate files are stored in a folder 
called `work/` that can be safely removed once the pipeline 
has finished. 

Create a conda environment from the yaml file 
`path/to/your/readharvester_run/environment.yaml`. The 
environment contains Nextflow and nf-core: 

```
conda env create -f environment.yaml
````

> Note that conda environments can also be created in 
directories other than `/home/USER/` by running 
`export CONDA_ENVS_PATH=/path/to/project/directory/conda_environments/` followed by 
`conda env create -f environment.yml -p /path/to/project/directory/conda_environments/read-harvester`. 
This is recommended if you run readharvester on the HPC 
cluster Dardel (PDC/KTH). 

Readharvester runs AMBER, a tool that is not available as 
container image or conda package yet. To install AMBER, 
clone the AMBER Github repository `https://github.com/tvandervalk/AMBER/` 
into a different directory than readharvester and copy the 
file `AMBER` to `read-harvester/bin`. 