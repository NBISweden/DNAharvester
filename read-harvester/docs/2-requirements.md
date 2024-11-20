# Pipeline requirements

This pipeline analyzes whole genome re-sequencing data and 
requires a compute system with at least 64 GB RAM to analyze 
data for a species with a genome size of around 3 Gbp. 

Conda/mamba has to be installed on the system to start a 
pipeline run. The pipeline processes can be run with Docker, 
Apptainer, or conda/mamba. Readharvester has been tested on 
a HPC system. 

Make sure there is storage space available as the pipeline 
produces a lot of large files. 

Clone this repository to the directory where you want to run 
the pipeline (hereafter called `path/to/your/readharvester_run`). 
All output files will be automatically written to subdirectories 
in `results/`. Intermediate files are stored in a folder called 
`work/` that can be safely removed once the pipeline has finished. 

Create a conda environment from the yaml file 
`path/to/your/readharvester_run/environment.yaml`. The environment 
contains Nextflow and nf-core: 

```
conda env create -f environment.yaml
````
