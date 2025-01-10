# Notes on bwa module implementation

Currently, the pipeline only works with merged paired-end or with single-end data.

bwa index
- Option for large (human-sized) genome assemblies available in configs/modules.config (commented out for test runs with small reference)

bwa aln
- Parameters for ancient DNA in configs/modules.config. To be placed into configs/custom.config?

bwa samse
- Readgroup is automatically extracted from the samplesheet