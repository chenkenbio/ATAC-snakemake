# ATAC-seq snakemake pipeline

## Usage  

1. Editing the following files to specify paths to references or programs:  
- `src/alignment_qc_v0.smk`  
- `src/fastq_screen.smk`  
- `src/blasn_screen.smk`  
- `src/picard_markdup.smk`  
- `src/filter_bam.smk`  
- `src/bamcoverage_atac.smk`  


2. Run the pipeline
```bash
snakemake -s pipeline_ATAC.smk -c [n_cores]
```




