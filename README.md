# A year in single-cell RNA sequencing data analysis

> This repo contains what I have learnt the last year working on single cell RNA sequencing data using R. 

The project I was involved with was to analyze CD19+ sorted cells derived from patients with a B-cell lymphoproliferative disorder and healthy donors (Droplet-based 10X genomics sequencing library preparation).
The Cell Ranger software (version 4.0.0) was used to demultiplex raw sequencing data of patients and the healthy donor and align reads to the GRCh38 human reference genome using default parameters.

The rscript_draft.R will guide you through the steps of analysing each sample individually and then integratng all datasets togehter. In the folder /functions you can find the functions I have created.

I have explored the following Bioconductor packages:
- Seurat (version 4.2.0) 
- harmony (version 0.1.0) 
- monocle (version 2.18)
- clusterProfiler (version 3.18.1)
- slingshot (version 2.10.0)
- inferCNV
R program environment (version 4.0.2)
