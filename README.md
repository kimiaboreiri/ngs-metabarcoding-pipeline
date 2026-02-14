# Metabarcoding Analysis Pipeline

This repository presents a metabarcoding and eDNA bioinformatics workflow for biodiversity assessment using next-generation sequencing (NGS) data.

##  Project Overview
Metabarcoding uses marker-gene sequencing to detect species from environmental DNA (eDNA).  
This pipeline reflects real-world metabarcoding analysis used in research and biodiversity monitoring.

##  Workflow
- Quality control of FASTQ files  
- Adapter and primer trimming  
- ASV generation  
- Taxonomic assignment  
- Diversity analysis  

##  Tools
- QIIME2  
- DADA2  
- Cutadapt
- Python
- R for visualization  

##  Applications
- Biodiversity monitoring  
- Environmental DNA studies  
- Ecological research  

##  Note
A metadata template is provided here. Please fill it in according to your sample names and negative blanks.
If you plan to use this metadata template, please follow the sample naming conventions below.
- Extraction blanks start with "EB"  
- PCR negative controls start with "neg"
  
This pipeline was tested on Fox River vertebrate and mussel metabarcoding data for biodiversity analysis. Depending on your target species, you may need to build a specific reference database.
