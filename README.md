# Integrative Transcriptome-Wide Analysis for Alzheimer’s Disease
This repository contains the code and methods used in the research paper "Integrative Transcriptome-Wide and Multimodal Deep Learning Analyses Reveal APOE-Associated Gene Expression in Alzheimer’s Disease". The study utilizes multi-tissue transcriptome-wide association studies (TWAS) and deep learning techniques to identify key genes and potential drug repositioning candidates related to Alzheimer's Disease (AD).
## Table of Contents
* Introduction
* Requirements
* Installation
* Code Structure
  * Data Preprocessing
  * TWAS Analysis
  * GO and KEGG Pathway Enrichment
  * Drug Repositioning Analysis
  * Visualization
* Contributing
* License
## Introduction
This project integrates transcriptome-wide association studies (TWAS) from multiple tissues using R-based pipelines to identify Alzheimer's Disease-related genes and conduct pathway and drug repositioning analyses. We highlight the roles of genes such as APOE, APOC1, and others in AD pathology, and investigate potential therapeutic drugs using gene expression profiles.
## Requirements
To run the R scripts in this repository, you need the following R packages:
* tidyverse 
* plink2R 
* qvalue
* TopmedPipeline 
* biomaRt 
* clusterProfiler 
* ggplot2
* RColorBrewer
* data.table
## Installation
Clone the repository to your local machine:
```
git clone https://github.com/BreadcrumbsMulti-Tissue_Analysis_DL_AD/AD-TWAS-Analysis-R.git
cd AD-TWAS-Analysis-R
```
## Code Structure
The code is divided into functional modules, each serving a specific purpose in the overall analysis workflow.

### 1. Data Preprocessing
Before running the main analysis, raw data must be preprocessed.
#### GWAS Summary Statistics
This project utilizes GWAS summary statistics from UK Biobank, analyzed using the fastGWA-GLMM method by Jiang et al. The data include 61 traits related to psychiatric disorders, neurological diseases, sleep patterns, and metabolic conditions, with sample sizes ranging from hundreds to hundreds of thousands. All samples are of European ancestry and have undergone rigorous quality control.
* Coordinate Conversion (hg19 to hg38)
GWAS data were converted from the human genome version 19 (hg19) to version 38 (hg38) for consistency with modern genomic resources.
* Genotype Imputation  
Missing genotype data were imputed using the 1000 Genomes Project reference panel, with a 100kb sliding window approach and standardized dosage information to improve imputation accuracy.
* Variant Filtering  
Variants were filtered based on:  
  * Minor Allele Frequency (MAF) > 0.01
  * Genotyping call rate > 0.8
* Tools Used
  * Liftover tool for genome coordinate conversion.
  * Imputation software for genotype imputation using the 1000 Genomes reference.
For more details, refer to the [Hakyimlab Summary-GWAS-Imputation repository](https://github.com/hakyimlab/summary-gwas-imputation/)
#### ROSMAP Dataset
This study uses data sets from Religious Orders Study and Memory and Aging Project (ROSMAP), which is managed by Rush Alzheimer's Disease Center in Chicago.  
Scripts:
Filtered out non-SNVs and multi-allelic loci, then converted the genotype data into an SNP matrix for association analysis.
* preprocess_filter_and_convert_SNP_matrix.R: Prepares genotype data for TWAS analysis by converting file formats, filtering SNPs, and handling missing values.
* preprocess_get_SNPAnnot.R:Annotated SNPs using the UCSC hg19 reference genome and dbSNP137, extracting RSIDs and adding them to the information table.
* preprocess_get_GeneAnnot.R：Retrieved gene chromosomal positions, IDs, and types from the Ensembl database, filtering out non-autosomal genes to retain only autosomal annotations.
