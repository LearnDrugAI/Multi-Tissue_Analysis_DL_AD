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
* VariantAnnotation 
* BSgenome.Hsapiens.UCSC.hg19
* SNPlocs.Hsapiens.dbSNP144.GRCh37
* GenomicRanges
* dplyr
* qvalue 
* biomaRt 
* clusterProfiler 
* ggplot2
* RColorBrewer
* data.table
* RColorBrewer
* tidyr
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
#### Drug-induced Gene Expression Data
Sourced from the LINCS L1000 project (CMap), measuring 978 landmark genes to infer genome-wide expression. These data were used for analyzing drug effects on gene expression and drug repositioning studies.  
* filt_trt_cp_for_cmap.R:  Filtered "trt_cp" (drug treatment) samples and used the cmapR package to extract and save the expression statistics matrix for further analysis.
### 2. Association between genetic prediction gene expression level and 61 brain-related characteristics
We used fastGWA-GLMM GWAS data from Jiang et al. and UK Biobank to analyze the association between predicted gene expression (from 49 GTEx tissues) and 61 brain-related traits using S-PrediXcan. Results were adjusted for multiple comparisons, with significant associations visualized via heatmaps and normalized significance counts.Methods are based on [SPrediXcan and MetaXcan](https://github.com/hakyimlab/MetaXcan).
#### Gene Expression and Significant Gene Count for AD Family History-Related Traits
After Q-value correction, 18 traits passed the 0.01 FDR threshold. Significant gene counts across tissues were visualized using Z-score normalization.
Make_Significant_Counts_Normalized_Ratio.R:After Q-value correction, 18 traits passed the 0.01 FDR threshold. Significant gene counts across tissues were visualized using Z-score normalization.
ADTraits_Tissue_GeneCount: Visualizing the distribution and total number of significant genes across 49 tissues for five key AD-related traits.

