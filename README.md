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
* reshape2
* grid
* ggtext
* stringr
## Clone the Repository
Clone the repository to your local machine:
```
git clone https://github.com/BreadcrumbsMulti-Tissue_Analysis_DL_AD/Multi-Tissue_Analysis_DL_AD.git
cd Multi-Tissue_Analysis_DL_AD
```
## Install Required R Packages
Install the necessary R packages to run the code in this repository:  
* CRAN Packages:
```
install.packages(c("dplyr", "ggplot2", "RColorBrewer", "data.table", "tidyr", "reshape2", "grid", "ggtext", "stringr"))
```
* Bioconductor Packages:
  * First, ensure BiocManager is installed:
    ```
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    ```
  * Then, install the required Bioconductor packages:
    ```
    BiocManager::install(c("VariantAnnotation", "BSgenome.Hsapiens.UCSC.hg19", "SNPlocs.Hsapiens.dbSNP144.GRCh37", "GenomicRanges", "qvalue", "biomaRt", "clusterProfiler"))
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
* Make_Significant_Counts_Normalized_Ratio.R:After Q-value correction, 18 traits passed the 0.01 FDR threshold. Significant gene counts across tissues were visualized using Z-score normalization.
* ADTraits_Tissue_GeneCount.R: Visualizing the distribution and total number of significant genes across 49 tissues for five key AD-related traits.
* ADTraits_GeneRegulationCount.R: The number of genes that are up-regulated, down-regulated, and ambiguous in each AD-related trait.
#### Gene Expression and Enrichment Analysis for AD Family History-Related Traits
We performed a comprehensive gene expression analysis of five Alzheimer’s disease (AD) family history-related traits based on the S-PrediXcan method.
* Make_AD_Traits_associated_gene_plots.R: Plots associated gene networks.
* GO_KEGG_Enrichment_Analysis: Performed GO and KEGG enrichment using the clusterProfiler package.
### 3. Individual-level PrediXcan Analysis of the ROSMAP Dataset and APOE Association Results
#### Training of gene expression prediction model and database construction
* Gene Expression Model Training and Database Construction: Using genotype and gene expression data from 553 ROSMAP samples, we trained gene expression prediction models and formatted them for PrediXcan using the PredictDB method. In accordance with GTEx guidelines, 60 PEER factors were applied. Full details can be found at [PredictDBPipeline](https://github.com/hakyimlab/PredictDBPipeline).
#### Model Training and PrediXcan Analysis at Individual Level
*Using trained gene expression models, we applied individual-level PrediXcan to predict gene expression and assess its association with key phenotypes such as APOE genotype, age at death, pathology scores, and cognitive diagnosis. This approach enabled us to evaluate the influence of gene expression on these traits. For detailed methods, see [PrediXcan](https://github.com/hakyimlab/PrediXcan).
#### Gene Expression Network Analysis of AD-Related Phenotypes
* Make_Brain_disorder_associated_gene_plots.R: Performs shared network analysis of gene expression for 9 Alzheimer’s disease (AD)-related phenotypes.
#### Functional Enrichment Analysis of APOE Genotype-Related Genes
* GO_KEGG_Enrichment_Analysis: Performed GO and KEGG enrichment using the clusterProfiler package.
### 4. Drug Repositioning Analysis and Potential Therapeutic Drugs for Alzheimer's Disease
#### Drug Repositioning via Gene Expression Correlation
Based on the method by Hon-Cheong So et al., we compared disease-related gene expression changes with drug-induced gene expression profiles to identify potential drug repositioning candidates.
* R codes used to compare drug and disease expression profiles, as well as ranking tests for evaluating the significance of drug scores, can be found online at:
https://sites.google.com/site/honcheongso/software/gwascmap.
* DrugIDTransform.R：Converts sig_id from the LINCS L1000 dataset to pert_iname to obtain standard drug names.
* Drug Repositioning Analysis for AD Family History-Related Traits
  * getTissueDrugHead50.R: Extracts the top 50 most potential drugs for each trait and tissue.
  * UpSetR.R: Performs intersection analysis to reveal common drugs and tissues among three AD family history traits.
### 5. Drug Repositioning Enrichment Analysis Method
To evaluate drug set enrichment in transcriptome-wide analysis, we used self-contained and competitive tests based on Hon-Cheong So et al.'s method, along with random drug sets for comparison.
* drug_set_enrich_KEGG.R: Evaluates drug repositioning results using ATC drug classifications from KEGG, including antiepileptic, antiparkinson, psychotropic, antipsychotic, and antidepressant/antianxiety drugs.
* drug_set_enrich_MEDI.R: Validates drug repositioning results using the MEDI-HPS dataset, which contains high-precision drug indication information.
* drug_set_enrich_ClinicalTrials.R: Assesses drug repositioning results based on a list of drugs from ClinicalTrials.gov.

### 6. AD-MIF

#### Requirement
- Pytorch --- 1.12.1
- Python --- 3.8.16
- Numpy --- 1.24.3
- Scipy --- 1.10.1
- Sklearn --- 1.2.2
- Munkres --- 1.1.4
- tqdm ---4.65.0
- Matplotlib ---3.7.1

#### Usage

#### Code structure
- ```data_loader.py```: loads the dataset and construct the gen graph
- ```opt.py```: defines parameters
- ```utils.py```: defines the utility functions
- ```encoder.py```: defines the AE and GAE
- ```AD-MIF.py```: defines the architecture of network
- ```main.py```: run the model
#### Example command

- Formal training model 
```
python main.py --name AD
```

## Contributing
We welcome contributions from the community. If you wish to contribute:
1. Fork the repository.  
2. Create a new branch (git checkout -b feature-branch).  
3. Make your changes and commit (git commit -am 'Add new feature').  
4. Push to the branch (git push origin feature-branch).  
5. Create a Pull Request.
## License
This project is licensed under the MIT License. See the LICENSE file for details.
