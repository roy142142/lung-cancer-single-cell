# lung-cancer-single-cell
Demo Workflow for single-cell and bulk transcriptomic analysis in lung adenocarcinoma.

## Overview

This repository presents a simplified demonstration workflow for single-cell and bulk transcriptomic analysis in lung adenocarcinoma. The project is inspired by my Master's research, where I worked on integrating single-cell RNA-seq and bulk RNA-seq data to investigate malignant epithelial cell states, copy-number variation patterns, differentially expressed genes, and candidate biomarkers in lung adenocarcinoma.

The aim of this repository is not to provide the full unpublished thesis pipeline, but to demonstrate the analytical logic, coding style, and reproducible structure behind the work.

## Research Motivation

Lung adenocarcinoma is a highly heterogeneous cancer, especially at early stages. Single-cell RNA-seq provides a powerful way to investigate tumor epithelial heterogeneity, identify malignant-like cell populations, and study gene expression programs at high resolution.

In this demo workflow, the main focus is on:

- Preprocessing single-cell RNA-seq data
- Annotating epithelial cell populations
- Comparing tumor and normal epithelial cells
- Prioritizing candidate genes using differential expression and CNV-related signals
- Connecting single-cell findings with bulk transcriptomic validation
- Demonstrating downstream survival analysis logic

## Main Analysis Steps

The full conceptual workflow includes:

1. **Single-cell data preprocessing**
   - Quality control
   - Normalization
   - Feature selection
   - Dimensionality reduction
   - Clustering

2. **Cell type annotation**
   - Marker-based epithelial cell identification
   - Refinement of tumor and normal epithelial populations

3. **Differential expression analysis**
   - Tumor versus normal epithelial comparison
   - Identification of candidate genes with altered expression

4. **CNV-related analysis**
   - Estimation of CNV-like signals from epithelial cells
   - Prioritization of high-CNV tumor-like populations
   - Integration of CNV-associated and expression-associated genes

5. **Bulk RNA-seq validation**
   - Integration with public bulk transcriptomic data
   - Expression validation in lung adenocarcinoma cohorts

6. **Survival analysis**
   - Kaplan-Meier survival analysis
   - Cox proportional hazards modelling
   - Candidate biomarker prioritization

## Repository Structure


```text
lung-cancer-single-cell/
│
├── README.md
├── scripts/
│   ├── 01_qc_preprocessing.R
│   ├── 02_normalization_clustering_umap.R
│   ├── 03_singleR_annotation.R
│   ├── 04_epithelial_subsetting_subtyping.R
│   ├── 05_dataset_integration_demo.R
│   ├── 06_infercnv_cnv_scoring_demo.R
│   ├── 07_scRNA_tumor_vs_normal_DEG.R
│   ├── 08_TCGA_stageI_DEG_limma_voom.R
│   ├── 09_gene_annotation_GO_mapping.R
│   ├── 10_ROC_biomarker_evaluation.R
│   ├── 11_PPI_network_centrality.R
│   └── 12_DepMap_lung_dependency_analysis.R
│
├── data/
│   └── README.md
│
├── figures/
│   └── README.md
│   └── Workflow_overview_public.png
│
└── requirements/
    └── package_versions.md
