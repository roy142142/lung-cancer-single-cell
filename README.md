# lung-cancer-single-cell
Workflow for single-cell and bulk transcriptomic analysis in lung adenocarcinoma.

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
lung-cancer-single-cell-demo/
│
├── README.md
├── scripts/
│   ├── 01_preprocessing.R
│   ├── 02_clustering_annotation.R
│   ├── 03_differential_expression.R
│   ├── 04_cnv_score_demo.R
│   └── 05_survival_analysis_demo.R
│
├── figures/
│   ├── workflow_overview.png
│   ├── example_umap.png
│   └── example_survival_plot.png
│
├── data/
│   └── README.md
│
└── requirements/
    └── packages_used.txt
