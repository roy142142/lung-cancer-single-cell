# R Package Versions

This file lists the main R packages used in the original analysis workflow.  
The repository contains a simplified demonstration version of the workflow, so not every package listed here is required for every script.
The versions below reflect the computational environment used during the analysis.

| Package | Version |
|---|---:|
| AnnotationDbi | 1.68.0 |
| biomaRt | 2.62.1 |
| celldex | 1.16.0 |
| cluster | 2.1.6 |
| clusterProfiler | 4.14.6 |
| cowplot | 1.2.0 |
| data.table | 1.17.8 |
| depmap | 1.20.0 |
| dplyr | 1.1.4 |
| edgeR | 4.4.2 |
| ggplot2 | 3.5.2 |
| ggrepel | 0.9.6 |
| glue | 1.8.0 |
| GO.db | 3.20.0 |
| gridExtra | 2.3 |
| gtexr | 0.2.0 |
| harmony | 1.2.3 |
| httr | 1.4.7 |
| infercnv | 1.22.0 |
| jsonlite | 1.9.1 |
| limma | 3.62.2 |
| lisi | 1 |
| magrittr | 2.0.3 |
| Matrix | 1.7.3 |
| openxlsx | 4.2.8 |
| org.Hs.eg.db | 3.20.0 |
| patchwork | 1.3.1 |
| pheatmap | 1.0.13 |
| pROC | 1.19.0.1 |
| purrr | 1.0.4 |
| readr | 2.1.5 |
| recount3 | 1.16.0 |
| Seurat | 5.3.0 |
| SeuratObject | 5.1.0 |
| SingleCellExperiment | 1.28.1 |
| SingleR | 2.8.0 |
| stats | 4.4.2 |
| stringr | 1.5.1 |
| SummarizedExperiment | 1.36.0 |
| survival | 3.7.0 |
| survminer | 0.5.0 |
| TCGAbiolinks | 2.34.1 |
| tibble | 3.2.1 |
| tidyr | 1.3.1 |
| tidyverse | 2.0.0 |

## Main Analysis Categories

### Single-cell RNA-seq analysis

- Seurat
- SeuratObject
- SingleCellExperiment
- SingleR
- celldex
- harmony
- lisi

### CNV inference from single-cell RNA-seq

- infercnv

### Bulk RNA-seq and TCGA-related analysis

- TCGAbiolinks
- SummarizedExperiment
- edgeR
- limma
- recount3

### Functional enrichment and annotation

- AnnotationDbi
- org.Hs.eg.db
- GO.db
- biomaRt
- clusterProfiler

### Biomarker and survival analysis

- survival
- survminer
- pROC

### Network and dependency analysis

- depmap
- cluster

### Data handling and visualization

- tidyverse
- dplyr
- tidyr
- tibble
- readr
- data.table
- stringr
- purrr
- ggplot2
- ggrepel
- cowplot
- patchwork
- pheatmap
- gridExtra
- openxlsx

## Notes

Package versions can affect object structure and function behavior, especially for Seurat v5, inferCNV, TCGAbiolinks, and Bioconductor packages.

For exact reproducibility, users should recreate a compatible R/Bioconductor environment or adapt the scripts to their local package versions.
