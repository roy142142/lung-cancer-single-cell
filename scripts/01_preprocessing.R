# 01_preprocessing.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates the initial quality-control workflow
# for a single-cell RNA-seq dataset using Seurat.
#
# Notes:
# - This is a demo version prepared for GitHub use.
# - The input object `data_matrix` should be a gene-by-cell count matrix.
# ------------------------------------------------------------

library(Seurat)
library(ggplot2)


# -----------------------------
# data_matrix <- readRDS("path/to/public_or_demo_count_matrix.rds")
# dataset_name <- "demo_dataset"
# output_dir <- "results"

dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Create Seurat object
# -----------------------------
seurat_obj <- CreateSeuratObject(
  counts = data_matrix,
  project = "LungCancerDemo",
  min.cells = 3,
  min.features = 200
)

total_cells_before_qc <- ncol(seurat_obj)

# -----------------------------
# Add QC metrics
# -----------------------------
# Mitochondrial gene patterns can differ depending on annotation style.
# Human datasets commonly use either "MT-" or occasionally "MT.".
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-|^MT\\."
)

seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^RPL|^RPS"
)

# -----------------------------
# QC visualization
# -----------------------------
qc_plot <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  ncol = 4,
  pt.size = 0.1
)

ggsave(
  filename = file.path(output_dir, "plots", paste0(dataset_name, "_QC_ViolinPlot.png")),
  plot = qc_plot,
  width = 10,
  height = 4,
  dpi = 300
)

# -----------------------------
# QC filtering
# -----------------------------
# These thresholds are examples and should be adjusted based on
# dataset-specific distributions and biological context.
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 3000 &
           percent.mt < 15 &
           percent.ribo < 100
)

total_cells_after_qc <- ncol(seurat_obj)

qc_summary <- data.frame(
  Dataset = dataset_name,
  Total_Cells_Before_QC = total_cells_before_qc,
  Total_Cells_After_QC = total_cells_after_qc,
  Cells_Removed = total_cells_before_qc - total_cells_after_qc
)

write.csv(
  qc_summary,
  file = file.path(output_dir, paste0(dataset_name, "_QC_summary.csv")),
  row.names = FALSE
)
