# 02_normalization_clustering_umap.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates a standard Seurat workflow for
# normalization, variable feature selection, scaling, PCA,
# neighbor graph construction, clustering, and UMAP visualization.
#
# Notes:
# - This script assumes that `seurat_obj` has already passed QC.
# - The variables `dataset_name` and `output_dir` should be defined.
# ------------------------------------------------------------

library(Seurat)
library(ggplot2)

dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Normalization and feature selection
# -----------------------------
seurat_obj <- NormalizeData(seurat_obj)

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# -----------------------------
# Scaling and PCA
# -----------------------------
seurat_obj <- ScaleData(seurat_obj)

seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(seurat_obj),
  npcs = 30
)

# -----------------------------
# Clustering and UMAP
# -----------------------------
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)

seurat_obj <- FindClusters(
  seurat_obj,
  resolution = 0.6
)

seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# -----------------------------
# UMAP plot
# -----------------------------
umap_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("UMAP Clustering")

ggsave(
  filename = file.path(output_dir, "plots", paste0(dataset_name, "_UMAP_Clusters.png")),
  plot = umap_plot,
  width = 7,
  height = 5,
  dpi = 300
)
