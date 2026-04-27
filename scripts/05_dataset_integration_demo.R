# 05_dataset_integration_demo.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates a Seurat anchor-based integration
# workflow for multiple single-cell RNA-seq datasets.
#
# Notes:
# - This is a simplified demo.
# - `seurat_object_list` should contain multiple preprocessed
#   Seurat objects, for example tumor and/or normal samples.
# - Parameter values such as k.filter, k.anchor, k.weight, and dims are examples
# ------------------------------------------------------------

library(Seurat)

# -----------------------------
# seurat_object_list <- list(sample1_obj, sample2_obj, sample3_obj)
# -----------------------------
# Select integration features
# -----------------------------
integration_features <- SelectIntegrationFeatures(
  object.list = seurat_object_list,
  nfeatures = 3000
)

# -----------------------------
# Find integration anchors
# -----------------------------
integration_anchors <- FindIntegrationAnchors(
  object.list = seurat_object_list,
  anchor.features = integration_features,
  k.filter = 10,
  k.anchor = 30,
  dims = 1:30
)

# -----------------------------
# Integrate datasets
# -----------------------------
integrated_obj <- IntegrateData(
  anchorset = integration_anchors,
  k.weight = 10,
  new.assay.name = "integrated",
  dims = 1:30
)

DefaultAssay(integrated_obj) <- "integrated"

# -----------------------------
# Downstream dimensionality reduction
# -----------------------------
integrated_obj <- ScaleData(integrated_obj)
integrated_obj <- RunPCA(integrated_obj, npcs = 30)
integrated_obj <- FindNeighbors(integrated_obj, dims = 1:30)
integrated_obj <- FindClusters(integrated_obj, resolution = 0.5)
integrated_obj <- RunUMAP(integrated_obj, dims = 1:30)
