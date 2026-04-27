# 03_singleR_annotation.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates automated cell-type annotation of a
# Seurat object using SingleR and the BlueprintEncodeData reference.
#
# Notes:
# - This is a demo workflow.
# - Annotation quality should always be checked manually using
#   known marker genes and biological context.
# ------------------------------------------------------------

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)

dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load reference dataset
# -----------------------------
ref <- BlueprintEncodeData()

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(seurat_obj)

# -----------------------------
# Run SingleR annotation
# -----------------------------
annotation_result <- SingleR(
  test = sce,
  ref = ref,
  labels = ref$label.main,
  fine.tune = TRUE
)

# Store pruned labels in Seurat metadata
seurat_obj$SingleR_Labels <- annotation_result$pruned.labels

# -----------------------------
# UMAP annotation plot
# -----------------------------
annotation_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "SingleR_Labels",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("UMAP Annotated by SingleR")

ggsave(
  filename = file.path(output_dir, "plots", paste0(dataset_name, "_UMAP_SingleR_Annotated.png")),
  plot = annotation_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# Annotation summary table
# -----------------------------
annotation_counts <- as.data.frame(table(seurat_obj$SingleR_Labels))
colnames(annotation_counts) <- c("Cell_Type", "Cell_Count")

write.csv(
  annotation_counts,
  file = file.path(output_dir, paste0(dataset_name, "_SingleR_annotation_counts.csv")),
  row.names = FALSE
)
