# 07_scRNA_tumor_vs_normal_DEG.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates differential expression analysis
# between tumor and normal cells in a merged Seurat object.
#
# Notes:
# - This script assumes that `seurat_merged`, `tumor_obj`, and
#   `normal_obj` or equivalent cell groups have already been created.
# - MAST is used as an example test method.
# ------------------------------------------------------------

library(Seurat)

# -----------------------------
# Define tumor/normal group labels
# -----------------------------
# Example:
# tumor_cells <- colnames(tumor_obj)
# normal_cells <- colnames(normal_obj)

seurat_merged$group <- ifelse(
  colnames(seurat_merged) %in% tumor_cells,
  "Tumor",
  "Normal"
)

seurat_merged$group <- factor(seurat_merged$group, levels = c("Normal", "Tumor"))

Idents(seurat_merged) <- seurat_merged$group

# -----------------------------
# Differential expression analysis
# -----------------------------
deg_result <- FindMarkers(
  object = seurat_merged,
  ident.1 = "Tumor",
  ident.2 = "Normal",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "MAST",
  slot = "data",
  assay = "RNA",
  only.pos = FALSE
)

deg_result$Gene <- rownames(deg_result)

# -----------------------------
# Save DEG table
# -----------------------------
# output_dir <- "results"
# dataset_name <- "demo_dataset"

write.csv(
  deg_result,
  file = file.path(output_dir, paste0(dataset_name, "_Tumor_vs_Normal_DEG_scRNA.csv")),
  row.names = FALSE
)
