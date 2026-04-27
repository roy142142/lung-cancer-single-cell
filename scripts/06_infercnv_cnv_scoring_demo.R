# 06_infercnv_cnv_scoring_demo.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates a simplified inferCNV-based workflow
# for estimating CNV-like signals in epithelial tumor cells and
# calculating a per-cell CNV score.
#
# Notes:
# - This is a demo template.

# ------------------------------------------------------------

library(infercnv)
library(Seurat)

# -----------------------------
# User-defined input placeholders
# -----------------------------
# raw_counts_matrix <- "path/to/raw_counts_matrix.txt"
# annotations_file <- "path/to/cell_annotations.txt"
# gene_order_file <- "path/to/gene_order_file.txt"
# infercnv_output_dir <- "results/infercnv_demo"
# reference_group_name <- "Normal"

# -----------------------------
# Create inferCNV object
# -----------------------------
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts_matrix,
  annotations_file = annotations_file,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = reference_group_name
)

# -----------------------------
# Run inferCNV
# -----------------------------
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = infercnv_output_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 4
)

# -----------------------------
# Calculate per-cell CNV score
# -----------------------------
cnv_matrix <- infercnv_obj@expr.data

cnv_scores <- apply(
  cnv_matrix,
  2,
  function(cell_values) {
    mad(cell_values, na.rm = TRUE)
  }
)

# Add CNV scores to a Seurat object.
# The Seurat object should contain cells matching the CNV matrix columns.
# tumor_obj <- AddMetaData(tumor_obj, metadata = cnv_scores, col.name = "CNV_Score")

# -----------------------------
# Optional gene-level CNV summary
# -----------------------------
# cnv_matrix_subset can be defined by selecting tumor cells or
# specific epithelial subpopulations.
# cnv_matrix_subset <- cnv_matrix[, selected_cells]

# gene_cnv_means <- rowMeans(cnv_matrix_subset, na.rm = TRUE)
# hist(gene_cnv_means, breaks = 50, col = "grey",
#      main = "Distribution of Mean Gene-Level CNV Signal",
#      xlab = "Mean CNV-like signal")
