# 04_epithelial_subsetting_subtyping.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates how epithelial cells can be extracted
# from a SingleR-annotated Seurat object and re-clustered for
# marker-based epithelial subtype annotation.
#
# Notes:
# - This is a simplified demo version.
# - Marker-based subtype calls should be interpreted carefully
#   and validated using multiple markers and biological context.
# ------------------------------------------------------------

library(Seurat)
library(ggplot2)

dir.create(file.path(output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Extract epithelial-like cells
# -----------------------------
epi_keywords <- c(
  "Epithelial",
  "Basal",
  "Luminal",
  "Club",
  "Alveolar",
  "Goblet"
)

is_epithelial <- grepl(
  pattern = paste(epi_keywords, collapse = "|"),
  x = seurat_obj$SingleR_Labels,
  ignore.case = TRUE
)

epi_obj <- subset(
  seurat_obj,
  cells = colnames(seurat_obj)[is_epithelial]
)

# -----------------------------
# Re-cluster epithelial cells
# -----------------------------
epi_obj <- NormalizeData(epi_obj)
epi_obj <- FindVariableFeatures(epi_obj, selection.method = "vst", nfeatures = 2000)
epi_obj <- ScaleData(epi_obj)
epi_obj <- RunPCA(epi_obj, features = VariableFeatures(epi_obj), npcs = 50)
epi_obj <- FindNeighbors(epi_obj, dims = 1:30)
epi_obj <- FindClusters(epi_obj, resolution = 0.5)
epi_obj <- RunUMAP(epi_obj, dims = 1:30)

# -----------------------------
# Marker-based epithelial subtype scoring
# -----------------------------
epithelial_markers <- list(
  Club = c("SCGB1A1", "SCGB3A1"),
  Basal = c("KRT5", "TP63", "KRT14"),
  Ciliated = c("FOXJ1", "DNAH5", "TUBA1A"),
  Goblet = c("MUC5AC", "SPDEF", "AGR2"),
  AT1 = c("AGER", "PDPN"),
  AT2 = c("SFTPC", "SFTPA1", "SFTPB")
)

epi_obj$Subtype <- "Unknown"

for (cell_type in names(epithelial_markers)) {
  marker_genes <- epithelial_markers[[cell_type]]
  marker_genes <- marker_genes[marker_genes %in% rownames(epi_obj)]

  if (length(marker_genes) > 0) {
    expression_matrix <- GetAssayData(epi_obj, slot = "data")
    subtype_score <- colMeans(expression_matrix[marker_genes, , drop = FALSE])
    epi_obj[[paste0("score_", cell_type)]] <- subtype_score
  }
}

score_columns <- paste0("score_", names(epithelial_markers))
score_columns <- score_columns[score_columns %in% colnames(epi_obj@meta.data)]

scores_df <- FetchData(epi_obj, vars = score_columns)

# Example threshold. This should be adjusted based on data distribution.
minimum_score_threshold <- 0.5

epi_obj$Subtype <- apply(scores_df, 1, function(score_vector) {
  if (all(is.na(score_vector)) || max(score_vector, na.rm = TRUE) < minimum_score_threshold) {
    return("Unknown")
  }

  best_score_name <- names(score_vector)[which.max(score_vector)]
  gsub("score_", "", best_score_name)
})

# -----------------------------
# UMAP plot of epithelial subtypes
# -----------------------------
subtype_plot <- DimPlot(
  epi_obj,
  reduction = "umap",
  group.by = "Subtype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Epithelial Subtypes by Marker Expression")

ggsave(
  filename = file.path(output_dir, "plots", paste0(dataset_name, "_UMAP_Epithelial_Subtypes.png")),
  plot = subtype_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# -----------------------------
# Cell count summary
# -----------------------------
cell_counts <- data.frame(
  Dataset = dataset_name,
  Total_Epithelial_Cells = ncol(epi_obj),
  Subtype_Club = sum(epi_obj$Subtype == "Club"),
  Subtype_Basal = sum(epi_obj$Subtype == "Basal"),
  Subtype_Ciliated = sum(epi_obj$Subtype == "Ciliated"),
  Subtype_Goblet = sum(epi_obj$Subtype == "Goblet"),
  Subtype_AT1 = sum(epi_obj$Subtype == "AT1"),
  Subtype_AT2 = sum(epi_obj$Subtype == "AT2"),
  Subtype_Unknown = sum(epi_obj$Subtype == "Unknown")
)

write.csv(
  cell_counts,
  file = file.path(output_dir, paste0(dataset_name, "_epithelial_subtype_counts.csv")),
  row.names = FALSE
)
