# 09_gene_annotation_GO_mapping.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates mapping gene symbols to Entrez IDs
# and retrieving Gene Ontology annotations using org.Hs.eg.db.
#
# Notes:
# - This is a general demo script.
# - The input data frame `deg_table` should contain a column named `Gene`.
# ------------------------------------------------------------

library(AnnotationDbi)
library(org.Hs.eg.db)

# -----------------------------
# User-defined input placeholder
# -----------------------------
# deg_table <- read.csv("path/to/demo_DEG_table.csv")

if (!"Gene" %in% colnames(deg_table)) {
  stop("Input DEG table must contain a column named 'Gene'.")
}

# -----------------------------
# Map gene symbols to Entrez IDs
# -----------------------------
entrez_ids <- mapIds(
  x = org.Hs.eg.db,
  keys = deg_table$Gene,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

deg_table$EntrezID <- entrez_ids

deg_table_mapped <- deg_table[!is.na(deg_table$EntrezID), ]

# -----------------------------
# Retrieve GO annotations
# -----------------------------
go_annotations <- AnnotationDbi::select(
  x = org.Hs.eg.db,
  keys = deg_table_mapped$EntrezID,
  columns = c("GO", "ONTOLOGY"),
  keytype = "ENTREZID"
)

# -----------------------------
# Save outputs
# -----------------------------
# output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  deg_table_mapped,
  file = file.path(output_dir, "DEG_with_EntrezID.csv"),
  row.names = FALSE
)

write.csv(
  go_annotations,
  file = file.path(output_dir, "DEG_GO_annotations.csv"),
  row.names = FALSE
)
