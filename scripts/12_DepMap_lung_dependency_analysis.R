# 12_DepMap_lung_dependency_analysis.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates how DepMap CRISPR dependency data
# can be filtered for lung cancer cell lines and summarized at
# the gene level.
#
# Notes:
# - This is a demo workflow using the depmap R package.
# - The threshold of mean_gene_effect <= -0.5 is used as an
#   example criterion for potential gene dependency.
# ------------------------------------------------------------

library(depmap)
library(dplyr)
library(readr)

# -----------------------------
# Load DepMap data
# -----------------------------
# Example objects often include:
# depmap::metadata
# depmap::crispr

meta <- depmap::metadata
crispr_data <- depmap::crispr

# -----------------------------
# Select lung lineage cell lines
# -----------------------------
lung_cell_line_ids <- meta %>%
  filter(grepl("lung", lineage, ignore.case = TRUE)) %>%
  pull(depmap_id)

lung_crispr <- crispr_data %>%
  filter(depmap_id %in% lung_cell_line_ids)

# Ensure dependency scores are numeric
lung_crispr <- lung_crispr %>%
  mutate(dependency = as.numeric(dependency))

# -----------------------------
# Summarize gene dependency across lung cell lines
# -----------------------------
gene_dependency_summary <- lung_crispr %>%
  group_by(gene_name) %>%
  summarise(
    mean_gene_effect = mean(dependency, na.rm = TRUE),
    n_cell_lines = sum(!is.na(dependency)),
    essential_lung = if_else(
      is.na(mean_gene_effect),
      NA,
      mean_gene_effect <= -0.5
    ),
    .groups = "drop"
  ) %>%
  arrange(mean_gene_effect)

# -----------------------------
# Save output
# -----------------------------
# output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  gene_dependency_summary,
  file = file.path(output_dir, "DepMap_lung_gene_dependency_summary.csv"),
  row.names = FALSE
)
