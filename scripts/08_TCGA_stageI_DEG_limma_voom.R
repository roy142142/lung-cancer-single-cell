# 08_TCGA_stageI_DEG_limma_voom.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates bulk RNA-seq differential expression
# analysis between Stage I lung adenocarcinoma tumor samples and
# normal tissue samples using limma-voom.
#
# Notes:
# - This is a demo template.
# - The input object `data_expr` should be a SummarizedExperiment-like
#   object containing count data and sample metadata.
# ------------------------------------------------------------

library(edgeR)
library(limma)
library(SummarizedExperiment)

# -----------------------------
# User-defined input placeholder
# -----------------------------
# data_expr <- readRDS("path/to/public_TCGA_or_demo_object.rds")

meta <- as.data.frame(colData(data_expr))

# -----------------------------
# Select Stage I tumor and normal samples
# -----------------------------
keep_stageI <- meta$shortLetterCode == "TP" &
  meta$ajcc_pathologic_stage %in% c("Stage I", "Stage IA", "Stage IB")

keep_normal <- meta$shortLetterCode == "NT"

data_stageI <- data_expr[, keep_stageI]
data_normal <- data_expr[, keep_normal]

data_filtered <- cbind(data_stageI, data_normal)

# -----------------------------
# Define sample groups
# -----------------------------
filtered_meta <- as.data.frame(colData(data_filtered))

group <- ifelse(
  filtered_meta$shortLetterCode == "TP",
  "Tumor",
  "Normal"
)

group <- factor(group, levels = c("Normal", "Tumor"))

# -----------------------------
# Extract count matrix
# -----------------------------
expr_counts <- assay(data_filtered)

# Optional low-expression filtering
keep_genes <- rowSums(expr_counts > 10) >= 10
expr_counts <- expr_counts[keep_genes, ]

# -----------------------------
# DEG analysis using limma-voom
# -----------------------------
dge <- DGEList(counts = expr_counts)
dge <- calcNormFactors(dge)

design <- model.matrix(~ group)

voom_obj <- voom(dge, design, plot = TRUE)

fit <- lmFit(voom_obj, design)
fit <- eBayes(fit)

deg_results <- topTable(
  fit,
  coef = "groupTumor",
  number = Inf,
  sort.by = "P"
)

deg_results$Gene <- rownames(deg_results)

# -----------------------------
# Save output
# -----------------------------
# output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  deg_results,
  file = file.path(output_dir, "TCGA_StageI_Tumor_vs_Normal_limma_voom_DEG.csv"),
  row.names = FALSE
)
