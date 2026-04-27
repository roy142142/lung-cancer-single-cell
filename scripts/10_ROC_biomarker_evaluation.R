# 10_ROC_biomarker_evaluation.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates ROC-based evaluation of candidate
# biomarker genes for binary classification, for example tumor
# versus normal samples.
#
# Notes:
# - This is a demo version.
# - `expr_mat` should be a gene-by-sample expression matrix.
# - `class_vec` should be a binary vector indicating sample class.
# - `genes_to_test` should contain candidate gene symbols.
# ------------------------------------------------------------

library(pROC)

# -----------------------------
# Function for ROC evaluation
# -----------------------------
evaluate_genes_by_roc <- function(expr_mat, class_vec, genes_to_test) {
  results <- data.frame(
    Gene = character(),
    AUC = numeric(),
    Best_Cutoff = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    Accuracy = numeric(),
    Precision = numeric(),
    stringsAsFactors = FALSE
  )

  for (gene in genes_to_test) {
    if (!gene %in% rownames(expr_mat)) {
      next
    }

    gene_expression <- as.numeric(expr_mat[gene, ])

    roc_obj <- tryCatch(
      {
        roc(class_vec, gene_expression, quiet = TRUE)
      },
      error = function(e) {
        NULL
      }
    )

    if (is.null(roc_obj)) {
      next
    }

    auc_value <- as.numeric(auc(roc_obj))

    best_coords <- tryCatch(
      {
        coords(
          roc_obj,
          x = "best",
          ret = c("threshold", "sensitivity", "specificity")
        )
      },
      error = function(e) {
        NULL
      }
    )

    if (is.null(best_coords) || any(is.na(best_coords["threshold"]))) {
      cutoff <- NA
      sensitivity <- NA
      specificity <- NA
      accuracy <- NA
      precision <- NA
    } else {
      cutoff <- as.numeric(best_coords[["threshold"]])
      sensitivity <- as.numeric(best_coords[["sensitivity"]])
      specificity <- as.numeric(best_coords[["specificity"]])

      predicted_class <- gene_expression > cutoff

      confusion_matrix <- table(
        Predicted = predicted_class,
        True = class_vec
      )

      TP <- ifelse(
        "TRUE" %in% rownames(confusion_matrix) &&
          "TRUE" %in% colnames(confusion_matrix),
        confusion_matrix["TRUE", "TRUE"],
        0
      )

      TN <- ifelse(
        "FALSE" %in% rownames(confusion_matrix) &&
          "FALSE" %in% colnames(confusion_matrix),
        confusion_matrix["FALSE", "FALSE"],
        0
      )

      FP <- ifelse(
        "TRUE" %in% rownames(confusion_matrix) &&
          "FALSE" %in% colnames(confusion_matrix),
        confusion_matrix["TRUE", "FALSE"],
        0
      )

      FN <- ifelse(
        "FALSE" %in% rownames(confusion_matrix) &&
          "TRUE" %in% colnames(confusion_matrix),
        confusion_matrix["FALSE", "TRUE"],
        0
      )

      accuracy <- ifelse(
        (TP + TN + FP + FN) == 0,
        NA,
        (TP + TN) / (TP + TN + FP + FN)
      )

      precision <- ifelse(
        (TP + FP) == 0,
        NA,
        TP / (TP + FP)
      )
    }

    results <- rbind(
      results,
      data.frame(
        Gene = gene,
        AUC = round(auc_value, 3),
        Best_Cutoff = round(cutoff, 3),
        Sensitivity = round(sensitivity, 3),
        Specificity = round(specificity, 3),
        Accuracy = round(accuracy, 3),
        Precision = round(precision, 3),
        stringsAsFactors = FALSE
      )
    )
  }

  results <- results[order(results$AUC, decreasing = TRUE), ]

  return(results)
}

# -----------------------------
# Example usage
# -----------------------------
# expr_mat <- readRDS("path/to/demo_expression_matrix.rds")
# class_vec <- readRDS("path/to/demo_class_vector.rds")
# genes_to_test <- c("GENE1", "GENE2", "GENE3")
#
# roc_results <- evaluate_genes_by_roc(expr_mat, class_vec, genes_to_test)
#
# output_dir <- "results"
# dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# write.csv(roc_results, file.path(output_dir, "ROC_biomarker_results.csv"), row.names = FALSE)
