# 11_PPI_network_centrality.R
# ------------------------------------------------------------
# Purpose:
# This script demonstrates basic protein-protein interaction
# network analysis using igraph. It calculates common centrality
# metrics such as degree, betweenness, eigenvector centrality,
# and closeness.
#
# Notes:
# - The input edge file should contain at least two columns:
#   `node1` and `node2`.
# - The edge file comes from resources such as STRING or
#   Cytoscape-exported interaction tables.
# ------------------------------------------------------------

library(igraph)

# -----------------------------
# User-defined input placeholder
# -----------------------------
# edge_file <- "path/to/demo_ppi_edges.tsv"
# output_dir <- "results"

edges <- read.delim(
  edge_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

required_columns <- c("node1", "node2")

if (!all(required_columns %in% colnames(edges))) {
  stop("The edge file must contain columns named 'node1' and 'node2'.")
}

# -----------------------------
# Build graph
# -----------------------------
ppi_graph <- graph_from_data_frame(
  d = edges[, c("node1", "node2")],
  directed = FALSE
)

# -----------------------------
# Calculate centrality metrics
# -----------------------------
centrality_table <- data.frame(
  Gene = V(ppi_graph)$name,
  Betweenness = betweenness(ppi_graph, normalized = TRUE),
  Degree = degree(ppi_graph, mode = "all"),
  Eigenvector = eigen_centrality(ppi_graph)$vector,
  Closeness = closeness(ppi_graph, normalized = TRUE)
)

centrality_table$Betweenness_Z <- as.numeric(scale(centrality_table$Betweenness))
centrality_table$Degree_Z <- as.numeric(scale(centrality_table$Degree))
centrality_table$Eigenvector_Z <- as.numeric(scale(centrality_table$Eigenvector))
centrality_table$Closeness_Z <- as.numeric(scale(centrality_table$Closeness))

centrality_table <- centrality_table[order(centrality_table$Degree, decreasing = TRUE), ]

# -----------------------------
# Save centrality table
# -----------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  centrality_table,
  file = file.path(output_dir, "PPI_network_centrality_metrics.csv"),
  row.names = FALSE
)
