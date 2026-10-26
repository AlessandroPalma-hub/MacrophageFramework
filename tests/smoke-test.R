
library(MacrophageFramework)

# Check that the default ontology can be loaded.
ontology <- load_ontology("v0.4")

stopifnot(
  is.data.frame(ontology),
  all(c("Gene", "Program", "Dimension") %in% names(ontology))
)

# Build a small synthetic expression matrix.
genes <- unique(
  ontology$Gene[
    !is.na(ontology$Gene) & nzchar(ontology$Gene)
  ]
)

set.seed(1)

expression <- matrix(
  sample(1:5, length(genes) * 6, replace = TRUE),
  nrow = length(genes),
  dimnames = list(
    genes,
    paste0("Cell", seq_len(6))
  )
)

clusters <- rep(c("Cluster1", "Cluster2"), each = 3)

# Check that the main analysis runs.
result <- analyze_macrophage_profile(
  expression = expression,
  clusters = clusters,
  ontology = ontology,
  detect_macrophages = FALSE,
  macrophage_clusters = unique(clusters),
  species = "Human",
  save_outputs = FALSE
)

stopifnot(
  is.list(result),
  all(paste0("Cell", seq_len(6)) %in% result$analyzed_cells)
)

# Check that the results can be summarized.
summary <- summarize_macrophage_profile(
  result = result,
  clusters = clusters
)

stopifnot(
  is.data.frame(summary),
  nrow(summary) > 0,
  all(c("Cluster", "Cells") %in% names(summary))
)
