# ============================================================
# Macrophage identity detection
# ============================================================

# Independent macrophage-core program.
#
# IMPORTANT:
# This program is intentionally kept outside the main Framework
# ontology. It is used only as an upstream gate to determine
# whether a population provides sufficient evidence of
# macrophage identity before Framework profiling.

.macrophage_core <- data.frame(
  Gene = c(
    "CD68",
    "CSF1R",
    "FCGR1A",
    "ADGRE1",
    "CD14"
  ),
  Weight = c(
    2,
    3,
    3,
    2,
    1
  ),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Internal helper: log-normalization
# ------------------------------------------------------------

.normalize_expression_for_detection <- function(
    expression,
    scale_factor = 10000
) {

  library_size <- Matrix::colSums(expression)

  if (any(library_size <= 0)) {
    stop("All cells must have positive library size.")
  }

  genes <- rownames(expression)

  expression_norm <- t(
    t(expression) / library_size
  ) * scale_factor

  expression_norm <- log1p(expression_norm)

  return(expression_norm)
}

# ------------------------------------------------------------
# Cluster-level macrophage detection
# ------------------------------------------------------------


.get_macrophage_core_ontology <- function() {

  data.frame(
    index = seq_len(5),
    Dimension = "Detector",
    Program = "Macrophage core",
    Gene = c(
      "CD68",
      "CSF1R",
      "FCGR1A",
      "ADGRE1",
      "CD14"
    ),
    Role = "Core",
    Weight = c(
      2,
      3,
      3,
      2,
      1
    ),
    Program.specificity = "High",
    Context.dependence = "Medium",
    Cross.dimension = "No",
    Evidence.type = "Framework detector",
    Species = "Human;Mouse",
    Source = "MacrophageFramework",
    Notes = "Independent macrophage identity detector",
    stringsAsFactors = FALSE
  )
}

.detect_macrophage_core <- function(
    expression,
    clusters,
    min_core_score = 0.1,
    species = "Human"
) {

  detector_ontology <- .get_macrophage_core_ontology()

  gene_mapping <- .resolve_ontology_genes(
    ontology = detector_ontology,
    expression_genes = rownames(expression),
    species = species
  )

  detector_ontology <- gene_mapping$ontology

  detector_scores <- score_programs(
    expression = expression,
    ontology = detector_ontology
  )

  detector_coverage <- calculate_detection_coverage(
    expression = expression,
    ontology = detector_ontology
  )

  detector_coherence <- calculate_program_coherence(
    expression = expression,
    ontology = detector_ontology
  )

  core_score <- as.numeric(
    detector_scores$scores[
      "Macrophage core",
    ]
  )

  core_coverage <- as.numeric(
    detector_coverage[
      "Macrophage core",
    ]
  )

  core_coherence <- as.numeric(
    detector_coherence[
      "Macrophage core",
    ]
  )

  cell_metrics <- data.frame(
    Cell = colnames(expression),
    macrophage_core_score = core_score,
    macrophage_core_coverage = core_coverage,
    macrophage_core_coherence = core_coherence,
    macrophage_detectable =
      !is.na(core_score) &
      core_score >= min_core_score,
    stringsAsFactors = FALSE
  )

  cluster_ids <- unique(as.character(clusters))

  cluster_metrics <- do.call(
    rbind,
    lapply(
      cluster_ids,
      function(cluster) {

        cells <- as.character(clusters) == cluster

        scores_cluster <- core_score[cells]
        coverage_cluster <- core_coverage[cells]
        coherence_cluster <- core_coherence[cells]
        detectable_cluster <-
          cell_metrics$macrophage_detectable[cells]

        data.frame(
          Cluster = cluster,
          n_cells = sum(cells),
          macrophage_core_score =
            mean(scores_cluster, na.rm = TRUE),
          macrophage_core_coverage =
            mean(coverage_cluster, na.rm = TRUE),
          macrophage_core_coherence =
            median(coherence_cluster, na.rm = TRUE),
          pct_detectable =
            mean(detectable_cluster, na.rm = TRUE) * 100,
          Status =
            ifelse(
              mean(scores_cluster, na.rm = TRUE) >=
                min_core_score,
              "included",
              "excluded"
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  included_clusters <- cluster_metrics$Cluster[
    cluster_metrics$Status == "included"
  ]

  excluded_clusters <- cluster_metrics$Cluster[
    cluster_metrics$Status == "excluded"
  ]

  return(
    list(
      cluster_metrics = cluster_metrics,
      included_clusters = included_clusters,
      excluded_clusters = excluded_clusters,
      cell_metrics = cell_metrics,
      threshold = min_core_score,
      ontology = detector_ontology,
      gene_resolution = gene_mapping$diagnostics
    )
  )
}
