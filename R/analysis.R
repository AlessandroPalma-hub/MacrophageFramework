#' Analyze a macrophage expression profile
#'
#' Runs the independent macrophage detector and applies the Macrophage
#' Framework to clusters passing the detector gate. The detector establishes
#' macrophage identity independently from the Framework ontology, whereas the
#' Framework characterizes lineage-associated and functional-state programs.
#'
#' The function accepts either a gene-by-cell expression matrix together with
#' cluster assignments or a Seurat object. For Seurat objects, expression data
#' and cluster identities are extracted automatically. The Framework ontology
#' is version 0.4 by default.
#'
#' @param expression A gene-by-cell expression matrix. Ignored when
#'   \code{object} is supplied.
#' @param clusters A vector assigning each cell to a cluster. Ignored when
#'   \code{object} is supplied.
#' @param ontology Optional Framework ontology. If \code{NULL}, the default
#'   ontology is loaded with \code{\link{load_ontology}}.
#' @param object Optional Seurat object. If supplied, expression data and
#'   cluster identities are extracted automatically.
#' @param assay Assay to use when \code{object} is a Seurat object.
#' @param layer Expression layer to use when \code{object} is a Seurat object.
#'   If \code{NULL}, available count or data layers are selected automatically.
#' @param reduction Reduction used to obtain the embedding from a Seurat
#'   object. If \code{NULL}, an available embedding is selected automatically.
#' @param detect_macrophages Logical indicating whether the independent
#'   macrophage detector should be run.
#' @param macrophage_clusters Optional vector specifying clusters to analyze
#'   directly. If \code{NULL}, clusters are selected using the macrophage
#'   detector.
#' @param min_core_score Minimum cluster-level macrophage-core score required
#'   for inclusion when macrophage clusters are selected by the detector.
#' @param species Species used for ontology gene mapping. Supported values
#'   include \code{"Human"} and \code{"Mouse"}.
#' @param output_dir Directory in which analysis outputs are saved when
#'   \code{save_outputs = TRUE}.
#' @param save_outputs Logical indicating whether analysis results and
#'   visualization outputs should be written to \code{output_dir}.
#'
#' @return A list containing macrophage detection results, Framework profiles,
#'   similarity, hierarchy, program scores, coverage and coherence diagnostics,
#'   ontology, annotations, embedding information, and the analyzed cells and
#'   clusters. When \code{save_outputs = TRUE}, the list also contains paths
#'   to the generated output files.
#'
#' @examples
#' \dontrun{
#' result <- analyze_macrophage_profile(
#'   object = seurat_object,
#'   species = "Human"
#' )
#' }
#'
#' @export
analyze_macrophage_profile <- function(
    expression = NULL,
    clusters = NULL,
    ontology = NULL,
    object = NULL,
    assay = "RNA",
    layer = NULL,
    reduction = NULL,
    detect_macrophages = TRUE,
    macrophage_clusters = NULL,
    min_core_score = 0.1,
    species = "Human",
    output_dir = "MacrophageFramework_output",
    save_outputs = TRUE
) {

  # ----------------------------------------------------------
  # Input mode
  # ----------------------------------------------------------

  if (!is.null(object)) {

    if (!inherits(object, "Seurat")) {
      stop("object must be a Seurat object.")
    }

    expression <- .extract_seurat_expression(
      object = object,
      assay = assay,
      layer = layer
    )

    clusters <- .extract_seurat_clusters(
      object = object,
      cells = colnames(expression)
    )

    embedding <- .get_embedding(
      object = object,
      reduction = reduction
    )

  } else {

    embedding <- NULL

    if (is.null(expression) || is.null(clusters)) {
      stop(
        "Provide either 'object' or both 'expression' and 'clusters'."
      )
    }
  }

  # ----------------------------------------------------------
  # Basic checks
  # ----------------------------------------------------------

  if (!is.matrix(expression) &&
      !inherits(expression, "Matrix")) {
    stop(
      "expression must be a matrix or sparse Matrix."
    )
  }

  if (is.null(rownames(expression))) {
    stop(
      "expression must have gene names as rownames."
    )
  }

  if (is.null(colnames(expression))) {
    stop(
      "expression must have cell names as colnames."
    )
  }

  if (length(clusters) != ncol(expression)) {
    stop(
      "clusters must contain one label for each cell."
    )
  }

  if (any(is.na(clusters))) {
    stop(
      "clusters must not contain NA values."
    )
  }

  clusters <- as.character(clusters)

  # ----------------------------------------------------------
  # Framework ontology
  # ----------------------------------------------------------

  if (is.null(ontology)) {
    ontology <- load_ontology("v0.4")
  }

  required_columns <- c(
    "Dimension",
    "Program",
    "Gene",
    "Weight"
  )

  if (!all(required_columns %in% colnames(ontology))) {
    stop(
      "Ontology is missing required columns."
    )
  }

  # Identity must not be part of the current Framework ontology.
  if ("Identity" %in% unique(ontology$Dimension)) {
    stop(
      "Current Framework ontology must not contain Identity. ",
      "Macrophage identity is established by the independent detector."
    )
  }

  # ----------------------------------------------------------
  # Independent macrophage detector
  # ----------------------------------------------------------

  if (isTRUE(detect_macrophages)) {

    macrophage_detection <- .detect_macrophage_core(
      expression = expression,
      clusters = clusters,
      min_core_score = min_core_score,
      species = species
    )

    included_clusters <-
      macrophage_detection$included_clusters

  } else {

    if (is.null(macrophage_clusters)) {
      stop(
        "When detect_macrophages = FALSE, ",
        "macrophage_clusters must be provided."
      )
    }

    macrophage_detection <- NULL

    included_clusters <-
      intersect(
        as.character(macrophage_clusters),
        unique(as.character(clusters))
      )

    if (length(included_clusters) == 0) {
      stop(
        "None of the provided macrophage_clusters ",
        "were found in the cluster labels."
      )
    }

    # Create a minimal detector-like table for downstream outputs
    detector_table <- data.frame(
      Cluster = sort(unique(as.character(clusters))),
      Status = ifelse(
        sort(unique(as.character(clusters))) %in% included_clusters,
        "included",
        "excluded"
      ),
      stringsAsFactors = FALSE
    )

  }

  # ----------------------------------------------------------
  # Stop cleanly if no macrophage clusters pass the gate
  # ----------------------------------------------------------

  if (length(included_clusters) == 0) {

    result <- list(
      macrophage_detection = macrophage_detection,
      profile = NULL,
      similarity = NULL,
      hierarchy = NULL,
      scores = NULL,
      detection = NULL,
      coherence = NULL,
      ontology = ontology,
      annotations = NULL,
      embedding = embedding,
      analyzed_cells = character(0),
      analyzed_clusters = character(0),
      output_dir = output_dir,
      saved_files = character(0)
    )

    if (isTRUE(save_outputs)) {
      result$saved_files <- .save_framework_outputs(
        result = result,
        embedding = embedding,
        clusters = clusters,
        output_dir = output_dir
      )
    }

    return(result)
  }

  # ----------------------------------------------------------
  # Restrict Framework analysis to included clusters
  # ----------------------------------------------------------

  keep_cells <- clusters %in% included_clusters

  expression_macrophages <- expression[
    ,
    keep_cells,
    drop = FALSE
  ]

  clusters_macrophages <- clusters[
    keep_cells
  ]

  # ----------------------------------------------------------
  # Run Framework
  # ----------------------------------------------------------

  framework <- .analyze_expression_profile(
    expression = expression_macrophages,
    clusters = clusters_macrophages,
    ontology = ontology,
    species = species
  )

  # ----------------------------------------------------------
  # Annotation
  # ----------------------------------------------------------

  annotations <- .generate_framework_annotations(
    profile = framework$profile,
    ontology = ontology
  )

  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------

  result <- list(
    macrophage_detection = macrophage_detection,
    profile = framework$profile,
    similarity = framework$similarity,
    hierarchy = framework$hierarchy,
    scores = framework$scores,
    detection = framework$detection,
    coherence = framework$coherence,
    ontology = ontology,
    annotations = annotations,
    embedding = embedding,
    analyzed_cells = colnames(expression_macrophages),
    analyzed_clusters = included_clusters,
    output_dir = output_dir,
    saved_files = character(0)
  )

  if (isTRUE(save_outputs)) {
    result$saved_files <- .save_framework_outputs(
      result = result,
      embedding = embedding,
      clusters = clusters,
      output_dir = output_dir
    )
  }

  return(result)
}

.extract_seurat_expression <- function(
    object,
    assay = "RNA",
    layer = NULL
) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required for Seurat input.")
  }

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }

  if (!assay %in% names(object@assays)) {
    stop("Assay '", assay, "' was not found in the Seurat object.")
  }

  assay_object <- object[[assay]]

  available_layers <- Layers(assay_object)

  if (is.null(layer)) {

    data_layers <- available_layers[
      grepl("^data", available_layers)
    ]

    count_layers <- available_layers[
      grepl("^counts", available_layers)
    ]

    if (length(count_layers) > 0) {
      layer_names <- count_layers
    } else if (length(data_layers) > 0) {
      layer_names <- data_layers
    } else {
      stop("No data or counts layer found in assay.")
    }

  } else {

    if (!layer %in% available_layers) {
      stop(
        "Layer '", layer,
        "' was not found in assay '", assay, "'."
      )
    }

    layer_names <- layer
  }

  matrices <- lapply(
    layer_names,
    function(x) {
      LayerData(
        object = assay_object,
        layer = x
      )
    }
  )

  if (length(matrices) == 1) {
    expression <- matrices[[1]]
  } else {

    common_genes <- Reduce(
      intersect,
      lapply(matrices, rownames)
    )

    if (length(common_genes) == 0) {
      stop("No common genes found across expression layers.")
    }

    matrices <- lapply(
      matrices,
      function(x) x[common_genes, , drop = FALSE]
    )

    expression <- do.call(
      cbind,
      matrices
    )
  }

  common_cells <- intersect(
    colnames(expression),
    colnames(object)
  )

  if (length(common_cells) == 0) {
    stop("No common cells between expression layer and Seurat object.")
  }

  expression <- expression[
    ,
    common_cells,
    drop = FALSE
  ]

  return(expression)
}


# ============================================================
# Internal: extract Seurat clusters
# ============================================================

.extract_seurat_clusters <- function(
    object,
    cells
) {

  clusters <- Seurat::Idents(object)

  if (is.null(names(clusters))) {
    stop("Seurat identities do not have cell names.")
  }

  clusters <- clusters[cells]

  if (any(is.na(clusters))) {
    stop("Some analyzed cells have missing cluster identities.")
  }

  return(as.character(clusters))
}

# ============================================================
# Internal: core Framework analysis
# ============================================================

.analyze_expression_profile <- function(
    expression,
    clusters,
    ontology,
    species = "Human"
) {

  # ----------------------------------------------------------
  # Cell-level scores
  # ----------------------------------------------------------
  gene_mapping <- .resolve_ontology_genes(
    ontology = ontology,
    expression_genes = rownames(expression),
    species = species
  )

  ontology <- gene_mapping$ontology

  scores <- score_programs(
    expression = expression,
    ontology = ontology
  )

  detection <- calculate_detection_coverage(
    expression = expression,
    ontology = ontology
  )

  # ----------------------------------------------------------
  # Cluster-level scores
  # ----------------------------------------------------------

  cluster_scores <- summarize_cluster_scores(
    scores$scores,
    clusters
  )

  cluster_coverage <- summarize_cluster_scores(
    detection,
    clusters
  )

  # ----------------------------------------------------------
  # Coherence
  # ----------------------------------------------------------

  cell_coherence <- calculate_program_coherence(
    expression = expression,
    ontology = ontology
  )

  cluster_coherence <- summarize_cluster_coherence(
    cell_coherence,
    clusters
  )

  # ----------------------------------------------------------
  # Multidimensional profile
  # ----------------------------------------------------------

  profile <- create_cluster_profile(
    cluster_scores = cluster_scores,
    cluster_coverage = cluster_coverage,
    cluster_coherence = cluster_coherence$median,
    cluster_coherence_IQR = cluster_coherence$IQR
  )

  # ----------------------------------------------------------
  # Similarity
  # ----------------------------------------------------------

  similarity <- calculate_similarity_matrix(
    profile
  )

  # ----------------------------------------------------------
  # Hierarchical ranking
  # ----------------------------------------------------------

  hierarchy <- rank_programs_by_dimension(
    profile,
    ontology = ontology
  )

  return(
    list(
      profile = profile,
      similarity = similarity,
      hierarchy = hierarchy,
      scores = scores$scores,
      detection = detection,
      coherence = cell_coherence,
      gene_resolution = gene_mapping$diagnostics
    )
  )
}

