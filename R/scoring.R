#' Score macrophage biological programs
#'
#' Calculates weighted expression scores for the biological programs defined
#' in the Macrophage Framework ontology. Expression values are normalized using
#' the transcriptome-wide library size of each cell, scaled by
#' \code{scale_factor}, and log-transformed before program scores are
#' calculated.
#'
#' Only ontology genes detected in the supplied expression matrix are used for
#' scoring. Program scores are calculated as the weighted mean of the
#' log-transformed expression values of the genes belonging to each program.
#' Gene weights are defined by the \code{Weight} column of the ontology.
#'
#' Program coverage is calculated as the proportion of ontology genes available
#' in the expression matrix for each program. Coverage therefore describes
#' gene availability and is not used as a probability or confidence measure.
#'
#' @param expression A gene-by-cell expression matrix. Dense matrices and
#'   sparse matrices from the \pkg{Matrix} package are supported. Gene names
#'   must be provided as row names.
#' @param ontology Optional Framework ontology data frame. If \code{NULL},
#'   the default ontology is loaded with \code{\link{load_ontology}}.
#' @param scale_factor Numeric scaling factor used for library-size
#'   normalization. Defaults to \code{10000}.
#' @param library_size Optional vector containing the library size for each
#'   cell. If \code{NULL}, library sizes are calculated as the column sums of
#'   \code{expression}.
#'
#' @return A list containing:
#' \describe{
#'   \item{scores}{A program-by-cell matrix of weighted program scores.}
#'   \item{genes_available}{The number of ontology genes available in the
#'   expression matrix for each program.}
#'   \item{genes_total}{The total number of genes defined for each program in
#'   the ontology.}
#'   \item{coverage}{The proportion of ontology genes available for each
#'   program.}
#' }
#'
#' @examples
#' ontology <- load_ontology()
#'
#' set.seed(1)
#' expression <- matrix(
#'   rpois(100, lambda = 2),
#'   nrow = 10,
#'   dimnames = list(
#'     ontology$Gene[seq_len(10)],
#'     paste0("Cell", seq_len(10))
#'   )
#' )
#'
#' result <- score_programs(
#'   expression = expression,
#'   ontology = ontology
#' )
#'
#' @export
#'

score_programs <- function(expression,
                           ontology = NULL,
                           scale_factor = 10000,
                           library_size = NULL) {

  # Check input
  if (!is.matrix(expression) && !inherits(expression, "Matrix")) {
    stop("expression must be a matrix or sparse Matrix.")
  }

  if (is.null(rownames(expression))) {
    stop("expression must have gene names as rownames.")
  }

  if (is.null(ontology)) {
    ontology <- load_ontology()
  }

  required_cols <- c("Dimension", "Program", "Gene", "Weight")
  if (!all(required_cols %in% colnames(ontology))) {
    stop("ontology is missing required columns.")
  }

  # Keep only ontology genes present in the expression matrix
  available_genes <- intersect(ontology$Gene, rownames(expression))

  if (length(available_genes) == 0) {
    stop("No ontology genes found in expression matrix.")
  }

  ontology_available <- ontology[
    ontology$Gene %in% available_genes,
    ,
    drop = FALSE
  ]

  # Library size normalization using the full count matrix
  if (is.null(library_size)) {
    library_size <- Matrix::colSums(expression)
  } else {
    if (length(library_size) != ncol(expression)) {
      stop("library_size must have one value per cell.")
    }

    if (any(library_size <= 0)) {
      stop("All cells must have positive library size.")
    }
  }

  if (any(library_size <= 0)) {
    stop("All cells must have positive library size.")
  }

  # IMPORTANT:
  # Normalize using transcriptome-wide library size,
  # then retain ontology genes for program scoring.
  # The calculation remains memory-safe for sparse scRNA-seq matrices.

  expression_ontology <- expression[
    available_genes,
    ,
    drop = FALSE
  ]

  # Apply transcriptome-wide normalization and log transformation
  # to the ontology genes.

  normalized <- expression_ontology %*%
    Matrix::Diagonal(
      x = scale_factor / library_size
    )

  log_expression <- log1p(normalized)

  # Program names
  programs <- unique(ontology$Program)

  # Initialize score matrix
  scores <- matrix(
    NA_real_,
    nrow = length(programs),
    ncol = ncol(expression),
    dimnames = list(programs, colnames(expression))
  )

  genes_available <- numeric(length(programs))
  genes_total <- numeric(length(programs))

  names(genes_available) <- programs
  names(genes_total) <- programs

  # Calculate weighted program scores
  for (program in programs) {

    ontology_program <- ontology[
      ontology$Program == program,
      ,
      drop = FALSE
    ]

    genes_total[program] <- nrow(ontology_program)

    genes_program <- intersect(
      ontology_program$Gene,
      rownames(log_expression)
    )

    genes_available[program] <- length(genes_program)

    if (length(genes_program) == 0) {
      next
    }

    weights <- ontology_program$Weight[
      match(genes_program, ontology_program$Gene)
    ]

    values <- log_expression[
      genes_program,
      ,
      drop = FALSE
    ]

    scores[program, ] <- Matrix::colSums(
      values * weights
    ) / sum(weights)
  }

  coverage <- genes_available / genes_total

  list(
    scores = scores,
    genes_available = genes_available,
    genes_total = genes_total,
    coverage = coverage
  )
}

#' Calculate program detection coverage
#'
#' Calculates the fraction of genes belonging to each Macrophage Framework
#' program that are detected in each cell. A gene is considered detected when
#' its expression value is greater than zero.
#'
#' Coverage is calculated independently for each program and cell using only
#' ontology genes that are present in the supplied expression matrix. It
#' therefore describes the extent to which the genes defining a program are
#' detected in a cell and should be interpreted as a diagnostic measure rather
#' than as a probability, confidence score, or measure of program activity.
#'
#' @param expression A gene-by-cell expression matrix. Dense matrices and
#'   sparse matrices from the \pkg{Matrix} package are supported. Gene names
#'   must be provided as row names.
#' @param ontology Optional Framework ontology data frame. If \code{NULL},
#'   the default ontology is loaded with \code{\link{load_ontology}}.
#'
#' @return A program-by-cell numeric matrix containing the proportion of
#'   available ontology genes detected in each cell. Rows correspond to
#'   Framework programs and columns correspond to cells.
#'
#' @examples
#' ontology <- load_ontology()
#'
#' set.seed(1)
#' genes <- unique(ontology$Gene)
#' expression <- matrix(
#'   rpois(length(genes) * 5, lambda = 1),
#'   nrow = length(genes),
#'   dimnames = list(
#'     genes,
#'     paste0("Cell", seq_len(5))
#'   )
#' )
#'
#' coverage <- calculate_detection_coverage(
#'   expression = expression,
#'   ontology = ontology
#' )
#'
#' head(coverage)
#'
#' @export

calculate_detection_coverage <- function(expression,
                                         ontology = NULL) {

  # Check input
  if (!is.matrix(expression) &&
      !inherits(expression, "Matrix")) {
    stop("expression must be a matrix or sparse Matrix.")
  }

  if (is.null(rownames(expression))) {
    stop("expression must have gene names as row names.")
  }

  if (is.null(ontology)) {
    ontology <- load_ontology()
  }

  required_columns <- c("Program", "Gene")

  if (!all(required_columns %in% colnames(ontology))) {
    stop("Ontology must contain Program and Gene columns.")
  }

  # Keep only ontology genes present in expression matrix
  available_genes <- intersect(
    ontology$Gene,
    rownames(expression)
  )

  if (length(available_genes) == 0) {
    stop("No ontology genes found in expression matrix.")
  }

  expression_ontology <- expression[
    available_genes,
    ,
    drop = FALSE
  ]

  programs <- unique(ontology$Program)

  coverage <- matrix(
    NA_real_,
    nrow = length(programs),
    ncol = ncol(expression),
    dimnames = list(
      programs,
      colnames(expression)
    )
  )

  for (program in programs) {

    program_genes <- ontology$Gene[
      ontology$Program == program
    ]

    genes_present <- intersect(
      program_genes,
      rownames(expression_ontology)
    )

    if (length(genes_present) == 0) {
      next
    }

    detected <- expression_ontology[
      genes_present,
      ,
      drop = FALSE
    ] > 0

    coverage[program, ] <-
      Matrix::colSums(detected) /
      length(genes_present)
  }

  coverage
}

summarize_cluster_scores <- function(scores,
                                     cluster_labels) {

  # Check that scores is a matrix
  if (!is.matrix(scores)) {
    stop("scores must be a matrix.")
  }

  # Check that cluster labels are provided for every cell
  if (length(cluster_labels) != ncol(scores)) {
    stop("cluster_labels must contain one label for each cell.")
  }

  # Calculate mean score for each program within each cluster
  cluster_ids <- unique(cluster_labels)

  cluster_scores <- sapply(
    cluster_ids,
    function(cluster) {

      cells <- cluster_labels == cluster

      cluster_mean <- rowMeans(
        scores[, cells, drop = FALSE],
        na.rm = TRUE
      )

      cluster_mean[is.nan(cluster_mean)] <- NA_real_

      cluster_mean
    }
  )

  # Make sure the result is a matrix
  cluster_scores <- as.matrix(cluster_scores)

  rownames(cluster_scores) <- rownames(scores)
  colnames(cluster_scores) <- cluster_ids

  return(cluster_scores)
}

create_cluster_profile <- function(cluster_scores,
                                   cluster_coverage,
                                   cluster_coherence = NULL,
                                   cluster_coherence_IQR = NULL) {

  # Check scores
  if (!is.matrix(cluster_scores)) {
    stop("cluster_scores must be a matrix.")
  }

  # Check coverage
  if (!is.matrix(cluster_coverage)) {
    stop("cluster_coverage must be a matrix.")
  }

  # Scores and coverage must have same dimensions
  if (!identical(dim(cluster_scores), dim(cluster_coverage))) {
    stop("cluster_scores and cluster_coverage must have the same dimensions.")
  }

  # Program names must match
  if (!identical(
    rownames(cluster_scores),
    rownames(cluster_coverage)
  )) {
    stop("Program names must match.")
  }

  # Cluster names must match
  if (!identical(
    colnames(cluster_scores),
    colnames(cluster_coverage)
  )) {
    stop("Cluster names must match.")
  }

  # Start profile
  profile <- list(
    scores = cluster_scores,
    coverage = cluster_coverage
  )

  # Add coherence if provided
  if (!is.null(cluster_coherence)) {

    if (!is.matrix(cluster_coherence)) {
      stop("cluster_coherence must be a matrix.")
    }

    if (!identical(
      dim(cluster_scores),
      dim(cluster_coherence)
    )) {
      stop("cluster_coherence must have the same dimensions as cluster_scores.")
    }

    if (!identical(
      rownames(cluster_scores),
      rownames(cluster_coherence)
    )) {
      stop("Program names must match.")
    }

    if (!identical(
      colnames(cluster_scores),
      colnames(cluster_coherence)
    )) {
      stop("Cluster names must match.")
    }

    profile$coherence <- cluster_coherence
  }

  # Add coherence IQR if provided
  if (!is.null(cluster_coherence_IQR)) {

    if (!is.matrix(cluster_coherence_IQR)) {
      stop("cluster_coherence_IQR must be a matrix.")
    }

    if (!identical(
      dim(cluster_scores),
      dim(cluster_coherence_IQR)
    )) {
      stop("cluster_coherence_IQR must have the same dimensions as cluster_scores.")
    }

    if (!identical(
      rownames(cluster_scores),
      rownames(cluster_coherence_IQR)
    )) {
      stop("Program names must match.")
    }

    if (!identical(
      colnames(cluster_scores),
      colnames(cluster_coherence_IQR)
    )) {
      stop("Cluster names must match.")
    }

    profile$coherence_IQR <- cluster_coherence_IQR
  }

  return(profile)
}
