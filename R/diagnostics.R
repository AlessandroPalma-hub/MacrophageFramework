#' Calculate program coherence
#'
#' Quantifies the internal coherence of each Macrophage Framework program
#' across individual cells. Coherence is calculated from pairwise similarity
#' among the genes belonging to each program that are detected in a given
#' cell.
#'
#' Expression values are normalized using transcriptome-wide library size,
#' scaled by \code{scale_factor}, and log-transformed before coherence is
#' calculated. Genes with zero expression in a cell are treated as
#' unobserved and are excluded from the pairwise calculation.
#'
#' For each program and cell, pairwise relative expression differences are
#' calculated among observed program genes and combined using the product of
#' their ontology weights. The resulting coherence score is bounded between
#' 0 and 1, with higher values indicating greater internal agreement among
#' the observed genes of a program.
#'
#' Coherence is a diagnostic measure of internal program consistency and
#' should not be interpreted as a probability, confidence score, or direct
#' measure of program activity.
#'
#' @param expression A gene-by-cell expression matrix. Dense matrices and
#'   sparse matrices from the \pkg{Matrix} package are supported. Gene names
#'   must be provided as row names.
#' @param ontology Optional Framework ontology data frame. If \code{NULL},
#'   the default ontology is loaded with \code{\link{load_ontology}}.
#'   The ontology must contain \code{Program}, \code{Gene}, and \code{Weight}
#'   columns.
#' @param scale_factor Numeric scaling factor used for library-size
#'   normalization. Defaults to \code{10000}.
#'
#' @return A program-by-cell numeric matrix containing coherence scores.
#'   Rows correspond to Framework programs and columns correspond to cells.
#'   Values range from 0 to 1. \code{NA} values indicate that fewer than two
#'   program genes were detected in the corresponding cell.
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
#' coherence <- calculate_program_coherence(
#'   expression = expression,
#'   ontology = ontology
#' )
#'
#' head(coherence)
#'
#' @export
calculate_program_coherence <- function(expression,
                                        ontology = NULL,
                                        scale_factor = 10000) {

  if (!is.matrix(expression) && !inherits(expression, "Matrix")) {
    stop("expression must be a matrix or sparse Matrix.")
  }

  if (is.null(rownames(expression))) {
    stop("expression must have gene names as rownames.")
  }

  if (is.null(ontology)) {
    ontology <- load_ontology()
  }

  required_cols <- c("Program", "Gene", "Weight")

  if (!all(required_cols %in% colnames(ontology))) {
    stop("ontology is missing required columns.")
  }

  available_genes <- intersect(
    ontology$Gene,
    rownames(expression)
  )

  if (length(available_genes) == 0) {
    stop("No ontology genes found in expression matrix.")
  }

  # Work only with ontology genes
  expression_ontology <- expression[
    available_genes,
    ,
    drop = FALSE
  ]

  # Library-size normalization
  library_size <- Matrix::colSums(expression)

  if (any(library_size <= 0)) {
    stop("All cells must have positive library size.")
  }

  normalized <- expression_ontology %*%
    Matrix::Diagonal(
      x = scale_factor / library_size
    )

  log_expression <- log1p(normalized)

  programs <- unique(ontology$Program)

  coherence <- matrix(
    NA_real_,
    nrow = length(programs),
    ncol = ncol(expression),
    dimnames = list(
      programs,
      colnames(expression)
    )
  )

  for (program in programs) {

    ontology_program <- ontology[
      ontology$Program == program,
      ,
      drop = FALSE
    ]

    genes_program <- intersect(
      ontology_program$Gene,
      rownames(log_expression)
    )

    if (length(genes_program) < 2) {
      next
    }

    weights <- ontology_program$Weight[
      match(
        genes_program,
        ontology_program$Gene
      )
    ]

    values <- as.matrix(
      log_expression[
        genes_program,
        ,
        drop = FALSE
      ]
    )

    for (cell in seq_len(ncol(values))) {

      x <- values[, cell]

      # Zero expression is treated as unobserved/dropout
      observed <- x > 0

      if (sum(observed) < 2) {
        next
      }

      x_obs <- x[observed]
      w_obs <- weights[observed]

      n <- length(x_obs)

      pair_distance <- numeric(0)
      pair_weight <- numeric(0)

      for (i in seq_len(n - 1)) {

        for (j in (i + 1):n) {

          distance <- abs(
            x_obs[i] - x_obs[j]
          ) / (
            x_obs[i] + x_obs[j] + 1e-8
          )

          weight <- w_obs[i] * w_obs[j]

          pair_distance <- c(
            pair_distance,
            distance
          )

          pair_weight <- c(
            pair_weight,
            weight
          )
        }
      }

      coherence[program, cell] <-
        1 - sum(pair_distance * pair_weight) /
        sum(pair_weight)
    }
  }

  # Numerical safety
  coherence[coherence < 0] <- 0
  coherence[coherence > 1] <- 1

  coherence
}

calculate_cell_coherence_diagnostics <- function(expression,
                                                 ontology,
                                                 clusters,
                                                 scale_factor = 10000) {

  if (!is.matrix(expression) && !inherits(expression, "Matrix")) {
    stop("expression must be a matrix or sparse Matrix.")
  }

  clusters <- as.character(clusters)

  genes_present <- intersect(
    ontology$Gene,
    rownames(expression)
  )

  expression_ontology <- expression[
    genes_present, ,
    drop = FALSE
  ]

  library_size <- Matrix::colSums(expression)
  library_size[library_size == 0] <- 1

  normalized <- expression_ontology %*%
    Matrix::Diagonal(
      x = scale_factor / library_size
    )

  log_expression <- log1p(normalized)

  results <- list()

  for (program in unique(ontology$Program)) {

    program_info <- ontology[
      ontology$Program == program &
        ontology$Gene %in% genes_present,
      ,
      drop = FALSE
    ]

    program_genes <- program_info$Gene
    weights <- program_info$Weight
    names(weights) <- program_genes

    if (length(program_genes) < 2) next

    for (cl in unique(clusters)) {

      cells <- which(clusters == cl)

      cell_coherence <- numeric(length(cells))
      n_observed <- integer(length(cells))

      for (i in seq_along(cells)) {

        cell <- cells[i]

        values <- as.numeric(
          log_expression[program_genes, cell]
        )

        observed <- values > 0
        n_observed[i] <- sum(observed)

        if (n_observed[i] < 2) {
          cell_coherence[i] <- NA_real_
          next
        }

        obs_values <- values[observed]
        obs_weights <- weights[observed]

        weighted_distance <- 0
        total_weight <- 0

        for (a in seq_len(length(obs_values) - 1)) {

          for (b in (a + 1):length(obs_values)) {

            pair_weight <- obs_weights[a] * obs_weights[b]

            relative_distance <-
              abs(obs_values[a] - obs_values[b]) /
              (obs_values[a] + obs_values[b] + 1e-8)

            weighted_distance <-
              weighted_distance +
              pair_weight * relative_distance

            total_weight <-
              total_weight + pair_weight
          }
        }

        cell_coherence[i] <-
          1 - weighted_distance / total_weight
      }

      eligible <- !is.na(cell_coherence)

      results[[length(results) + 1]] <- data.frame(
        Program = program,
        Cluster = cl,
        MedianCoherence = median(
          cell_coherence[eligible]
        ),
        IQRCoherence = IQR(
          cell_coherence[eligible]
        ),
        MedianGenesObserved = median(
          n_observed[eligible]
        ),
        IQRGenesObserved = IQR(
          n_observed[eligible]
        ),
        CoherenceEligibleFraction =
          mean(eligible),
        CellsWithCoherence =
          sum(eligible),
        TotalCells =
          length(cells),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, results)
}

#' Summarize program coherence by cluster
#'
#' @param cell_coherence Program-by-cell coherence matrix.
#' @param cluster_labels Cluster assignment for each cell.
#'
#' @return A list containing cluster median and IQR coherence.
summarize_cluster_coherence <- function(cell_coherence,
                                        cluster_labels) {

  if (!is.matrix(cell_coherence)) {
    stop("cell_coherence must be a matrix.")
  }

  if (length(cluster_labels) != ncol(cell_coherence)) {
    stop(
      "cluster_labels must contain one label for each cell."
    )
  }

  clusters <- unique(as.character(cluster_labels))

  median_coherence <- matrix(
    NA_real_,
    nrow = nrow(cell_coherence),
    ncol = length(clusters),
    dimnames = list(
      rownames(cell_coherence),
      clusters
    )
  )

  IQR_coherence <- median_coherence

  for (cluster in clusters) {

    cells <- which(
      as.character(cluster_labels) == cluster
    )

    for (program in rownames(cell_coherence)) {

      values <- cell_coherence[
        program,
        cells
      ]

      values <- values[
        is.finite(values)
      ]

      if (length(values) == 0) {
        next
      }

      median_coherence[program, cluster] <-
        median(values)

      IQR_coherence[program, cluster] <-
        IQR(values)
    }
  }

  list(
    median = median_coherence,
    IQR = IQR_coherence
  )
}
