#' Retrieve a two-dimensional embedding
#'
#' Retrieves a precomputed two-dimensional embedding from a Seurat object.
#' If no reduction is specified, UMAP, t-SNE, or the first two PCA
#' dimensions are searched in this order.
#'
#' @param object A Seurat object.
#' @param reduction Optional name of the dimensional reduction to use.
#'
#' @return A two-column matrix containing the embedding coordinates.
#' @keywords internal
.get_embedding <- function(object, reduction = NULL) {

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }

  available_reductions <- names(object@reductions)

  # User-specified reduction
  if (!is.null(reduction)) {

    if (!reduction %in% available_reductions) {
      stop(
        "Reduction '",
        reduction,
        "' was not found in the Seurat object. ",
        "Available reductions: ",
        paste(available_reductions, collapse = ", "),
        "."
      )
    }

    embedding <- Seurat::Embeddings(
      object,
      reduction = reduction
    )

    if (ncol(embedding) < 2) {
      stop(
        "Reduction '",
        reduction,
        "' does not contain at least two dimensions."
      )
    }

    embedding <- embedding[, 1:2, drop = FALSE]

    return(embedding)
  }

  # Automatic search: UMAP > t-SNE > PCA
  reduction_priority <- c(
    "umap",
    "tsne",
    "pca"
  )

  selected_reduction <- reduction_priority[
    reduction_priority %in% available_reductions
  ]

  if (length(selected_reduction) == 0) {
    stop(
      "No suitable dimensional reduction was found. ",
      "Please provide a Seurat object containing UMAP, t-SNE, or PCA, ",
      "or specify a reduction explicitly."
    )
  }

  selected_reduction <- selected_reduction[1]

  embedding <- Seurat::Embeddings(
    object,
    reduction = selected_reduction
  )

  if (ncol(embedding) < 2) {
    stop(
      "Reduction '",
      selected_reduction,
      "' does not contain at least two dimensions."
    )
  }

  embedding <- embedding[, 1:2, drop = FALSE]

  return(embedding)
}
