#' Summarize macrophage profiles in a human-readable table
#'
#' @param result Output from analyze_macrophage_profile().
#' @param clusters Optional vector of cell-level cluster assignments.
#'
#' @return A data.frame summarizing the multidimensional macrophage profiles.
#' @export
summarize_macrophage_profile <- function(
    result,
    clusters = NULL
) {

  # --------------------------------------------------
  # Input checks
  # --------------------------------------------------

  if (!is.list(result)) {
    stop(
      "result must be the output of analyze_macrophage_profile()."
    )
  }

  if (is.null(result$profile) ||
      is.null(result$profile$scores) ||
      is.null(result$profile$coverage) ||
      is.null(result$profile$coherence)) {

    stop(
      "result does not contain a complete macrophage profile."
    )
  }

  scores <- result$profile$scores
  coverage <- result$profile$coverage
  coherence <- result$profile$coherence

  cluster_ids <- colnames(scores)

  # --------------------------------------------------
  # Cluster sizes
  # --------------------------------------------------

  cluster_size <- rep(
    NA_integer_,
    length(cluster_ids)
  )

  names(cluster_size) <- cluster_ids

  if (!is.null(clusters)) {

    clusters <- as.character(clusters)

    cluster_size <- table(clusters)

    cluster_size <- cluster_size[cluster_ids]

    cluster_size <- as.integer(cluster_size)
  }

  # --------------------------------------------------
  # Hierarchy
  # --------------------------------------------------

  hierarchy <- result$hierarchy

  if (is.null(hierarchy$Lineage) ||
      is.null(hierarchy$`Functional State`)) {

    stop(
      "result$hierarchy must contain 'Lineage' and 'Functional State'."
    )
  }

  lineage <- hierarchy$Lineage

  functional <- hierarchy$`Functional State`

  # --------------------------------------------------
  # Build summary
  # --------------------------------------------------

  summary <- data.frame(
    Cluster = cluster_ids,
    Cells = cluster_size,
    Lineage = NA_character_,
    FunctionalState = NA_character_,
    ProfileLabel = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(cluster_ids)) {

    cl <- cluster_ids[i]

    lineage_rank <- lineage[[cl]]
    functional_rank <- functional[[cl]]

    summary$Lineage[i] <-
      if (length(lineage_rank) > 0) {
        lineage_rank[1]
      } else {
        NA_character_
      }

    summary$FunctionalState[i] <-
      if (length(functional_rank) > 0) {
        functional_rank[1]
      } else {
        NA_character_
      }

    summary$ProfileLabel[i] <- paste(
      summary$Lineage[i],
      summary$FunctionalState[i],
      sep = " / "
    )
  }

  # --------------------------------------------------
  # Add program intensities
  # --------------------------------------------------

  for (program in rownames(scores)) {

    column_name <- paste0(
      "Score_",
      gsub("[^A-Za-z0-9]+", "_", program)
    )

    summary[[column_name]] <-
      as.numeric(
        scores[program, cluster_ids]
      )
  }

  # --------------------------------------------------
  # Add coverage
  # --------------------------------------------------

  for (program in rownames(coverage)) {

    column_name <- paste0(
      "Coverage_",
      gsub("[^A-Za-z0-9]+", "_", program)
    )

    summary[[column_name]] <-
      as.numeric(
        coverage[program, cluster_ids]
      )
  }

  # --------------------------------------------------
  # Add coherence
  # --------------------------------------------------

  for (program in rownames(coherence)) {

    column_name <- paste0(
      "Coherence_",
      gsub("[^A-Za-z0-9]+", "_", program)
    )

    summary[[column_name]] <-
      as.numeric(
        coherence[program, cluster_ids]
      )
  }

  # --------------------------------------------------
  # Add dominant functional program
  # --------------------------------------------------

  summary$DominantFunctionalProgram <-
    vapply(
      cluster_ids,
      function(cl) {

        functional_rank <- functional[[cl]]

        if (length(functional_rank) > 0) {
          functional_rank[1]
        } else {
          NA_character_
        }
      },
      character(1)
    )

  return(summary)
}

#' Plot multidimensional macrophage profiles
#'
#' Visualizes the relative intensity of Macrophage Framework programs across
#' analyzed macrophage clusters. Program scores are normalized independently
#' for each cluster to unit Euclidean norm for visualization, allowing the
#' relative composition of multidimensional profiles to be compared.
#'
#' This normalization is applied only for visualization and does not modify
#' the original program scores stored in the analysis result.
#'
#' @param result Output list returned by
#'   \code{\link{analyze_macrophage_profile}}.
#'
#' @return A \code{ggplot} object showing normalized program-intensity
#'   profiles for the analyzed clusters.
#'
#' @examples
#' \dontrun{
#' result <- analyze_macrophage_profile(
#'   expression = expression_matrix,
#'   clusters = cluster_labels
#' )
#'
#' plot_macrophage_profile(result)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_classic theme element_text
#' @export
plot_macrophage_profile <- function(result) {

  if (!is.list(result)) {
    stop("result must be the output of analyze_macrophage_profile().")
  }

  if (is.null(result$profile$scores)) {
    stop("result does not contain program scores.")
  }

  scores <- result$profile$scores

  # --------------------------------------------------
  # Normalize profiles for visualization only
  # --------------------------------------------------

  normalized_scores <- apply(
    scores,
    2,
    function(x) {

      denominator <- sqrt(sum(x^2, na.rm = TRUE))

      if (denominator == 0) {
        return(rep(NA_real_, length(x)))
      }

      x / denominator
    }
  )

  rownames(normalized_scores) <- rownames(scores)
  colnames(normalized_scores) <- colnames(scores)

  # --------------------------------------------------
  # Convert to long format
  # --------------------------------------------------

  profile_long <- do.call(
    rbind,
    lapply(
      colnames(normalized_scores),
      function(cl) {

        data.frame(
          Program = rownames(normalized_scores),
          Cluster = cl,
          Intensity = as.numeric(
            normalized_scores[, cl]
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  # --------------------------------------------------
  # Preserve biological program order
  # --------------------------------------------------

  program_order <- c(
    "Monocyte-associated",
    "Resident-associated",
    "Inflammatory",
    "Repair-associated",
    "Lipid-associated",
    "ECM-remodeling",
    "IFN-responsive"
  )

  profile_long$Program <- factor(
    profile_long$Program,
    levels = program_order
  )
  # --------------------------------------------------
  # Plot
  # --------------------------------------------------

  p <- ggplot2::ggplot(
    profile_long,
    ggplot2::aes(
      x = Program,
      y = Intensity,
      group = Cluster,
      color = Cluster
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = NULL,
      y = "Relative program intensity",
      color = "Cluster",
      title = "Multidimensional macrophage profiles",
      subtitle =
        "Program-intensity profiles normalized to unit Euclidean norm for visualization"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
      legend.position = "right"
    )

  return(p)
}

#' Plot macrophage profile similarity
#'
#' Visualizes pairwise similarity between macrophage profiles as a heatmap.
#' Similarity values are obtained from the profile similarity matrix generated
#' by \code{\link{analyze_macrophage_profile}} and represent cosine similarity
#' between multidimensional program-score profiles.
#'
#' Higher values indicate greater similarity between profiles. The diagonal
#' represents the similarity of each profile with itself.
#'
#' @param result Output list returned by
#'   \code{\link{analyze_macrophage_profile}}.
#'
#' @return A \code{ggplot} object containing a heatmap of pairwise macrophage
#'   profile similarity.
#'
#' @examples
#' \dontrun{
#' result <- analyze_macrophage_profile(
#'   expression = expression_matrix,
#'   clusters = cluster_labels
#' )
#'
#' plot_macrophage_similarity(result)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_viridis_c labs coord_equal theme_classic
#' @export
plot_macrophage_similarity <- function(result) {

  if (!is.list(result)) {
    stop("result must be the output of analyze_macrophage_profile().")
  }

  if (is.null(result$similarity)) {
    stop("result does not contain a similarity matrix.")
  }

  similarity <- result$similarity

  # --------------------------------------------------
  # Convert to long format
  # --------------------------------------------------

  similarity_long <- expand.grid(
    Cluster1 = rownames(similarity),
    Cluster2 = colnames(similarity),
    stringsAsFactors = FALSE
  )

  similarity_long$Similarity <-
    as.numeric(similarity)

  # --------------------------------------------------
  # Preserve cluster order on both axes
  # --------------------------------------------------

  cluster_order <- colnames(similarity)

  similarity_long$Cluster1 <- factor(
    similarity_long$Cluster1,
    levels = cluster_order
  )

  similarity_long$Cluster2 <- factor(
    similarity_long$Cluster2,
    levels = rev(cluster_order)
  )

  # --------------------------------------------------
  # Plot
  # --------------------------------------------------

  p <- ggplot2::ggplot(
    similarity_long,
    ggplot2::aes(
      x = Cluster1,
      y = Cluster2,
      fill = Similarity
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.2f", Similarity)
      ),
      size = 3
    ) +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "Cosine similarity",
      title = "Biological profile similarity"
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_classic()

  return(p)
}

#' Plot macrophage Framework annotations on a dimensional reduction
#'
#' Visualizes Macrophage Framework profile annotations on a dimensional
#' reduction stored in a Seurat object. The function produces two plots:
#' the original computational cluster assignments and the corresponding
#' multidimensional Macrophage Framework profile annotations.
#'
#' Cells are matched between the dimensional reduction and the supplied
#' cluster assignments using their cell barcodes. Framework profile labels
#' are obtained from \code{\link{summarize_macrophage_profile}}.
#'
#' Although the function is named \code{plot_macrophage_umap}, any Seurat
#' dimensional reduction containing at least two dimensions can be used.
#'
#' @param result Output list returned by
#'   \code{\link{analyze_macrophage_profile}}.
#' @param object Seurat object containing the dimensional reduction to plot.
#' @param clusters A vector of cell-level cluster assignments corresponding
#'   to \code{cell_names}.
#' @param cell_names Character vector of cell barcodes corresponding to the
#'   elements of \code{clusters}.
#' @param reduction Name of the Seurat dimensional reduction to use.
#'   Defaults to \code{"umap"}.
#'
#' @return A \code{patchwork} object containing two \code{ggplot} panels:
#'   the original computational clusters and the corresponding Macrophage
#'   Framework profile annotations.
#'
#' @examples
#' \dontrun{
#' result <- analyze_macrophage_profile(
#'   object = seurat_object,
#'   species = "Human"
#' )
#'
#' plot_macrophage_umap(
#'   result = result,
#'   object = seurat_object,
#'   clusters = seurat_object$seurat_clusters,
#'   cell_names = colnames(seurat_object),
#'   reduction = "umap"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point labs theme_classic theme
#' @export
plot_macrophage_umap <- function(result,
                                 object,
                                 clusters,
                                 cell_names,
                                 reduction = "umap") {

  # --------------------------------------------------
  # Input checks
  # --------------------------------------------------

  if (!is.list(result)) {
    stop("result must be the output of analyze_macrophage_profile().")
  }

  if (is.null(result$profile) ||
      is.null(result$profile$scores)) {
    stop("result does not contain a macrophage profile.")
  }

  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }

  if (!reduction %in% SeuratObject::Reductions(object)) {
    stop(
      paste0(
        "Reduction '", reduction,
        "' not found in the Seurat object."
      )
    )
  }

  if (length(clusters) != length(cell_names)) {
    stop(
      "clusters and cell_names must have the same length."
    )
  }

  # --------------------------------------------------
  # Extract dimensional reduction
  # --------------------------------------------------

  embeddings <- SeuratObject::Embeddings(
    object,
    reduction = reduction
  )

  if (ncol(embeddings) < 2) {
    stop(
      "The selected reduction must contain at least two dimensions."
    )
  }

  umap_df <- data.frame(
    Cell = rownames(embeddings),
    Dim1 = embeddings[, 1],
    Dim2 = embeddings[, 2],
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------
  # Create cluster table using supplied cell names
  # --------------------------------------------------

  cluster_df <- data.frame(
    Cell = cell_names,
    Cluster = as.character(clusters),
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------
  # Align clusters to dimensional reduction
  # --------------------------------------------------

  umap_df <- merge(
    umap_df,
    cluster_df,
    by = "Cell"
  )

  # --------------------------------------------------
  # Create framework profile labels
  # --------------------------------------------------

  summary_profile <- summarize_macrophage_profile(
    result = result,
    clusters = clusters
  )

  label_map <- setNames(
    summary_profile$ProfileLabel,
    summary_profile$Cluster
  )

  umap_df$ProfileLabel <- unname(
    label_map[as.character(umap_df$Cluster)]
  )

  # --------------------------------------------------
  # Original cluster plot
  # --------------------------------------------------

  p_clusters <- ggplot2::ggplot(
    umap_df,
    ggplot2::aes(
      x = Dim1,
      y = Dim2,
      color = Cluster
    )
  ) +
    ggplot2::geom_point(
      size = 0.35,
      alpha = 0.7
    ) +
    ggplot2::labs(
      x = paste0(toupper(reduction), "_1"),
      y = paste0(toupper(reduction), "_2"),
      color = "Cluster",
      title = "Original computational clusters"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right"
    )

  # --------------------------------------------------
  # Macrophage framework profile plot
  # --------------------------------------------------

  p_profiles <- ggplot2::ggplot(
    umap_df,
    ggplot2::aes(
      x = Dim1,
      y = Dim2,
      color = ProfileLabel
    )
  ) +
    ggplot2::geom_point(
      size = 0.35,
      alpha = 0.7
    ) +
    ggplot2::labs(
      x = paste0(toupper(reduction), "_1"),
      y = paste0(toupper(reduction), "_2"),
      color = "Macrophage profile",
      title = "Macrophage Framework profiles"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right"
    )

  # --------------------------------------------------
  # Combine plots
  # --------------------------------------------------

  combined_plot <- p_clusters + p_profiles

  return(combined_plot)
}

# ----------------------------------------------------------
# Generic heatmap helper
# ----------------------------------------------------------

.plot_heatmap <- function(
    matrix,
    filename,
    title,
    output_dir,
    zscore_rows = FALSE,
    limits = NULL,
    midpoint = NULL
) {

  if (is.null(matrix)) {
    return(character(0))
  }

  mat <- as.matrix(matrix)

  if (zscore_rows) {
    mat <- t(
      apply(
        mat,
        1,
        function(x) {
          s <- stats::sd(x, na.rm = TRUE)
          if (is.na(s) || s == 0) {
            rep(0, length(x))
          } else {
            as.numeric(
              (x - mean(x, na.rm = TRUE)) / s
            )
          }
        }
      )
    )
    rownames(mat) <- rownames(matrix)
    colnames(mat) <- colnames(matrix)
  }

  df <- expand.grid(
    Program = rownames(mat),
    Cluster = colnames(mat),
    stringsAsFactors = FALSE
  )

  df$Value <- as.vector(mat)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = Cluster,
      y = Program,
      fill = Value
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.2f", Value)
      ),
      size = 5
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 16
      ),
      axis.text.y = ggplot2::element_text(
        size = 17
      ),
      plot.title = ggplot2::element_text(
        size = 20,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 13
      ),
      legend.title = ggplot2::element_text(
        size = 14
      )
    )
  ggplot2::labs(
    title = title,
    x = NULL,
    y = NULL,
    fill = NULL
  )

  if (is.null(midpoint)) {
    midpoint <- if (!is.null(limits)) {
      mean(limits)
    } else {
      0
    }
  }

  if (!is.null(limits)) {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = midpoint,
        limits = limits
      )
  } else {
    p <- p +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = midpoint
      )
  }

  out_png <- file.path(
    output_dir,
    paste0(filename, ".png")
  )

  out_pdf <- file.path(
    output_dir,
    paste0(filename, ".pdf")
  )

  ggplot2::ggsave(
    out_png,
    p,
    width = 8,
    height = 6,
    dpi = 300
  )

  ggplot2::ggsave(
    out_pdf,
    p,
    width = 8,
    height = 6
  )

  c(out_png, out_pdf)
}


# ============================================================
# Internal: save tables and figures
# ============================================================

.save_framework_outputs <- function(
    result,
    embedding = NULL,
    clusters = NULL,
    output_dir = "MacrophageFramework_output"
) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  saved_files <- character(0)

  # ----------------------------------------------------------
  # Tables
  # ----------------------------------------------------------

  if (!is.null(result$macrophage_detection)) {

    detector_table <- result$macrophage_detection$cluster_metrics

  } else {

    detector_table <- data.frame(
      Cluster = sort(unique(as.character(clusters))),
      Status = ifelse(
        sort(unique(as.character(clusters))) %in%
          as.character(result$analyzed_clusters),
        "included",
        "excluded"
      ),
      stringsAsFactors = FALSE
    )

  }

  f <- file.path(output_dir, "01_macrophage_detector_clusters.csv")
  write.csv(detector_table, f, row.names = FALSE)
  saved_files <- c(saved_files, f)


  f <- file.path(output_dir, "02_framework_annotations.csv")
  write.csv(
    result$annotations,
    f,
    row.names = FALSE
  )
  saved_files <- c(saved_files, f)

  if (!is.null(result$scores)) {
    f <- file.path(output_dir, "03_program_intensity.csv")
    write.csv(
      as.data.frame(result$scores),
      f,
      row.names = TRUE
    )
    saved_files <- c(saved_files, f)
  }

  if (!is.null(result$detection)) {
    f <- file.path(output_dir, "04_program_coverage_cell_level.csv")
    write.csv(
      as.data.frame(result$detection),
      f,
      row.names = TRUE
    )
    saved_files <- c(saved_files, f)
  }

  if (!is.null(result$coherence)) {
    f <- file.path(output_dir, "05_program_coherence_cell_level.csv")
    write.csv(
      as.data.frame(result$coherence),
      f,
      row.names = TRUE
    )
    saved_files <- c(saved_files, f)
  }

  if (!is.null(result$profile)) {

    f <- file.path(output_dir, "06_cluster_profile_scores.csv")
    write.csv(
      as.data.frame(result$profile$scores),
      f,
      row.names = TRUE
    )
    saved_files <- c(saved_files, f)

    f <- file.path(output_dir, "07_cluster_profile_coverage.csv")
    write.csv(
      as.data.frame(result$profile$coverage),
      f,
      row.names = TRUE
    )
    saved_files <- c(saved_files, f)

    f <- file.path(output_dir, "08_cluster_profile_coherence.csv")
    write.csv(
      as.data.frame(result$profile$coherence),
      f,
      row.names = TRUE
    )
    saved_files <- c(saved_files, f)
  }

  if (!is.null(result$similarity)) {
    f <- file.path(output_dir, "09_profile_similarity.csv")
    write.csv(
      as.data.frame(result$similarity),
      f,
      row.names = TRUE
    )
    saved_files <- c(saved_files, f)
  }

  saveRDS(
    result,
    file = file.path(output_dir, "MacrophageFramework_result.rds")
  )
  saved_files <- c(
    saved_files,
    file.path(output_dir, "MacrophageFramework_result.rds")
  )

  # ----------------------------------------------------------
  # Figures
  # ----------------------------------------------------------

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning(
      "ggplot2 is not installed. Tables and RDS were saved, ",
      "but figures were not generated."
    )
    return(saved_files)
  }

  # ----------------------------------------------------------
  # UMAP 1: original computational clusters
  # ----------------------------------------------------------

  if (!is.null(embedding) && !is.null(clusters)) {

    umap_df <- data.frame(
      Cell = rownames(embedding),
      UMAP_1 = embedding[, 1],
      UMAP_2 = embedding[, 2],
      Cluster = as.character(clusters),
      stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(
      umap_df,
      ggplot2::aes(
        x = UMAP_1,
        y = UMAP_2,
        color = Cluster
      )
    ) +
      ggplot2::geom_point(
        size = 0.25,
        alpha = 0.6
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = "Original computational clusters",
        color = "Cluster"
      )

    f <- file.path(output_dir, "10_UMAP_original_clusters.png")
    ggplot2::ggsave(
      f,
      p,
      width = 8,
      height = 6,
      dpi = 300
    )
    saved_files <- c(saved_files, f)

    f <- file.path(output_dir, "10_UMAP_original_clusters.pdf")
    ggplot2::ggsave(
      f,
      p,
      width = 8,
      height = 6
    )
    saved_files <- c(saved_files, f)

    # --------------------------------------------------------
    # UMAP 2: detector status
    # --------------------------------------------------------

    status_lookup <- detector_table[, c("Cluster", "Status")]
    status_lookup$Cluster <- as.character(status_lookup$Cluster)

    umap_df$Status <- status_lookup$Status[
      match(umap_df$Cluster, status_lookup$Cluster)
    ]

    umap_df$Status[is.na(umap_df$Status)] <- "excluded"

    p <- ggplot2::ggplot(
      umap_df,
      ggplot2::aes(
        x = UMAP_1,
        y = UMAP_2,
        color = Status
      )
    ) +
      ggplot2::geom_point(
        size = 0.25,
        alpha = 0.6
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = "Macrophage detector",
        color = "Detector status"
      )

    f <- file.path(output_dir, "11_UMAP_macrophage_detector.png")
    ggplot2::ggsave(
      f,
      p,
      width = 8,
      height = 6,
      dpi = 300
    )
    saved_files <- c(saved_files, f)

    f <- file.path(output_dir, "11_UMAP_macrophage_detector.pdf")
    ggplot2::ggsave(
      f,
      p,
      width = 8,
      height = 6
    )
    saved_files <- c(saved_files, f)

    # --------------------------------------------------------
    # UMAP 3: Framework annotation
    # --------------------------------------------------------

    annotation_lookup <- result$annotations[
      ,
      c("Cluster", "Full_annotation"),
      drop = FALSE
    ]

    annotation_lookup$Cluster <-
      as.character(annotation_lookup$Cluster)

    umap_df$Framework_annotation <-
      annotation_lookup$Full_annotation[
        match(
          umap_df$Cluster,
          annotation_lookup$Cluster
        )
      ]

    umap_df$Framework_annotation[
      is.na(umap_df$Framework_annotation)
    ] <- "Excluded"

    p <- ggplot2::ggplot(
      umap_df,
      ggplot2::aes(
        x = UMAP_1,
        y = UMAP_2,
        color = Framework_annotation
      )
    ) +
      ggplot2::geom_point(
        size = 0.25,
        alpha = 0.6
      ) +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = "MacrophageFramework annotation",
        color = "Framework annotation"
      )

    f <- file.path(
      output_dir,
      "12_UMAP_framework_annotation.png"
    )
    ggplot2::ggsave(
      f,
      p,
      width = 10,
      height = 6,
      dpi = 300
    )
    saved_files <- c(saved_files, f)

    f <- file.path(
      output_dir,
      "12_UMAP_framework_annotation.pdf"
    )
    ggplot2::ggsave(
      f,
      p,
      width = 10,
      height = 6
    )
    saved_files <- c(saved_files, f)
  }



  if (!is.null(result$profile)) {

    saved_files <- c(
      saved_files,
      .plot_heatmap(
        result$profile$scores,
        "13_heatmap_program_scores",
        "Program scores",
        output_dir = output_dir,
        zscore_rows = FALSE,
        limits = c(
          0,
          max(result$profile$scores, na.rm = TRUE)
        )
      )
    )

    saved_files <- c(
      saved_files,
      .plot_heatmap(
        result$profile$coverage,
        "14_heatmap_program_coverage",
        "Program coverage",
        output_dir = output_dir,
        zscore_rows = FALSE,
        limits = c(0, 1)
      )
    )

    saved_files <- c(
      saved_files,
      .plot_heatmap(
        result$profile$coherence,
        "15_heatmap_program_coherence",
        "Program coherence",
        output_dir = output_dir,
        zscore_rows = FALSE,
        limits = c(0, 1)
      )
    )

    saved_files <- c(
      saved_files,
      .plot_heatmap(
        result$similarity,
        "16_heatmap_profile_similarity",
        "Profile similarity",
        output_dir = output_dir,
        zscore_rows = FALSE,
        limits = c(0, 1)
      )
    )
  }

  return(saved_files)
}
