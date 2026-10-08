#' Calculate similarity between macrophage profiles
#'
#' Calculates the cosine similarity between two Macrophage Framework profiles
#' based on their program scores. Only biological programs shared by both
#' profiles and with non-missing scores in both profiles are used.
#'
#' Cosine similarity measures the similarity in the overall pattern of
#' program scores independently of the absolute magnitude of the two
#' profiles. Values closer to 1 indicate more similar program profiles,
#' whereas values closer to 0 indicate less similar profiles.
#'
#' @param profile_1 A Macrophage Framework profile list containing a
#'   \code{scores} matrix.
#' @param profile_2 A second Macrophage Framework profile list containing a
#'   \code{scores} matrix.
#'
#' @return A numeric value representing the cosine similarity between the
#'   two profiles. \code{NA} is returned when fewer than two programs can be
#'   compared or when one of the profiles has zero magnitude.
#'
#' @examples
#' profile_1 <- list(
#'   scores = matrix(
#'     c(1, 0.5, 0.2),
#'     ncol = 1,
#'     dimnames = list(
#'       c("Program_A", "Program_B", "Program_C"),
#'       "Profile_1"
#'     )
#'   )
#' )
#'
#' profile_2 <- list(
#'   scores = matrix(
#'     c(0.8, 0.6, 0.1),
#'     ncol = 1,
#'     dimnames = list(
#'       c("Program_A", "Program_B", "Program_C"),
#'       "Profile_2"
#'     )
#'   )
#' )
#'
#' calculate_profile_similarity(profile_1, profile_2)
#'
#' @export
calculate_profile_similarity <- function(profile_1,
                                         profile_2) {

  # Check profiles
  if (!is.list(profile_1)) {
    stop("profile_1 must be a profile list.")
  }

  if (!is.list(profile_2)) {
    stop("profile_2 must be a profile list.")
  }

  # Check that scores exist
  if (is.null(profile_1$scores)) {
    stop("profile_1 must contain scores.")
  }

  if (is.null(profile_2$scores)) {
    stop("profile_2 must contain scores.")
  }

  scores_1 <- profile_1$scores
  scores_2 <- profile_2$scores

  # Check matrices
  if (!is.matrix(scores_1) || !is.matrix(scores_2)) {
    stop("Profile scores must be matrices.")
  }

  # Identify common programs
  common_programs <- intersect(
    rownames(scores_1),
    rownames(scores_2)
  )

  if (length(common_programs) < 2) {
    stop("At least two common programs are required.")
  }

  # Extract common programs
  vector_1 <- scores_1[common_programs, 1]
  vector_2 <- scores_2[common_programs, 1]

  # Keep programs with observed values in both profiles
  observed <- !is.na(vector_1) & !is.na(vector_2)

  vector_1 <- vector_1[observed]
  vector_2 <- vector_2[observed]

  # Need at least two comparable programs
  if (length(vector_1) < 2) {
    return(NA_real_)
  }

  # Cosine similarity
  denominator <-
    sqrt(sum(vector_1^2)) *
    sqrt(sum(vector_2^2))

  # Avoid division by zero
  if (denominator == 0) {
    return(NA_real_)
  }

  similarity <-
    sum(vector_1 * vector_2) / denominator

  return(similarity)
}

calculate_similarity_matrix <- function(profile) {

  # Check profile
  if (!is.list(profile)) {
    stop("profile must be a profile list.")
  }

  # Check scores
  if (is.null(profile$scores)) {
    stop("profile must contain scores.")
  }

  scores <- profile$scores

  if (!is.matrix(scores)) {
    stop("profile$scores must be a matrix.")
  }

  if (is.null(rownames(scores))) {
    stop("profile$scores must have program names as rownames.")
  }

  if (is.null(colnames(scores))) {
    stop("profile$scores must have cluster names as colnames.")
  }

  # Number of clusters
  n_clusters <- ncol(scores)

  # Initialize similarity matrix
  similarity <- matrix(
    NA_real_,
    nrow = n_clusters,
    ncol = n_clusters,
    dimnames = list(
      colnames(scores),
      colnames(scores)
    )
  )

  # Calculate pairwise similarity
  for (i in seq_len(n_clusters)) {

    for (j in seq_len(n_clusters)) {

      profile_i <- list(
        scores = scores[, i, drop = FALSE]
      )

      profile_j <- list(
        scores = scores[, j, drop = FALSE]
      )

      similarity[i, j] <- calculate_profile_similarity(
        profile_i,
        profile_j
      )
    }
  }

  return(similarity)
}
