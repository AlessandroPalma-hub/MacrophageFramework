rank_programs_by_dimension <- function(profile,
                                       ontology = NULL) {

  # Check profile
  if (!is.list(profile)) {
    stop("profile must be a profile list.")
  }

  if (is.null(profile$scores)) {
    stop("profile must contain scores.")
  }

  scores <- profile$scores

  if (!is.matrix(scores)) {
    stop("profile$scores must be a matrix.")
  }

  # Load ontology if not provided
  if (is.null(ontology)) {
    ontology <- load_ontology()
  }

  # Check ontology
  required_columns <- c(
    "Dimension",
    "Program"
  )

  if (!all(required_columns %in% colnames(ontology))) {
    stop("Ontology is missing required columns.")
  }

  # Check that all score programs are in ontology
  missing_programs <- setdiff(
    rownames(scores),
    ontology$Program
  )

  if (length(missing_programs) > 0) {
    stop(
      paste(
        "The following programs are missing from ontology:",
        paste(missing_programs, collapse = ", ")
      )
    )
  }

  # Program -> dimension mapping
  program_dimension <- ontology[
    match(rownames(scores), ontology$Program),
    c("Program", "Dimension")
  ]

  # Store results
  result <- list()

  # Unique dimensions
  dimensions <- unique(program_dimension$Dimension)

  for (dimension in dimensions) {

    programs <- program_dimension$Program[
      program_dimension$Dimension == dimension
    ]

    dimension_scores <- scores[
      programs,
      ,
      drop = FALSE
    ]

    # Rank programs within each cluster
    ranking <- lapply(
      seq_len(ncol(dimension_scores)),
      function(i) {

        values <- dimension_scores[, i]

        observed <- !is.na(values)

        if (!any(observed)) {
          return(character(0))
        }

        # Explicitly preserve program names
        program_names <- rownames(dimension_scores)[observed]
        values <- values[observed]

        # Rank programs by score
        program_names[
          order(values, decreasing = TRUE)
        ]
      }
    )

    names(ranking) <- colnames(dimension_scores)

    result[[dimension]] <- ranking
  }

  return(result)
}
