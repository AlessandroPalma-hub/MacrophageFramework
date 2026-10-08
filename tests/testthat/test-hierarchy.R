test_that("hierarchy ranking separates dimensions", {

  ontology <- load_ontology()

  hierarchy_scores <- matrix(
    c(
      # Cluster_1
      5, 8, 2, 7, 1,

      # Cluster_2
      5, 2, 8, 1, 7
    ),
    nrow = 5,
    ncol = 2,
    dimnames = list(
      c(
        "Core",
        "Monocyte-associated",
        "Resident-associated",
        "Inflammatory",
        "Lipid-associated"
      ),
      c(
        "Cluster_1",
        "Cluster_2"
      )
    )
  )

  hierarchy_profile <- list(
    scores = hierarchy_scores
  )

  result <- rank_programs_by_dimension(
    hierarchy_profile,
    ontology = ontology
  )

  # Check dimensions are present
  expect_true("Identity" %in% names(result))
  expect_true("Lineage" %in% names(result))
  expect_true("Functional State" %in% names(result))

  # Identity
  expect_equal(
    result$Identity$Cluster_1,
    "Core"
  )

  expect_equal(
    result$Identity$Cluster_2,
    "Core"
  )

  # Lineage ranking
  expect_equal(
    result$Lineage$Cluster_1,
    c(
      "Monocyte-associated",
      "Resident-associated"
    )
  )

  expect_equal(
    result$Lineage$Cluster_2,
    c(
      "Resident-associated",
      "Monocyte-associated"
    )
  )

  # Functional-state ranking
  expect_equal(
    result$`Functional State`$Cluster_1,
    c(
      "Inflammatory",
      "Lipid-associated"
    )
  )

  expect_equal(
    result$`Functional State`$Cluster_2,
    c(
      "Lipid-associated",
      "Inflammatory"
    )
  )

})


test_that("hierarchy handles missing scores", {

  ontology <- load_ontology()

  scores <- matrix(
    c(
      5, NA,
      3, NA
    ),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c(
        "Core",
        "Monocyte-associated"
      ),
      c(
        "Cluster_1",
        "Cluster_2"
      )
    )
  )

  profile <- list(
    scores = scores
  )

  result <- rank_programs_by_dimension(
    profile,
    ontology = ontology
  )

  expect_equal(
    result$Identity$Cluster_1,
    "Core"
  )

  expect_length(
    result$Lineage$Cluster_2,
    0
  )

})
