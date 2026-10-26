
test_that("hierarchy ranking separates dimensions", {
  
  ontology <- load_ontology("v0.4")
  
  programs <- c(
    "Monocyte-associated",
    "Resident-associated",
    "Inflammatory",
    "Repair-associated",
    "Lipid-associated",
    "ECM-remodeling",
    "IFN-responsive"
  )
  
  hierarchy_scores <- matrix(
    c(
      # Cluster_1
      8, 2, 7, 1, 6, 3, 4,
      
      # Cluster_2
      2, 8, 1, 7, 3, 6, 5
    ),
    nrow = length(programs),
    ncol = 2,
    dimnames = list(
      programs,
      c("Cluster_1", "Cluster_2")
    )
  )
  
  result <- rank_programs_by_dimension(
    list(scores = hierarchy_scores),
    ontology = ontology
  )
  
  expect_setequal(
    names(result),
    c("Lineage", "Functional State")
  )
  
  # Lineage ranking
  expect_equal(
    result$Lineage$Cluster_1,
    c("Monocyte-associated", "Resident-associated")
  )
  
  expect_equal(
    result$Lineage$Cluster_2,
    c("Resident-associated", "Monocyte-associated")
  )
  
  # Functional-state ranking
  expect_equal(
    result$`Functional State`$Cluster_1,
    c(
      "Inflammatory",
      "Lipid-associated",
      "IFN-responsive",
      "ECM-remodeling",
      "Repair-associated"
    )
  )
  
  expect_equal(
    result$`Functional State`$Cluster_2,
    c(
      "Repair-associated",
      "ECM-remodeling",
      "IFN-responsive",
      "Lipid-associated",
      "Inflammatory"
    )
  )
  
  expect_false("Identity" %in% names(result))
})


test_that("hierarchy handles missing scores", {
  
  ontology <- load_ontology("v0.4")
  
  scores <- matrix(
    c(
      5, NA,
      NA, NA
    ),
    nrow = 2,
    ncol = 2,
    dimnames = list(
      c(
        "Monocyte-associated",
        "Resident-associated"
      ),
      c("Cluster_1", "Cluster_2")
    )
  )
  
  result <- rank_programs_by_dimension(
    list(scores = scores),
    ontology = ontology
  )
  
  expect_equal(
    result$Lineage$Cluster_1,
    "Monocyte-associated"
  )
  
  expect_length(
    result$Lineage$Cluster_2,
    0
  )
  
  expect_false("Identity" %in% names(result))
})
