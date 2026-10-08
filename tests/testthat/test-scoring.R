test_that("program scoring works", {

  ontology <- load_ontology()

  result <- score_programs(
    test_expression,
    ontology = ontology
  )

  expect_true(is.list(result))

  expect_true(
    all(c(
      "scores",
      "genes_available",
      "genes_total",
      "coverage"
    ) %in% names(result))
  )

  expect_true(is.matrix(result$scores))

  expect_equal(
    nrow(result$scores),
    length(unique(ontology$Program))
  )

})


test_that("program scores identify expected toy programs", {

  ontology <- load_ontology()

  result <- score_programs(
    test_expression,
    ontology = ontology
  )

  # Cell_A: Core
  expect_gt(
    result$scores["Core", "Cell_A"],
    0
  )

  # Cell_B: Lipid-associated
  expect_gt(
    result$scores["Lipid-associated", "Cell_B"],
    0
  )

  # Cell_C: Inflammatory
  expect_gt(
    result$scores["Inflammatory", "Cell_C"],
    0
  )

})


test_that("detection coverage works", {

  ontology <- load_ontology()

  detection <- calculate_detection_coverage(
    test_expression,
    ontology = ontology
  )

  expect_true(is.matrix(detection))

  expect_equal(
    detection["Core", "Cell_A"],
    1
  )

  expect_equal(
    detection["Lipid-associated", "Cell_B"],
    1
  )

  expect_equal(
    detection["Inflammatory", "Cell_C"],
    1
  )

})


test_that("cluster score summarization works", {

  ontology <- load_ontology()

  result <- score_programs(
    test_expression,
    ontology = ontology
  )

  cluster_scores <- summarize_cluster_scores(
    result$scores,
    cluster_labels
  )

  expect_true(is.matrix(cluster_scores))

  expect_equal(
    ncol(cluster_scores),
    length(unique(cluster_labels))
  )

})


test_that("coherence is bounded and handles dropout", {

  ontology <- load_ontology()

  coherence <- calculate_program_coherence(
    test_log,
    ontology = ontology
  )

  # All finite coherence values must lie between 0 and 1
  finite_values <- coherence[
    is.finite(coherence)
  ]

  expect_true(
    all(finite_values >= 0 & finite_values <= 1)
  )

})
