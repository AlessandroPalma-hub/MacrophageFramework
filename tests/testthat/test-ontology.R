test_that("ontology loads correctly", {

  ontology <- load_ontology()

  expect_true(is.data.frame(ontology))

  expect_equal(
    nrow(ontology),
    29
  )

  expect_true(
    all(
      c(
        "Dimension",
        "Program",
        "Gene",
        "Weight",
        "Program_specificity",
        "Context_dependence"
      ) %in% colnames(ontology)
    )
  )

  expect_true(
    all(ontology$Weight > 0)
  )

  expect_true(
    all(!is.na(ontology$Program))
  )

})
