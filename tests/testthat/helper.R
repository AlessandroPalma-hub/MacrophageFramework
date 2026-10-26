
# Shared test fixtures for MacrophageFramework

test_ontology <- MacrophageFramework::load_ontology("v0.4")

test_genes <- unique(
  test_ontology$Gene[
    !is.na(test_ontology$Gene) &
      nzchar(test_ontology$Gene)
  ]
)

test_programs <- unique(test_ontology$Program)

test_cells <- c("Cell_A", "Cell_B", "Cell_C")

test_expression <- matrix(
  0,
  nrow = length(test_genes),
  ncol = length(test_cells),
  dimnames = list(test_genes, test_cells)
)

# Assign each program to a synthetic cell.
# Three programs are assigned explicitly for targeted checks.
program_cell <- setNames(
  rep(test_cells, length.out = length(test_programs)),
  test_programs
)

program_cell["Monocyte-associated"] <- "Cell_A"
program_cell["Lipid-associated"] <- "Cell_B"
program_cell["Inflammatory"] <- "Cell_C"

for (program in test_programs) {
  program_genes <- unique(
    test_ontology$Gene[
      test_ontology$Program == program
    ]
  )
  
  program_genes <- intersect(program_genes, test_genes)
  
  cell <- program_cell[[program]]
  
  if (length(program_genes) > 0) {
    test_expression[program_genes, cell] <- 10
  }
}

test_log <- log1p(test_expression)

cluster_labels <- c(
  "Cluster1",
  "Cluster1",
  "Cluster2"
)
