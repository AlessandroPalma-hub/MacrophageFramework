# ============================================================
# Species-aware ontology gene resolution
# ============================================================

.resolve_ontology_genes <- function(
    ontology,
    expression_genes,
    species = "Human"
) {

  expression_genes <- as.character(expression_genes)

  # Explicit ortholog mappings that cannot be resolved
  # by case-insensitive matching alone.
  mouse_orthologs <- c(
    "FCGR1A" = "Fcgr1",
    "ADGRE1" = "Emr1"
  )

  resolved_gene <- vapply(
    ontology$Gene,
    function(gene) {

      # 1. Exact match
      if (gene %in% expression_genes) {
        return(gene)
      }

      # 2. Case-insensitive match
      idx <- which(
        toupper(expression_genes) == toupper(gene)
      )

      if (length(idx) == 1) {
        return(expression_genes[idx])
      }

      # 3. Explicit Mouse ortholog
      if (
        identical(species, "Mouse") &&
        gene %in% names(mouse_orthologs)
      ) {

        target <- mouse_orthologs[[gene]]

        idx <- which(
          toupper(expression_genes) == toupper(target)
        )

        if (length(idx) == 1) {
          return(expression_genes[idx])
        }
      }

      # 4. Unresolved
      NA_character_
    },
    character(1)
  )

  status <- ifelse(
    is.na(resolved_gene),
    "unresolved",
    ifelse(
      resolved_gene == ontology$Gene,
      "exact",
      "resolved"
    )
  )

  # Replace only genes that were successfully resolved.
  resolved_ontology <- ontology[
    !is.na(resolved_gene),
    ,
    drop = FALSE
  ]

  resolved_ontology$Gene <- resolved_gene[
    !is.na(resolved_gene)
  ]

  diagnostics <- data.frame(
    Ontology_gene = ontology$Gene,
    Resolved_gene = resolved_gene,
    Status = status,
    stringsAsFactors = FALSE
  )

  list(
    ontology = resolved_ontology,
    diagnostics = diagnostics
  )
}
