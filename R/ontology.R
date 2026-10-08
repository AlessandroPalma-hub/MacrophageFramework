#' Load the macrophage ontology
#'
#' Loads a predefined version of the Macrophage Framework ontology from the
#' package's installed \code{inst/extdata} directory.
#'
#' The current ontology version is \code{"v0.4"}. In this version, macrophage
#' identity is established independently by the macrophage detector and is
#' therefore not represented as an ontology dimension. The ontology defines
#' lineage-associated and functional-state programs used for Framework
#' scoring and annotation.
#'
#' @param version Character string specifying the ontology version to load.
#'   Defaults to \code{"v0.4"}.
#'
#' @return A data frame containing the selected Macrophage Framework ontology.
#'
#' @examples
#' ontology <- load_ontology()
#' head(ontology)
#'
#' @export
#'
load_ontology <- function(version = "v0.4") {

  ontology_path <- system.file(
    "extdata",
    paste0("ontology_", version, ".csv"),
    package = "MacrophageFramework"
  )

  if (ontology_path == "") {
    stop(
      "Ontology version '", version,
      "' was not found in the installed MacrophageFramework package."
    )
  }

  ontology <- utils::read.csv(
    ontology_path,
    stringsAsFactors = FALSE
  )

  return(ontology)
}
