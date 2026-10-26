## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE
)

library(MacrophageFramework)

## ----load-data----------------------------------------------------------------
demo_file <- system.file(
  "extdata",
  "GSE114725_demo.rds",
  package = "MacrophageFramework"
)

if (!nzchar(demo_file) || !file.exists(demo_file)) {
  stop("The demonstration Seurat object was not found in the installed package.")
}

demo_seurat <- readRDS(demo_file)

dim(demo_seurat)
colnames(demo_seurat[[]])
Seurat::Reductions(demo_seurat)
SeuratObject::Layers(demo_seurat[["RNA"]])

## ----select-clusters----------------------------------------------------------
Seurat::Idents(demo_seurat) <- "seurat_clusters"

table(Seurat::Idents(demo_seurat))

## ----framework-analysis-------------------------------------------------------
demo_result <- analyze_macrophage_profile(
  object = demo_seurat,
  assay = "RNA",
  layer = "counts",
  reduction = "umap",
  detect_macrophages = TRUE,
  min_core_score = 0.1,
  species = "Human",
  save_outputs = FALSE
)

names(demo_result)
demo_result$analyzed_clusters
length(demo_result$analyzed_cells)

## ----detector-summary---------------------------------------------------------
detector_table <- demo_result$macrophage_detection$cluster_metrics

if (!is.null(detector_table)) {
  knitr::kable(
    detector_table[
      order(as.numeric(as.character(detector_table$Cluster))),
      c(
        "Cluster",
        "n_cells",
        "macrophage_core_score",
        "macrophage_core_coverage",
        "macrophage_core_coherence",
        "pct_detectable",
        "Status"
      )
    ],
    digits = 3,
    caption = "Cluster-level macrophage detector metrics for the demonstration dataset."
  )
}

## ----manual-clusters, eval=FALSE----------------------------------------------
# manual_result <- analyze_macrophage_profile(
#   object = demo_seurat,
#   assay = "RNA",
#   layer = "counts",
#   reduction = "umap",
#   detect_macrophages = FALSE,
#   macrophage_clusters = c("2", "3", "5"),
#   species = "Human",
#   save_outputs = FALSE
# )

## ----program-scores-----------------------------------------------------------
dim(demo_result$profile$scores)
round(demo_result$profile$scores, 3)

## ----program-coverage---------------------------------------------------------
dim(demo_result$profile$coverage)
round(demo_result$profile$coverage, 3)

## ----program-coherence--------------------------------------------------------
dim(demo_result$profile$coherence)
round(demo_result$profile$coherence, 3)

## ----coherence-variability----------------------------------------------------
round(demo_result$profile$coherence_IQR, 3)

## ----profile-similarity-------------------------------------------------------
dim(demo_result$similarity)
round(demo_result$similarity, 3)

## ----matrix-input, eval=FALSE-------------------------------------------------
# counts <- Seurat::GetAssayData(
#   demo_seurat,
#   assay = "RNA",
#   layer = "counts"
# )
# 
# cluster_labels <- Seurat::Idents(demo_seurat)[colnames(counts)]
# 
# matrix_result <- analyze_macrophage_profile(
#   expression = counts,
#   clusters = cluster_labels,
#   detect_macrophages = TRUE,
#   min_core_score = 0.1,
#   species = "Human",
#   save_outputs = FALSE
# )

## ----inspect-outputs, eval=FALSE----------------------------------------------
# str(demo_result, max.level = 1)

## ----session-info-------------------------------------------------------------
sessionInfo()

