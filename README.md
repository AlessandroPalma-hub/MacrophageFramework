# MacrophageFramework

## Multidimensional macrophage profiling

MacrophageFramework is a multidimensional computational framework for profiling macrophage populations using an ontology of biological programs.

The framework is designed to capture continuous and combinatorial macrophage states rather than imposing binary activation categories. Macrophage populations are independently detected and subsequently characterized across multiple biological dimensions.

## Framework overview

MacrophageFramework evaluates macrophage profiles across seven biological programs:

- **Monocyte-associated**
- **Resident-associated**
- **Inflammatory**
- **Repair-associated**
- **Lipid-associated**
- **ECM-remodeling**
- **IFN-responsive**

The framework provides complementary quantitative measures of:

- program-level expression intensity;
- detection coverage;
- gene-program coherence;
- multidimensional macrophage profiles;
- profile similarity between populations.

Macrophage identity detection is intentionally separated from biological-state profiling. This separation allows the framework to characterize macrophage heterogeneity without using the same markers both to define macrophage identity and to assign functional states.

## Installation

The development version can currently be installed from GitHub:

```r
install.packages("remotes")
remotes::install_github("AlessandroPalma-hub/MacrophageFramework")
```

The package is under active development and is intended for submission to Bioconductor.

## Basic usage

For a Seurat object:

```r
library(MacrophageFramework)

result <- analyze_macrophage_profile(
  object = seu_object,
  species = "Human"
)
```

The resulting object contains macrophage detection results, biological program scores, detection coverage, coherence measures, multidimensional profiles, profile similarity, annotations, and optional dimensional-reduction information.

The main public functions are:

- `analyze_macrophage_profile()`
- `load_ontology()`
- `score_programs()`
- `calculate_detection_coverage()`
- `calculate_program_coherence()`
- `calculate_profile_similarity()`
- `summarize_macrophage_profile()`
- `plot_macrophage_profile()`
- `plot_macrophage_similarity()`
- `plot_macrophage_umap()`

## Ontology

The current ontology is version 0.4 and is distributed with the package in `inst/extdata/`.

The ontology distinguishes two biological dimensions:

- **Lineage**
- **Functional State**

Macrophage identity is handled independently by a dedicated macrophage detector and is not treated as a biological-state dimension of the framework.

## Development status

MacrophageFramework is currently under active development.

Version 0.99.0 represents the development version being prepared for Bioconductor submission.

## Citation

A formal citation will be added when the associated manuscript and software release are finalized.

## License

GPL-3
