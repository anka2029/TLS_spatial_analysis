# Package Requirements

## R Version
- R ≥ 4.0.0

## Required Packages

### Core Data Manipulation
- dplyr (≥ 1.0.0)
- tidyr (≥ 1.0.0)
- purrr (≥ 0.3.0)
- tibble (≥ 3.0.0)

### Spatial & Matrix Operations
- Matrix (≥ 1.3.0) - Sparse matrices
- sp (≥ 1.4.0) - Spatial data structures
- GWmodel (≥ 2.2.0) - Geographically weighted regression

### Network Analysis
- igraph (≥ 1.2.0) - Network graphs
- circlize (≥ 0.4.0) - Chord diagrams

### Single-cell & Spatial Transcriptomics
- Seurat (≥ 4.0.0)
- spacexr (optional for 04)

### Visualization
- ggplot2 (≥ 3.3.0)
- ggpubr
- ggsignif
- ggrepel
- patchwork
- viridis
- pheatmap
- ComplexHeatmap
- RColorBrewer

### Functional Enrichment
- clusterProfiler (≥ 4.0.0)
- org.Hs.eg.db
- enrichplot
- msigdbr

### Statistics
- broom
- rstatix
- limma

### Utilities
- scales
- stringr
- conflicted
- hdf5r
- biomaRt
- compositions

## Installation
```r
# Install CRAN packages
install.packages(c(
  "dplyr", "tidyr", "purrr", "tibble",
  "Matrix", "sp", "GWmodel",
  "igraph", "circlize",
  "ggplot2", "ggpubr", "ggsignif", "ggrepel", "patchwork",
  "viridis", "pheatmap", "RColorBrewer",
  "broom", "rstatix", "scales", "stringr",
  "conflicted", "hdf5r", "compositions"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "limma",
  "ComplexHeatmap",
  "biomaRt",
  "msigdbr"
))

# Install Seurat
install.packages("Seurat")

# Optional: spacexr (for script 04)
# devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
```

## Package-Script Mapping

| Script | Key Packages |
|--------|-------------|
| 01_TLSvsNonTLS_EcoTypern.R | dplyr, tidyr, ggplot2, ggpubr, ggsignif |
| 02_GWRmodel_cellstate.R | GWmodel, sp, limma, clusterProfiler, patchwork |
| 03_LR_Network_Analysis.R | igraph, circlize, pheatmap, clusterProfiler |
| 04_PMN_validation.R | Seurat, ComplexHeatmap, spacexr (optional) |

## Session Info Example
```r
sessionInfo()
#R version 4.4.2 (2024-10-31)
#Platform: x86_64-pc-linux-gnu
#Running under: Ubuntu 18.04.6 LTS
```