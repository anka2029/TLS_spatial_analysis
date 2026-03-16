# Package Requirements

## R Version

- R >= 4.0.0

## Required Packages

### Core Data Manipulation
- dplyr (>= 1.0.0)
- tidyr (>= 1.0.0)
- purrr (>= 0.3.0)
- tibble (>= 3.0.0)
- stringr (>= 1.4.0)
- forcats

### Matrix Operations
- Matrix (>= 1.3.0) - Sparse matrices

### Statistical Modeling
- ranger (>= 0.14.0) - Random Forest
- hstats (>= 0.3.0) - H-statistic for interaction detection
- lme4 (>= 1.1.0) - Mixed-effect models (glmer)
- broom.mixed - Tidy model summaries for mixed models
- car - VIF and regression diagnostics
- rstatix - Pipe-friendly statistical tests

### Network Analysis
- igraph (>= 1.2.0) - Network graphs and community detection
- circlize (>= 0.4.0) - Chord diagrams
- ggraph - Network visualization

### Single-cell & Spatial Transcriptomics
- Seurat (>= 4.0.0) - Required for Script 04 only

### Visualization
- ggplot2 (>= 3.3.0)
- ggpubr
- ggsignif
- ggrepel
- patchwork
- viridis
- pheatmap
- ComplexHeatmap (Bioconductor)
- RColorBrewer
- scales

### Functional Enrichment
- clusterProfiler (>= 4.0.0) (Bioconductor)
- org.Hs.eg.db (Bioconductor)
- AnnotationDbi (Bioconductor)

### Utilities
- conflicted - Namespace conflict resolution
- readr

## Installation

```r
# Install CRAN packages
install.packages(c(
  # Core
  "dplyr", "tidyr", "purrr", "tibble", "stringr", "forcats",
  "Matrix", "conflicted", "readr",
  # Statistical modeling
  "ranger", "hstats", "lme4", "broom.mixed", "car", "rstatix",
  # Network
  "igraph", "circlize", "ggraph",
  # Visualization
  "ggplot2", "ggpubr", "ggsignif", "ggrepel", "patchwork",
  "viridis", "pheatmap", "RColorBrewer", "scales"
))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "ComplexHeatmap"
))

# Install Seurat (required for Script 04 only)
install.packages("Seurat")
```

## Package-Script Mapping

| Script | Key Packages |
|--------|-------------|
| 00_EcoTyper_preparation.R | Matrix, dplyr, purrr, tidyr |
| 01_TLSvsNonTLS_EcoTyper.R | dplyr, tidyr, ggplot2, ggpubr, ggsignif, rstatix |
| 02_Logistic.R | ranger, hstats, lme4, broom.mixed, car, ggrepel, patchwork |
| 02_CellState_DIstanceGradient.R | igraph, dplyr, tidyr, ggplot2, patchwork, stringr |
| 02_CellState_Functional_Characterization.R | ComplexHeatmap, clusterProfiler, org.Hs.eg.db, circlize, patchwork |
| 03_LR_Network_Analysis.R | igraph, circlize, ggraph, pheatmap, clusterProfiler, viridis, ggrepel |
| 03_Abundance_vs.Centrality.R | dplyr, tidyr, ggplot2, ggrepel, purrr |
| 03_sensitivity_check_directionality.R | igraph, dplyr, ggplot2, patchwork, stringr, scales |
| 04_PMN_analysis.R | Seurat, dplyr, tidyr, ggplot2, pheatmap, patchwork |

## Session Info Example

```r
sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 18.04.6 LTS
```