# Spatial Transcriptomics Analysis of Tertiary Lymphoid Structures

Systematic spatial profiling of tertiary lymphoid structures (TLS) across multiple cancer types using Visium spatial transcriptomics reveals conserved organizational principles and immune cell coordination patterns.

## Overview

This repository contains analysis scripts and example data for characterizing TLS spatial organization through:

- **Cell state deconvolution** using EcoTyper framework
- **Spatial enrichment modeling** via mixed-effect logistic regression with Random Forest pre-selection
- **Chemokine network analysis** across three key axes (CXCL13-CXCR5, CXCL12-CXCR4, CCL19/21-CCR7)
- **Functional characterization** of cell states with divergent TLS associations
- **Single-cell validation** of spatial findings

## Repository Structure
```
TLS_spatial_analysis/
├── README.md
├── REQUIREMENTS.md
│
├── data/
│   └── example_data/
│       ├── DATA_STRUCTURE.md
│       ├── all_marker_gene.rds
│       ├── combined_obj.RData         (DOWNLOAD SEPARATELY, 4.3 GB)
│       ├── multi_regression.RData
│       ├── LR_analysis.RData
│       └── EcoTyper_Output/
│           └── Ecotypes/
│               └── ecotype_abundance.txt
│
├── scripts/
│   ├── 00_EcoTyper_preparation.R
│   ├── 01_TLSvsNonTLS_EcoTyper.R
│   ├── 02_Logistic.R
│   ├── 02_CellState_DIstanceGradient.R
│   ├── 02_CellState_Functional_Characterization.R
│   ├── 03_LR_Network_Analysis.R
│   ├── 03_Abundance_vs.Centrality.R
│   ├── 03_sensitivity_check_directionality.R
│   └── 04_PMN_analysis.R
│
└── output/
    ├── 00_gene_avg_by_sample.txt
    ├── 00_sample_histology.txt
    ├── 01_CE9_CE10_MannWhitney_boxplot.pdf
    ├── 01_CE_abundance_from_TLS.pdf
    ├── 01_Ecotype_TLS_comparison_with_pvalues.pdf
    ├── 02_volcano_main_effects.pdf
    ├── 02_volcano_plot.pdf
    ├── 02_forest_interactions.pdf
    ├── 02b_cellstate_distance_gradient.pdf
    ├── 02b_cellstate_composition.pdf
    ├── 03_fig3d_CD4T_S01S02_Core_Edge_correlation.pdf
    ├── 03_CXCL13-CXCR5_Network_Prominence.pdf
    ├── 03_CXCL13-CXCR5_Hub_community.pdf
    ├── 03_CXCL12_CXCR4_Network_Prominence.pdf
    ├── 03_CXCL12_CXCR4_Hub_community.pdf
    ├── 03_CCL19_21_CCR7_Network_Prominence.pdf
    ├── 03_CCL19_21_CCR7_Hub_community.pdf
    ├── 03_C2_TLS_vs_NonTLS.pdf
    ├── 03_Hub_Lineage_Diversity_C1.pdf
    ├── 03_Hub_Lineage_Diversity_C2.pdf
    ├── 03b_abundance_vs_centrality.pdf
    ├── 03b_abundance_vs_centrality_per_axis.pdf
    ├── 03c_hub_directionality_combined.pdf
    ├── 03c_hub_network_role.pdf
    ├── 04_PMN_validation_combined.pdf
    └── SFig.3c_hub_community_gradient.pdf
```

## Quick Start

### 1. Clone Repository
```bash
git clone https://github.com/anka2029/TLS_spatial_analysis.git
cd TLS_spatial_analysis
```

### 2. Install Dependencies

See [REQUIREMENTS.md](REQUIREMENTS.md) for complete installation guide.

**Essential packages:**
```r
# Core data manipulation
install.packages(c("dplyr", "tidyr", "purrr", "tibble", "stringr"))

# Statistical modeling
install.packages(c("ranger", "hstats", "lme4", "broom.mixed", "car"))

# Network analysis
install.packages(c("igraph", "circlize", "ggraph"))

# Visualization
install.packages(c("ggplot2", "ggrepel", "ggpubr", "pheatmap",
                    "patchwork", "viridis", "scales"))

# Single-cell analysis (for Script 04 only)
install.packages("Seurat")

# Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap"))
```

### 3. Download Data

Due to file size limitations, `combined_obj.RData` (4.3 GB) is not included in this repository. It contains all processed Visium samples (breast, kidney, cervical, prostate cancer).

Download from Zenodo: https://doi.org/10.5281/zenodo.18490685

### 4. Run Analysis
```r
setwd("path/to/TLS_spatial_analysis/scripts")

# Run in order
source("00_EcoTyper_preparation.R")       # Optional: prepare EcoTyper input
source("01_TLSvsNonTLS_EcoTyper.R")       # Ecotype-TLS association
source("02_Logistic.R")                    # Cell state-TLS association
source("02_CellState_DIstanceGradient.R")  # Spatial gradient analysis
source("02_CellState_Functional_Characterization.R")  # Functional analysis
source("03_LR_Network_Analysis.R")         # Chemokine network analysis
source("03_Abundance_vs.Centrality.R")     # Abundance vs centrality
source("03_sensitivity_check_directionality.R")  # Sensitivity & directionality
source("04_PMN_analysis.R")               # Optional: requires scRNA-seq data
```

## Analysis Pipeline

### Script 00: Prepare Input for EcoTyper Bulk Mode (Optional)

**Purpose**: Generate pseudo-bulk gene expression profiles by averaging Visium spot-level expression within TLS and non-TLS regions for each sample.

**Input**: `combined_obj.RData`

**Methods**: Average gene expression calculated separately for TLS and non-TLS areas. TLS spots are labeled by a pathologist based on H&E images.

**Key outputs**: `00_gene_avg_by_sample.txt`, `00_sample_histology.txt` (direct input for EcoTyper bulk mode)

**Runtime**: ~2 minutes

---

### Script 01: TLS vs Non-TLS Ecotype Comparison

**Purpose**: Identify carcinoma ecotypes (CEs) enriched in TLS regions, with focus on CE9 and CE10.

**Input**: `EcoTyper_Output/Ecotypes/ecotype_abundance.txt`, `combined_obj.RData`

**Methods**: Wilcoxon rank-sum test for TLS vs non-TLS comparison at both bulk and spot levels, distance gradient analysis using array coordinates.

**Key outputs**:
- `01_Ecotype_TLS_comparison_with_pvalues.pdf` — CE abundance comparison across all ecotypes
- `01_CE9_CE10_MannWhitney_boxplot.pdf` — CE9/CE10 spot-level association
- `01_CE_abundance_from_TLS.pdf` — Spatial CE abundance gradient from TLS

**Runtime**: ~2-5 minutes

---

### Script 02: Cell State Association with TLS (Logistic Regression)

**Purpose**: Identify cell states and their interactions associated with TLS presence using a two-stage approach: Random Forest variable selection followed by mixed-effect logistic regression.

**Input**: `combined_obj.RData`

**Methods**:
- Random Forest (ranger) with permutation importance for variable pre-selection (top 30)
- H-statistic (hstats) for pairwise interaction detection
- Multivariable mixed-effect logistic regression (glmer) with sample-level random intercept
- VIF check for multicollinearity

**Key outputs**:
- `02_volcano_main_effects.pdf` — Volcano plot of cell state associations (FDR < 0.05)
- `02_forest_interactions.pdf` — Forest plot of interaction effects
- `multi_regression.RData` — Saved model results for downstream analyses

**Runtime**: ~10-20 minutes

---

### Script 02 (continued): Cell State Distance Gradient

**Purpose**: Analyze spatial distribution of TLS-associated cell states relative to TLS centroids, cell state composition by cell type, and hub community spatial gradients.

**Input**: `combined_obj.RData`, `multi_regression.RData`

**Methods**:
- Connected component analysis to identify discrete TLS regions and compute geometric centroids
- Distance binning: TLS Core/Edge (median split) and non-TLS quartile categories
- Gradient statistics: median abundance per sample per bin, then mean +/- SE across samples

**Key outputs**:
- `02b_cellstate_distance_gradient.pdf` — TLS-enriched and depleted cell state gradients
- `02b_cellstate_composition.pdf` — Cell state composition by cell type with TLS association
- `SFig.3c_hub_community_gradient.pdf` — C1 vs C2 hub community normalized abundance gradient

**Runtime**: ~5-10 minutes

---

### Script 02 (continued): Functional Characterization of Cell States

**Purpose**: Characterize functional differences among cell states with divergent TLS associations within the same cell type (Macrophage and CD4.T cell states).

**Input**: `combined_obj.RData`, `multi_regression.RData`, `all_marker_gene.rds`

**Methods**:
- DEG heatmaps of marker gene expression per state
- GO Biological Process enrichment via compareCluster
- Spot-level Spearman correlation for spatial co-occurrence of TLS-enriched states
- CD4.T S01 vs S02 functional program correlation in TLS Core vs Edge zones

**Key outputs**:
- DEG heatmaps per cell type (Macrophage, CD4.T)
- GO enrichment dot plots per cell type
- Co-occurrence correlation heatmap of TLS-enriched states
- `03_fig3d_CD4T_S01S02_Core_Edge_correlation.pdf` — Functional program heatmap

**Runtime**: ~5-10 minutes

---

### Script 03: Ligand-Receptor Network Analysis

**Purpose**: Characterize chemokine-mediated cell state coordination across three axes: CXCL13-CXCR5, CXCL12-CXCR4, and CCL19/21-CCR7.

**Input**: `combined_obj.RData`, `LR_analysis.RData` (pre-computed networks), `all_marker_gene.rds`

**Methods**:
- Spatially-weighted interaction score calculation with exponential decay (lambda = 2 array units)
- Cross-sample hubness analysis (median weighted degree and edge frequency)
- Hub community identification via Jaccard-corrected Louvain clustering
- GO enrichment for hub communities
- Directional interaction heatmaps (sender-receiver relationships)
- Hub cell state interaction partner lineage diversity analysis

**Key outputs**:
- `03_*_Network_Prominence.pdf` (3 files) — Network prominence bubble plots
- `03_*_Hub_community.pdf` (3 files) — Hub community chord diagrams
- `03_C2_TLS_vs_NonTLS.pdf` — C2 hub members TLS vs non-TLS comparison
- `03_Hub_Lineage_Diversity_C1.pdf`, `03_Hub_Lineage_Diversity_C2.pdf`

**Runtime**: With pre-computed networks ~10-15 minutes; from scratch ~1-2 hours

---

### Script 03 (continued): Abundance vs Centrality

**Purpose**: Test whether chemokine network hub status is driven by cell state abundance in TLS or reflects genuine interaction specificity.

**Input**: `combined_obj.RData`, `LR_analysis.RData`

**Methods**: Scatter plot of mean TLS abundance (z-score) vs median weighted degree (z-score) across all cell states, with Spearman correlation. PMNs_S03 (low abundance, high centrality) is highlighted as a key outlier.

**Key outputs**:
- `03b_abundance_vs_centrality.pdf` — Combined scatter across axes
- `03b_abundance_vs_centrality_per_axis.pdf` — Per-axis supplementary scatter

**Runtime**: ~2-5 minutes

---

### Script 03 (continued): Sensitivity Analysis and Directionality

**Purpose**: Assess robustness of hub membership across edge weight thresholds (85th/90th/95th percentile) and decompose hub cell state interactions into sender vs receiver components.

**Input**: `combined_obj.RData`, `LR_analysis.RData`

**Methods**:
- Sensitivity: re-run network construction at three thresholds, compare hub membership via Jaccard similarity
- Hub threshold justification: degree distribution drop-off analysis
- Directionality: weighted out-degree vs in-degree decomposition per hub per axis

**Key outputs**:
- `03c_hub_directionality_combined.pdf` — Stacked bar + dot plot of sender/receiver roles
- `03c_hub_network_role.pdf` — Out-degree proportion summary

**Note**: Requires `analyze_LR_network()` and `get_hub_modules()` from `03_LR_Network_Analysis.R`. Run that script first or source it.

**Runtime**: Sensitivity analysis ~1-2 hours (re-computes networks); directionality ~5 minutes

---

### Script 04: PMN Validation (Optional)

**Purpose**: Validate PMNs_S03 marker gene expression in single-cell RNA-seq data from TLS-positive kidney cancer samples.

**Input**: User-provided scRNA-seq Seurat object, `all_marker_gene.rds`

**Methods**:
- Marker gene projection and module score calculation
- Average expression analysis across cell types
- Identification of common and neutrophil-specific marker genes

**Key outputs**: `04_PMN_validation_combined.pdf` — UMAP with cell types and module score

**Note**: Requires user-provided scRNA-seq data. The script contains placeholder paths that must be updated. Cell type annotation column name ("Subset") may need adjustment.

**Runtime**: ~10-20 minutes (depends on scRNA-seq data size)

## Data Format

See [data/example_data/DATA_STRUCTURE.md](data/example_data/DATA_STRUCTURE.md) for detailed specifications of all data files.

## Citation

If you use this code or data, please cite:
```bibtex

```

## Related Resources

- **EcoTyper**: [GitHub](https://github.com/digitalcytometry/ecotyper)
- **Visium**: [10x Genomics](https://www.10xgenomics.com/products/spatial-gene-expression)

## Contact

- **Author**: Ange Yan
- **Institution**: University of Tokyo
- **Email**: 7949284134@edu.k.u-tokyo.ac.jp
- **Issues**: [GitHub Issues](https://github.com/anka2029/TLS_spatial_analysis/issues)

## License


## Acknowledgments


---

**Last updated**: March 2026
**Status**: Active development | Manuscript in revision