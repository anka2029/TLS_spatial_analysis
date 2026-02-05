# Spatial Transcriptomics Analysis of Tertiary Lymphoid Structures

Systematic spatial profiling of tertiary lymphoid structures (TLS) across multiple cancer types using Visium spatial transcriptomics reveals conserved organizational principles and immune cell coordination patterns.

## 📋 Overview

This repository contains analysis scripts and example data for characterizing TLS spatial organization through:

- **Cell state deconvolution** using EcoTyper framework
- **Spatial enrichment modeling** via Geographically Weighted Regression (GWR)
- **Chemokine network analysis** across three key axes (CXCL13-CXCR5, CXCL12-CXCR4, CCL19/21-CCR7)
- **Single-cell validation** of spatial findings

## 🗂️ Repository Structure
```
TLS_spatial_analysis/
├── README.md                          # Project overview (you are here)
├── REQUIREMENTS.md                    # Package dependencies and installation guide
│
├── data/
│   └── example_data/
│       ├── DATA_STRUCTURE.md          # Detailed data format specifications
│       ├── DOWNLOAD.md                # Instructions for full dataset download
│       ├── all_marker_gene.rds        # EcoTyper cell state markers (56 KB)
│       ├── combined_obj.RData         # Full dataset (DOWNLOAD SEPARATELY, 4.3 GB)
│       ├── gwr_results.RData          # Pre-computed GWR models (91 MB)
│       ├── LR_analysis.RData          # Pre-computed L-R networks (5.7 MB)
│       └── EcoTyper_Output/           # Ecotype deconvolution results
│           └── Ecotypes/
│               └── ecotype_abundance.txt
│           └── B.cells
│               └── state_abundances.txt
│               └── state_assignment_heatmap.pdf
│               └── state_assignment_heatmap.png
│               └── state_assignment.txt
│               ...
├── scripts/
│   ├── 00_EcoTyper_preparation.R    # prepare pseudo-bulk input for EcoTyper 
│   ├── 01_TLSvsNonTLS_EcoTyper.R    # TLS vs non-TLS ecotype comparison
│   ├── 02_GWRmodel_cellstate.R       # Spatial enrichment analysis with GWR
│   ├── 03_LR_Network_Analysis.R      # Chemokine network characterization
│   └── 04_PMN_validation.R           # scRNA-seq validation of PMNs_S03
│
└── output/                            # Generated figures (created by scripts)
    ├── 00_gene_avg_by_sample.txt
    ├── 00_sample_histology.txt
    ├── 01_CE9_CE10_MannWhitney_boxplot.pdf
    ├── 01_CE_abundance_from_TLS.pdf
    ├── 01_Ecotype_TLS_comparison_with_pvalues.pdf
    ├── 02_Impact_of_cell_states.pdf
    ├── 02_macrophage_spatial_associations.pdf
    ├── 02_spatial_distribution.pdf
    ├── 02_Tcell_Plasmacell_analysis.pdf
    ├── 03_CCL19_21_CCR7_Hub_community.pdf
    ├── 03_CCL19_21_CCR7_Network_Prominence.pdf
    ├── 03_CXCL12_CXCR4_Hub_community.pdf
    ├── 03_CXCL12_CXCR4_Network_Prominence.pdf
    ├── 03_CXCL13-CXCR5_Hub_community.pdf
    ├── 03_CXCL13-CXCR5_Network_Prominence.pdf
    └── 04_PMN_validation_combined.pdf
```

## 🚀 Quick Start

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
install.packages(c("dplyr", "tidyr", "purrr", "tibble"))

# Spatial analysis
install.packages(c("Matrix", "sp", "GWmodel"))

# Network analysis
install.packages(c("igraph", "circlize"))

# Single-cell analysis
install.packages("Seurat")

# Visualization
install.packages(c("ggplot2", "ggpubr", "pheatmap", "patchwork"))

# Bioconductor packages
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap"))
```

### 3. Download Data

**Full dataset (DOWNLOAD REQUIRED)**
- All cancer samples (breast, kidney, cervical, prostate, melanoma)
- Complete analysis and validation
- File: `combined_obj.RData` (4.3 GB)
- **Download instructions**: See [data/example_data/DOWNLOAD.md](data/example_data/DOWNLOAD.md)

### 4. Run Analysis
```r
# Set working directory
setwd("path/to/TLS_spatial_analysis")

# Execute scripts in order
source("scripts/01_TLSvsNonTLS_EcoTypern.R")
source("scripts/02_GWRmodel_cellstate.R")
source("scripts/03_LR_Network_Analysis.R")
source("scripts/04_PMN_validation.R")  # Optional: requires scRNA-seq data
```

## 📊 Analysis Pipeline

### Script 00: Prepare input for EcoTyper bulk mode (optional)

**Purpose**: 
This script generates pseudo-bulk gene expression profiles by averaging Visium spot-level expression within TLS and Non-TLS regions for each sample.

**Input**:
- `combined_obj.RData`

**Methods**:
- Average gene expression for TLS/NonTLS area separately
- TLS spots are labeled by a pathologist

**Key outputs**:
- `gene_avg_by_sample.txt`
    - Gene expression matrix (genes x samples)
    - Columns: sample_TLS / sample_NonTLS
- `sample_histology.txt`
    - Metadata file mapping each column to:
    - TLS or NonTLS
These two files are used as direct input for EcoTyper bulk mode.

**Runtime**: ~2 minutes

---

### Script 01: TLS vs Non-TLS Ecotype Comparison

**Purpose**: Identify cell states and ecotypes enriched in TLS regions

**Input**:
- `EcoTyper_Output/Ecotypes/ecotype_abundance.txt`

**Methods**:
- Wilcoxon rank-sum test for TLS vs non-TLS comparison
- Multiple testing correction (FDR)
- Effect size calculation

**Key outputs**:
- `01_Ecotype_TLS_comparison_with_pvalues.pdf`
    - Figures comparing CE abundance in TLS vs NonTLS 
- `01_CE9_CE10_MannWhitney_boxplot.pdf`
    - CE9 / CE10 TLS association statistics
- `01_CE_abundance_from_TLS.pdf`
    - Spatial CE abundance gradient plots relative to TLS distance 

**Runtime**: ~2-5 minutes

---

### Script 02: Geographically Weighted Regression Analysis

**Purpose**: Model spatial enrichment patterns of cell states related to TLS presence

**Input**:
- `combined_obj.RData`
- `all_marker_gene.rds`

**Methods**:
- Two steps variable screening(cell state) via univariate/multivariate logistic regression
- GWR with adaptive bandwidth selection
- Spatial autocorrelation analysis
- Gene Ontology enrichment for TLS-enriched states

**Key outputs**:
- `02_Impact_of_cell_states.pdf`
    - Barplot of cell state impact on TLS (median |β| and fraction significant spots)   
- `02_macrophage_spatial_associations.pdf`
    - Macrophages S01 vs S08 local regression coefficients across samples shows consistent positive association of S01 with TLS
- `02_Tcell_Plasmacell_analysis.pdf`
    - Violin plots of local coefficient distributions (CD8+ T S01, CD4+ T S02, PC S01)
    - Statistical significance vs effect size scatter plots
- `02_spatial_distribution.pdf`
    - Spatial visualization of CD4+ T S02, CD8+ T S01, PC S01, and cancer cells
    - B26 sample with TLS regions outlined in red
- `02_go_plot.pdf`
    - GO enrichment analysis (Biological Process) for high-abundance spot DEGs
    - Functional characterization of CD8+ T S01, CD4+ T S02, and PC S01
  
**Runtime**: 
- With pre-computed results: ~5-10 minutes
- From scratch: ~30-120 minutes

---

### Script 03: Ligand-Receptor Network Analysis

**Purpose**: Characterize chemokine-mediated cell state coordination

**Input**:
- `combined_obj.RData` 
- `LR_analysis.RData` (pre-computed networks)

**Methods**:
- Three chemokine axes analysis:
  - CXCL13-CXCR5 
  - CXCL12-CXCR4 
  - CCL19/21-CCR7
- Network centrality metrics (degree, betweenness, closeness)
- Community detection for hub identification

**Key outputs**:
- `03_*_Network_Prominence.pdf`, 3 files
  - Network prominence plots
- `03_*_Hub_community.pdf`, 3 files 
  - Hub community chord diagrams
- Chord diagrams for L-R interactions
- GO enrichment for network communities

**Runtime**:
- With pre-computed networks: ~10-15 minutes
- From scratch: ~1-2 hours

---

### Script 04: PMN Validation (Optional)

**Purpose**: Validate PMNs_S03 marker genes in scRNA-seq data

**Input**:
- Your own scRNA-seq data from TLS+ samples
- `all_marker_gene.rds`

**Methods**:
- Marker gene projection onto UMAP
- Cell-type-specific expression analysis

**Key outputs**:
- UMAP with marker expression
- Violin plots per cell type

**Runtime**: ~10-20 minutes (depends on scRNA-seq data size)

**Note**: Requires user-provided scRNA-seq data. Template included for adaptation.

## 🔍 Data Format

See [data/example_data/DATA_STRUCTURE.md](data/example_data/DATA_STRUCTURE.md) for detailed specifications.

**combined_obj** structure (Seurat-based):
```r
combined_obj <- list(
  "Sample1" = list(
    cellstate_raw  = data.frame(ID, pixel_x, pixel_y, Label, cell_states...),
    cellstate_norm = data.frame(normalized abundances),
    ecotype        = data.frame(CE1-CE10 abundances),
    hs             = dgCMatrix(genes × spots expression)
  ),
  "Sample2" = ...
)
```

## 📝 Citation

If you use this code or data, please cite:
```bibtex

```

## 📚 Related Resources

- **EcoTyper**: [GitHub](https://github.com/digitalcytometry/ecotyper)
- **GWmodel**: [CRAN](https://cran.r-project.org/package=GWmodel)



## 📧 Contact

- **Author**: Ange Yan
- **Institution**: University of Tokyo
- **Email**: 7949284134@edu.k.u-tokyo.ac.jp
- **Issues**: [GitHub Issues](https://github.com/anka2029/TLS_spatial_analysis/issues)

## 📄 License


## 🙏 Acknowledgments


---

**Last updated**: February 2026  
**Status**: Active development | Manuscript in revision
