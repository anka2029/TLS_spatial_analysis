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
│       ├── combined_obj_small.RData   # Example dataset: 2-3 samples (~300 MB)
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
│   ├── 01_TLSvsNonTLS_EcoTyper.R    # TLS vs non-TLS ecotype comparison
│   ├── 02_GWRmodel_cellstate.R       # Spatial enrichment analysis with GWR
│   ├── 03_LR_Network_Analysis.R      # Chemokine network characterization
│   └── 04_PMN_validation.R           # scRNA-seq validation of PMNs_S03
│
└── output/                            # Generated figures (created by scripts)
    ├── 01_ecotype_comparison.pdf
    ├── 02_GWR_spatial_patterns.pdf
    ├── 03_network_analysis.pdf
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

### Script 01: TLS vs Non-TLS Ecotype Comparison

**Purpose**: Identify cell states and ecotypes enriched in TLS regions

**Input**:
- `EcoTyper_Output/Ecotypes/ecotype_abundance.txt`

**Methods**:
- Wilcoxon rank-sum test for TLS vs non-TLS comparison
- Multiple testing correction (FDR)
- Effect size calculation

**Key outputs**:
- Ecotype abundance heatmap
- Statistical comparison (p-values, fold-changes)
- Identifies CE9/CE10 as TLS-associated ecotypes

**Runtime**: ~2-5 minutes

---

### Script 02: Geographically Weighted Regression Analysis

**Purpose**: Model spatial enrichment patterns of cell states around TLS

**Input**:
- `combined_obj.RData` or `combined_obj_small.RData`
- `all_marker_gene.rds`

**Methods**:
- GWR with adaptive bandwidth selection
- Spatial autocorrelation analysis
- Gene Ontology enrichment for TLS-enriched states

**Key outputs**:
- GWR coefficient heatmaps
- Spatial distribution maps
- GO enrichment results for hub states

**Runtime**: 
- With pre-computed results: ~5-10 minutes
- From scratch: ~30-60 minutes

---

### Script 03: Ligand-Receptor Network Analysis

**Purpose**: Characterize chemokine-mediated cell state coordination

**Input**:
- `combined_obj.RData` or `combined_obj_small.RData`
- `LR_analysis.RData` (pre-computed networks)

**Methods**:
- Three chemokine axes analysis:
  - CXCL13-CXCR5 (B cell recruitment)
  - CXCL12-CXCR4 (plasma cell homing)
  - CCL19/21-CCR7 (T cell migration)
- Network centrality metrics (degree, betweenness, closeness)
- Community detection for hub identification

**Key outputs**:
- Chord diagrams for L-R interactions
- Network centrality plots
- Hub state characterization
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
- Spatial deconvolution validation (optional: requires RCTD/spacexr)

**Key outputs**:
- UMAP with marker expression
- Violin plots per cell type
- Validation of spatial findings

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
@article{YourName2025,
  title={Systematic spatial profiling reveals conserved organizational principles of tertiary lymphoid structures},
  author={Your Name and Collaborators},
  journal={Journal Name},
  year={2025},
  volume={XX},
  pages={XXX-XXX},
  doi={10.xxxx/xxxxx}
}
```

## 📚 Related Resources

- **EcoTyper**: [GitHub](https://github.com/digitalcytometry/ecotyper)
- **Visium Spatial**: [10x Genomics](https://www.10xgenomics.com/products/spatial-gene-expression)
- **GWmodel**: [CRAN](https://cran.r-project.org/package=GWmodel)
- **CellChat database**: [GitHub](https://github.com/sqjin/CellChat)

## 🤝 Contributing

Contributions are welcome! Please feel free to:
- Report bugs via [Issues](https://github.com/anka2029/TLS_spatial_analysis/issues)
- Suggest improvements
- Submit pull requests

## 📧 Contact

- **Author**: Anka
- **Institution**: University of Tokyo
- **Email**: 7949284134@edu.k.u-tokyo.ac.jp
- **Issues**: [GitHub Issues](https://github.com/anka2029/TLS_spatial_analysis/issues)

## 📄 License

This project is licensed under the MIT License - see below for details:
```
MIT License

Copyright (c) 2025 Anka

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## 🙏 Acknowledgments


---

**Last updated**: February 2026  
**Status**: Active development | Manuscript in revision
