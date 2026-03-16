# Example Data Structure Documentation

This document describes the structure and contents of example data files used in the TLS spatial transcriptomics analysis pipeline.

---

## Overview

The `data/example_data/` directory contains processed spatial transcriptomics data, EcoTyper cell state classifications, and pre-computed analysis results from a multi-cancer TLS study. These files enable reproducible execution of analysis scripts (00-04) without requiring access to raw sequencing data or time-consuming computation steps.

**Data Source**: Visium spatial transcriptomics from multiple cancer types (breast, kidney, cervical, prostate) with pathologist-annotated TLS regions.

**Processing Pipeline**: Raw Visium data → Space Ranger → EcoTyper (bulk & Visium modes) → Custom TLS annotation → Downstream analysis

---

## File Inventory

| File/Directory | Type | Approx. Size | Required For | Description |
|----------------|------|--------------|--------------|-------------|
| `all_marker_gene.rds` | RDS | ~56 KB | Scripts 02, 03, 04 | EcoTyper-defined cell state marker genes |
| `combined_obj.RData` | RData | ~4.3 GB | Scripts 00-04 | Multi-sample Visium data with annotations (download separately) |
| `EcoTyper_Output/` | Directory | ~10-50 MB | Script 01 | EcoTyper bulk mode results |
| `multi_regression.RData` | RData | ~5-20 MB | Scripts 02b, 02_Functional | Logistic regression model results |
| `LR_analysis.RData` | RData | ~5.7 MB | Scripts 03, 03b, 03c | Pre-computed L-R network analysis |

---

## Detailed File Structures

### 1. all_marker_gene.rds

**Purpose**: Marker genes defining each EcoTyper cell state for functional analysis

**Type**: R data frame (tibble)

**Structure**:
```r
# Data frame with 2 columns
Columns:
  - cell_state (character): Cell state identifier from EcoTyper
    Format: "[CellType]_S[##]"
    Examples: "CD8.T_S01", "CD4.T_S02", "PMNs_S03", "PCs_S01",
              "Monocytes.and.Macrophages_S01"

  - Gene (character): Gene symbol (HGNC nomenclature)
    Examples: "CD3D", "CD8A", "S100A8", "CXCL13"
```

**Dimensions**: ~2,000-10,000 rows (varies by number of cell states x markers per state)

**Example Usage**:
```r
gene_tbl <- readRDS("all_marker_gene.rds")
head(gene_tbl)
#   cell_state                         Gene
#   CD8.T_S01                          CD3D
#   CD8.T_S01                          CD8A
#   CD8.T_S01                          GZMA
#   PMNs_S03                           S100A8
#   PMNs_S03                           S100A9
#   Monocytes.and.Macrophages_S01      CD68
```

**Used In**:
- `02_CellState_Functional_Characterization.R`: DEG heatmaps, GO enrichment, co-occurrence analysis
- `03_LR_Network_Analysis.R`: Hub community functional characterization
- `04_PMN_analysis.R`: PMNs_S03 marker projection

**Notes**:
- Defined by EcoTyper
- Genes are deduplicated within each cell state

---

### 2. combined_obj.RData

**Purpose**: Integrated multi-sample Visium spatial transcriptomics data with TLS annotations and EcoTyper results (Visium mode)

**Type**: Named list of lists (nested structure)

**Top-level Structure**:
```r
combined_obj <- list(
  "B26" = list(...),      # Sample 1 (Breast cancer)
  "C01" = list(...),      # Sample 2 (Cervical cancer)
  ...                     # Additional samples
)
```

**Per-sample Structure**:
Each sample is a list containing 4 elements:

#### a) `$cellstate_raw` (Data frame)
Raw, non-normalized cell state abundance scores from EcoTyper Visium mode
```r
Columns:
  - ID (character): Spot barcode
    Example: "AAACAAGTATCTCCCA-1"

  - array_x (numeric): X-coordinate in array space
  - array_y (numeric): Y-coordinate in array space

  - Label (character): TLS annotation
    Values: "TLS" or "NonTLS"
    Source: Pathologist manual annotation

  - [CellState]_S## (numeric): Abundance for each cell state
    Examples:
      CD8.T_S01, CD8.T_S02
      CD4.T_S01, CD4.T_S02
      PCs_S01 (Plasma Cells)
      PMNs_S01, PMNs_S03 (Neutrophils)
      Monocytes.and.Macrophages_S01, Monocytes.and.Macrophages_S08
      B_S01
      NK_S01
      Fibroblasts_S01
      Endothelial_S01
      ... (71 cell states total)

Dimensions: ~2,000-5,000 rows (spots) x ~78 columns
```

#### b) `$cellstate_norm` (Data frame)
Normalized cell state abundances (same structure as cellstate_raw)
```r
# Normalization: Proportional abundance across all cell states within each spot (0-1)
# Used for: Most downstream analyses requiring comparable scales across spots
```

#### c) `$ecotype` (Data frame)
Carcinoma ecotype (CE) abundance per spot from EcoTyper Visium mode
```r
Columns:
  - ID (character): Spot barcode (matches cellstate tables)
  - CE1, CE2, CE3, ..., CE10 (numeric): Ecotype abundance

Dimensions: Same number of rows as cellstate_raw x 11 columns (ID + 10 CEs)
```

#### d) `$hs` (Sparse matrix - dgCMatrix)
Gene expression data for each spot
```r
Class: dgCMatrix (from Matrix package)
Dimensions: genes (rows) x spots (columns)
  Rows: 15,000-30,000 genes
  Cols: 2,000-5,000 spots (matches cellstate_raw rows)

Row names: Gene symbols (HGNC format)
  Examples: "CD3D", "EPCAM", "CXCL13", "KRT8"

Column names: Spot barcodes
  Example: "AAACAAGTATCTCCCA-1"

Values: Log-normalized expression
  Normalization: log1p(CPM) or log1p(counts/size_factor)
  Range: Typically 0-10 for most genes

Storage: Sparse format (only non-zero values stored)
```

**Example Access**:
```r
load("combined_obj.RData")

# Access sample B26
sample_b26 <- combined_obj$B26

# View cell state data
head(sample_b26$cellstate_norm[, 1:8])
#   ID                    array_x  array_y  Label   CD8.T_S01  CD4.T_S02  PCs_S01  PMNs_S03
#   AAACAAGTATCTCCCA-1    12       34       TLS     0.25       0.18       0.42     0.08
#   AAACACCAATAACTGC-1    14       36       NonTLS  0.05       0.03       0.01     0.12

# Check ecotype distribution
summary(sample_b26$ecotype[, c("CE9", "CE10")])

# Gene expression for specific genes
expr <- sample_b26$hs[c("CXCL13", "CD3D", "CD8A"), ]
```

**Used In**: All scripts (00-04)

**Quality Metrics**:
- Spots per sample: 2,000-5,000 (depends on tissue size)
- TLS spots: Typically 5-30% of total spots
- Genes detected: ~15,000-25,000 protein-coding genes
- Cell states: ~71 states

---

### 3. EcoTyper_Output/

**Purpose**: EcoTyper bulk mode deconvolution results for TLS vs NonTLS comparison

**Directory Structure**:
```
EcoTyper_Output/
├── B.cells/
│   ├── state_abundances.txt
│   ├── state_assignment_heatmap.pdf
│   ├── state_assignment_heatmap.png
│   └── state_assignment.txt
├── CD4.T.cells/
│   └── ...
├── CD8.T.cells/
│   └── ...
├── Dendritic.cells/
│   └── ...
├── Ecotypes/
│   ├── ecotype_abundance.txt
│   ├── ecotype_assignment.txt
│   ├── heatmap_assigned_samples_viridis.pdf
│   └── heatmap_assigned_samples_viridis.png
...
```

**File: ecotype_abundance.txt**

**Format**: Tab-separated text file

**Structure**:
```
Sample_ID        CE1      CE2      CE3      ...  CE10
Sample1_TLS      0.123    0.045    0.234    ...  0.456
Sample1_NonTLS   0.087    0.156    0.098    ...  0.023
Sample2_TLS      0.145    0.034    0.267    ...  0.489
Sample2_NonTLS   0.091    0.178    0.112    ...  0.019
...

Columns:
  - Sample_ID: Sample identifier with TLS/NonTLS suffix
    Format: "[SampleName]_TLS" or "[SampleName]_NonTLS"

  - CE1-CE10: Ecotype abundance
    Interpretation: Fraction of tissue composition for each ecotype
```

**Dimensions**:
- Rows: 2 x N_samples (TLS + NonTLS per sample)
- Columns: 11 (Sample_ID + 10 ecotypes)

**Used In**: `01_TLSvsNonTLS_EcoTyper.R` (Part 1)

---

### 4. multi_regression.RData

**Purpose**: Pre-computed results from the Random Forest + mixed-effect logistic regression analysis of cell state associations with TLS presence

**Type**: Multiple R objects saved together

**Contents**:
```r
# Objects saved in this file:
multi_results   # Tidy results from glmer (data frame)
rf_fit          # Fitted ranger Random Forest object
vi              # Variable importance vector (named numeric)
H               # hstats object (interaction statistics)
h2_pair         # Pairwise H-statistics
h2_df           # H-statistic data frame (interaction pairs)
```

**multi_results structure** (primary object used downstream):
```r
# Data frame from broom.mixed::tidy() with added columns
Columns:
  - term (character): Cell state name or interaction term
    Examples: "CD4.T_S02", "B_S01", "CD4.T_S02:Dendritic_S01"

  - estimate (numeric): Exponentiated coefficient (Odds Ratio)
  - std.error (numeric): Standard error
  - statistic (numeric): Wald z-statistic
  - p.value (numeric): p-value from Wald test
  - conf.low (numeric): Lower 95% CI for OR
  - conf.high (numeric): Upper 95% CI for OR
  - FDR (numeric): FDR-adjusted p-value
  - is_interaction (logical): Whether term is an interaction

Dimensions: ~40 rows (30 main effects + interactions) x 9 columns
```

**Example Usage**:
```r
load("multi_regression.RData")

# View significant main effects
multi_results %>%
  filter(!is_interaction, FDR < 0.05) %>%
  arrange(p.value)

# View significant interactions
multi_results %>%
  filter(is_interaction, p.value < 0.05)

# Check variable importance from Random Forest
head(sort(vi, decreasing = TRUE), 10)
```

**Used In**:
- `02_CellState_DIstanceGradient.R`: Select top enriched/depleted states for gradient analysis
- `02_CellState_Functional_Characterization.R`: Classify TLS association direction per cell state

---

### 5. LR_analysis.RData

**Purpose**: Pre-computed ligand-receptor interaction networks for three chemokine axes

**Type**: Nested list structure

**Top-level Structure**:
```r
LR_analysis <- list(
  CXCL13_CXCR5 = list(
    "Sample1" = list(...),
    "Sample2" = list(...),
    ...
  ),
  CXCL12_CXCR4 = list(
    "Sample1" = list(...),
    "Sample2" = list(...),
    ...
  ),
  CCL19_21_CCR7 = list(
    "Sample1" = list(...),
    "Sample2" = list(...),
    ...
  )
)
```

**Per-axis, per-sample Structure**:
Each element contains:

#### a) `$graph` (igraph object)
Directed interaction network
```r
Class: igraph (from igraph package)

Vertices:
  - Names: Cell state names (e.g., "CD8.T_S01", "PCs_S01")
  - Attributes:
    $name: Cell state identifier
    $community: Community ID from Louvain clustering
    $degree: Node degree
    $betweenness: Betweenness centrality
    $pagerank: PageRank score

Edges:
  - Directed: Sender -> Receiver
  - Attributes:
    $score: Interaction strength
      Calculation: Spatially-weighted product of ligand expression
                   in sender and receptor expression in receiver,
                   with exponential decay (lambda = 2 array units)
```

#### b) `$score` (Matrix)
Interaction score matrix (sender x receiver)
```r
Class: Numeric matrix
Dimensions: N_cell_states x N_cell_states

Row names: Sender cell states
Column names: Receiver cell states

Values: Interaction scores (>= 0)
  0: No interaction (or below threshold)
  >0: Interaction strength

Interpretation:
  score[i, j] = strength of i -> j interaction
```

#### c) `$node_metrics` (Data frame)
Node-level network statistics
```r
Columns:
  - cell_state (character): Cell state name
  - degree (integer): Total number of connections (in + out)
  - betweenness (numeric): Betweenness centrality
  - pagerank (numeric): PageRank centrality
  - community (factor): Community assignment from Louvain clustering

Dimensions: N_cell_states x 5 columns
```

#### d) `$edge_list` (Data frame)
Filtered edge list
```r
Columns:
  - sender (character): Sender cell state
  - receiver (character): Receiver cell state
  - score (numeric): Interaction score (above 90th percentile threshold)
```

**Example Usage**:
```r
load("LR_analysis.RData")

# Access CXCL13-CXCR5 network for sample B26
cxcl13_b26 <- LR_analysis$CXCL13_CXCR5$B26

# View graph summary
summary(cxcl13_b26$graph)

# Check interaction scores
cxcl13_b26$score["CD8.T_S01", "PCs_S01"]

# Node metrics
head(cxcl13_b26$node_metrics)
#   cell_state  degree  betweenness  pagerank  community
#   B_S01       72      380.56       0.083     1
#   PCs_S01     47      11.02        0.043     1
```

**Used In**:
- `03_LR_Network_Analysis.R`: Hubness analysis, community detection, GO enrichment, heatmaps
- `03_Abundance_vs.Centrality.R`: Weighted degree computation
- `03_sensitivity_check_directionality.R`: Directionality decomposition

**Computation Time** (if recomputing from scratch):
- Per axis: ~20-40 minutes
- All three axes: ~1-2 hours
- Recommendation: Use pre-computed for initial analysis

---

## Data Provenance

### Original Data Collection
- **Technology**: 10x Genomics Visium Spatial Gene Expression
- **Tissue types**: FFPE sections from surgical resections
- **TLS annotation**: Manual review by pathologists using H&E morphology

### Processing Steps
1. **Sequencing**: Illumina NovaSeq (paired-end reads)
2. **Alignment**: Space Ranger (10x Genomics) -> gene x spot count matrix
3. **Quality control**: Filter low-quality spots, doublets, high mitochondrial content
4. **Normalization**: Log-normalization
5. **TLS annotation**: Pathologist review -> Label column (TLS/NonTLS)
6. **EcoTyper analysis**:
   - Bulk mode: Deconvolve pseudo-bulk TLS/NonTLS profiles
   - Visium mode: Assign cell state abundances per spot
7. **Post-processing**: Combine into `combined_obj` structure

### Data Availability
- **Processed data**: Available from Zenodo (https://doi.org/10.5281/zenodo.18490685)
- **Raw data**: Not included

---

## Coordinate Systems

### Array Coordinates (used in all analyses)
- **array_x, array_y**: Array-space coordinates from Visium grid
  - Standardized spacing: 100 um between spot centers
  - Used for all distance calculations (consistent across samples)
  - Distance formula: `distance_um = sqrt((ax1 - ax2)^2 + (ay1 - ay2)^2) * 100`

### TLS Labels
- **Source**: Pathologist manual annotation using H&E morphology
- **Criteria for TLS**:
  - Dense lymphoid aggregates
  - Presence of B cell follicle-like structures
  - High endothelial venules (HEVs) when visible
  - Germinal center formation (mature TLS)
- **NonTLS**: All other tissue regions

---

## Troubleshooting

### Common Issues

**Problem**: "Sample X has no TLS spots"
```r
# Check TLS availability
tls_counts <- sapply(combined_obj, function(s)
  sum(s$cellstate_raw$Label == "TLS"))
valid_samples <- names(tls_counts[tls_counts > 0])
```

**Problem**: "Matrix operation failed - sparse matrix issue"
```r
# Convert to dense if needed
expr_dense <- as.matrix(combined_obj$B26$hs)

# Or work with sparse matrix directly
library(Matrix)
expr_sparse <- combined_obj$B26$hs  # Keep as dgCMatrix
```

**Problem**: "Cell state names don't match"
```r
# Check naming conventions
colnames(combined_obj$B26$cellstate_norm)
# Some functions use "." instead of " " in names
# Use gsub() to standardize if needed
```

---

## Version Information

**Data version**: 2.0
**Last updated**: 2026-03
**Compatible with**: R >= 4.0, Scripts 00-04

**Package requirements**:
- Matrix (sparse matrices)
- igraph (network analysis)
- lme4, ranger (statistical modeling)

---
*End of DATA_STRUCTURE documentation*