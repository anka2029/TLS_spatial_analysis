# ============================================================
# 04_PMN_validation.R
#
# Purpose:
#   Validate PMNs_S03 marker gene expression in single-cell RNA-seq data
#   from TLS-positive kidney cancer samples.
#
#   This script includes:
#   1) Projection of PMNs_S03 marker genes onto single-cell data
#   2) Average expression analysis across cell types
#   3) Module score calculation and UMAP visualization
#   4) Analysis of PMNs_S03 gene expression in Macrophage, Fibroblast, and Neutrophil
#   5) Identification of common and cell-type-specific marker genes
#
# Input:
#   (1) your_sc_data: Single-cell RNA-seq Seurat object
#       (Replace with your own scRNA-seq data)
#   
#   (2) Marker gene table:
#       ../data/example_data/all_marker_gene.rds
#       (Contains EcoTyper defined cell_state information and Gene columns)
#
# Output:
#   - Combined UMAP plots (04_PMN_validation_combined.pdf)
#   - DotPlots of common and neutrophil-specific genes
#
# ============================================================

## ── Required packages ──────────────────────────────────────────
library(spacexr)
library(Seurat)
library(ggplot2)
library(hdf5r)
library(tidyverse)
library(ggrepel)
library(compositions)
library(ComplexHeatmap)
library(biomaRt)
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(org.Hs.eg.db)   # For human genes
library(enrichplot)     # For dotplot, barplot, etc.

## ── Load data ──────────────────────────────────────────────────
# Load your single-cell data (replace with your own data)
your_sc_data <- readRDS("path/to/your_sc_data.RDS")

# Load marker genes
gene_tbl <- readRDS("../data/example_data/all_marker_gene.rds")


# ============================================================================
# PART 1: Project PMNs_S03 Marker Genes onto Single-Cell Data
# ============================================================================

## ── 1.1 Extract PMNs_S03 marker genes ────────────────────────────────────
PMNs03_genes <- gene_tbl %>% 
  filter(cell_state == "PMNs_S03") %>% 
  pull(Gene) %>% 
  unique()

## ── 1.2 Set cell type identity ───────────────────────────────────────────
# Set Subset as identity
Idents(your_sc_data) <- your_sc_data$Subset

## ── 1.3 Calculate average expression ─────────────────────────────────────
# Extract average expression (using RNA slot)
avg_exp <- AverageExpression(
  your_sc_data, 
  features = PMNs03_genes, 
  return.seurat = FALSE
)$RNA

# Create heatmap
pheatmap(
  avg_exp,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Average expression of PMNs S03 markers per subset"
)

## ── 1.4 Calculate and visualize summary scores ───────────────────────────
# Calculate summed expression score per cell type
subset_score <- colSums(avg_exp)

# Format data
score_df <- data.frame(
  Subset = names(subset_score),
  Score = as.numeric(subset_score)
) |>
  dplyr::arrange(desc(Score)) |>
  dplyr::mutate(Subset = factor(Subset, levels = Subset))

# Remove parentheses and their contents from Subset in dataframe
score_df <- score_df %>%
  mutate(Subset_clean = gsub("\\s*\\(.*?\\)", "", Subset))

# Create barplot
ggplot(score_df, aes(x = reorder(Subset_clean, Score), y = Score)) +
  geom_col(width = 0.7, fill = "#4B9CD3") +
  coord_flip() +  # Flip axes
  theme_minimal(base_size = 18) +
  labs(
    title = "PMNs_S03 Marker Expression \n (per cell type)",
    x = NULL, 
    y = "Summed Average Expression"
  ) +
  theme(
    axis.text.x = element_text(size = 16),  # Numeric side
    axis.text.y = element_text(size = 16, hjust = 1),  # Cell type name side
    axis.title.x = element_text(size = 18),  # "Summed Average Expression"
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    panel.grid.major.y = element_blank()  # Remove y-axis grid lines
  )


# ============================================================================
# PART 2: Add Module Score and Create UMAP Visualization
# ============================================================================

## ── 2.1 Calculate module score ───────────────────────────────────────────
your_sc_data <- AddModuleScore(
  your_sc_data, 
  features = list(PMNs03_genes),
  name = "PMNs03_marker_genes"
)

## ── 2.2 Create UMAP plots ────────────────────────────────────────────────
library(patchwork)

# p1: Cell Type UMAP (top)
p1 <- DimPlot(
  your_sc_data, 
  group.by = "Subset", 
  label = TRUE, 
  label.size = 4,
  repel = TRUE
) +
  ggtitle("Cell Type UMAP") + 
  NoLegend() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

# p2: PMNs_S03 Marker Expression (bottom)
p2 <- FeaturePlot(
  your_sc_data, 
  features = "PMNs03_marker_genes1", 
  reduction = "umap", 
  min.cutoff = 0
) + 
  ggtitle("PMNs_S03 Marker Genes") +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

# Save two plots stacked vertically
combined_plot <- p1 / p2

ggsave(
  "../output/04_PMN_validation_combined.pdf", 
  combined_plot, 
  width = 10, 
  height = 16  # Double height for vertical stacking
)


# ============================================================================
# PART 3: PMNs_S03 Gene Expression in Macrophage, Fibroblast & Neutrophil
# ============================================================================

library(dplyr)
library(tidyr)

## ── 3.1 Define cell types and extract expression ─────────────────────────
# Define three cell types of interest
cell_types_of_interest <- c("Monocyte / Macrophage", "Fibroblast", "Neutrophil")
Idents(your_sc_data) <- your_sc_data$Subset

# Get PMNs_S03 marker gene expression data
exp_data <- FetchData(
  your_sc_data, 
  vars = c(PMNs03_genes, "Subset"),
  layer = "data"
)

## ── 3.2 Calculate statistics per cell type ───────────────────────────────
# Calculate mean expression and percentage of expressing cells for each gene in each cell type
gene_stats <- exp_data %>%
  filter(Subset %in% cell_types_of_interest) %>%
  pivot_longer(cols = -Subset, names_to = "Gene", values_to = "Expression") %>%
  group_by(Subset, Gene) %>%
  summarise(
    mean_exp = mean(Expression, na.rm = TRUE),
    pct_exp = sum(Expression > 0) / n() * 100,
    .groups = "drop"
  )

## ── 3.3 Define thresholds and filter genes ───────────────────────────────
# Set thresholds (adjustable)
mean_threshold <- 0.1  # Threshold for mean expression
pct_threshold <- 1     # Threshold for percentage of expressing cells (%)

# Identify genes expressed in each cell type
expressed_genes <- gene_stats %>%
  filter(mean_exp > mean_threshold & pct_exp > pct_threshold) %>%
  select(Subset, Gene) %>%
  distinct()

# Genes expressed in all three cell types (common genes)
common_genes <- expressed_genes %>%
  group_by(Gene) %>%
  filter(n() == 3) %>%  # Expressed in all three
  pull(Gene) %>%
  unique()

# Neutrophil-specific genes
neutrophil_specific <- expressed_genes %>%
  group_by(Gene) %>%
  filter(n() == 1 & Subset == "Neutrophil") %>%  # Expressed only in Neutrophil
  pull(Gene) %>%
  unique()

## ── 3.4 Visualize gene expression patterns ───────────────────────────────
# Plot 1: DotPlot of commonly expressed genes
if(length(common_genes) > 0) {
  p1 <- DotPlot(
    your_sc_data, 
    features = common_genes,
    group.by = "Subset",
    idents = cell_types_of_interest
  ) + 
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    ggtitle("Common Genes Across 3 cell types")
  
  print(p1)
}

# Plot 2: DotPlot of neutrophil-specific genes
if(length(neutrophil_specific) > 0) {
  p2 <- DotPlot(
    your_sc_data, 
    features = neutrophil_specific,
    group.by = "Subset",
    idents = cell_types_of_interest
  ) + 
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    ggtitle("Neutrophil-Specific PMNs_S03 Marker Genes")
  
  print(p2)
}

# ============================================================================
# END OF SCRIPT
# ============================================================================