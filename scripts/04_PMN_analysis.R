# ============================================================
# 04_PMN_analysis.R
#
# Purpose:
#   Validate PMNs_S03 marker gene expression in single-cell
#   RNA-seq data from TLS-positive kidney cancer samples.
#
#   This script includes:
#   1) Projection of PMNs_S03 marker genes onto single-cell data
#   2) Average expression analysis across cell types
#   3) Module score calculation and UMAP visualization
#   4) PMNs_S03 gene expression in Macrophage, Fibroblast,
#      and Neutrophil
#   5) Identification of common and cell-type-specific markers
#
# Input:
#   (1) Single-cell RNA-seq Seurat object
#       Replace the path below with your own scRNA-seq data.
#       The object must contain:
#         - RNA assay with normalized expression
#         - UMAP reduction
#         - A metadata column with cell type annotations
#           (default: "Subset"; adjust if your column differs)
#
#   (2) Marker gene table:
#       ../data/example_data/all_marker_gene.rds
#       Contains EcoTyper-defined cell_state and Gene columns
#
# Output:
#   04_PMN_marker_expression_barplot.pdf
#   04_PMN_validation_combined.pdf
#   04_PMN_common_genes_dotplot.pdf
#   04_PMN_neutrophil_specific_dotplot.pdf
#
# Dependencies:
#   Seurat, dplyr, tidyr, ggplot2, pheatmap, patchwork, conflicted
# ============================================================

library(conflicted)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(patchwork)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(base::intersect)

# ════════════════════════════════════════════════════════════
# DATA LOADING
# ════════════════════════════════════════════════════════════

# --- User-specific scRNA-seq data ---
# Replace this path with your own Seurat object.
# The object should have a metadata column for cell type
# annotations. The default column name used below is "Subset";
# change it throughout the script if yours differs.
your_sc_data <- readRDS("path/to/your_sc_data.RDS")

# --- Marker genes ---
gene_tbl <- readRDS("../data/example_data/all_marker_gene.rds")


# ════════════════════════════════════════════════════════════
# PART 1: PROJECT PMNs_S03 MARKER GENES ONTO SINGLE-CELL DATA
# ════════════════════════════════════════════════════════════

## ── 1.1 Extract PMNs_S03 marker genes ────────────────────
PMNs03_genes <- gene_tbl %>%
  filter(cell_state == "PMNs_S03") %>%
  pull(Gene) %>%
  unique()

cat("PMNs_S03 marker genes:", length(PMNs03_genes), "\n")

## ── 1.2 Set cell type identity ───────────────────────────
# NOTE: Replace "Subset" with your cell type annotation column
Idents(your_sc_data) <- your_sc_data$Subset

## ── 1.3 Calculate average expression ─────────────────────
avg_exp <- AverageExpression(
  your_sc_data,
  features = PMNs03_genes,
  return.seurat = FALSE
)$RNA

# Heatmap of average expression per cell type
pheatmap(
  avg_exp,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Average expression of PMNs_S03 markers per cell type"
)

## ── 1.4 Summarized expression score per cell type ────────
subset_score <- colSums(avg_exp)

score_df <- data.frame(
  Subset = names(subset_score),
  Score  = as.numeric(subset_score)
) %>%
  arrange(desc(Score)) %>%
  mutate(
    Subset       = factor(Subset, levels = Subset),
    Subset_clean = gsub("\\s*\\(.*?\\)", "", Subset)
  )

p_bar <- ggplot(score_df, aes(x = reorder(Subset_clean, Score),
                              y = Score)) +
  geom_col(width = 0.7, fill = "#4B9CD3") +
  coord_flip() +
  theme_minimal(base_size = 18) +
  labs(
    title = "PMNs_S03 Marker Expression\n(per cell type)",
    x     = NULL,
    y     = "Summed Average Expression"
  ) +
  theme(
    axis.text.x    = element_text(size = 16),
    axis.text.y    = element_text(size = 16, hjust = 1),
    axis.title.x   = element_text(size = 18),
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 20),
    panel.grid.major.y = element_blank()
  )

ggsave(p_bar,
       file   = "../output/04_PMN_marker_expression_barplot.pdf",
       width  = 10,
       height = 8,
       dpi    = 300)


# ════════════════════════════════════════════════════════════
# PART 2: MODULE SCORE AND UMAP VISUALIZATION
# ════════════════════════════════════════════════════════════

## ── 2.1 Calculate module score ───────────────────────────
your_sc_data <- AddModuleScore(
  your_sc_data,
  features = list(PMNs03_genes),
  name     = "PMNs03_marker_genes"
)

## ── 2.2 Create UMAP plots ───────────────────────────────

# Cell type UMAP
p_umap_type <- DimPlot(
  your_sc_data,
  group.by   = "Subset",
  label      = TRUE,
  label.size = 4,
  repel      = TRUE
) +
  ggtitle("Cell Type UMAP") +
  NoLegend() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 14)
  )

# PMNs_S03 module score UMAP
p_umap_score <- FeaturePlot(
  your_sc_data,
  features   = "PMNs03_marker_genes1",
  reduction  = "umap",
  min.cutoff = 0
) +
  ggtitle("PMNs_S03 Marker Genes") +
  theme(
    plot.title   = element_text(size = 20, face = "bold"),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    legend.text  = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

combined_plot <- p_umap_type / p_umap_score

ggsave(combined_plot,
       file   = "../output/04_PMN_validation_combined.pdf",
       width  = 10,
       height = 16,
       dpi    = 300)


# ════════════════════════════════════════════════════════════
# PART 3: PMNs_S03 GENES IN MACROPHAGE, FIBROBLAST & NEUTROPHIL
# ════════════════════════════════════════════════════════════

## ── 3.1 Define cell types and extract expression ─────────
# NOTE: Adjust these names to match your annotation labels
cell_types_of_interest <- c("Monocyte / Macrophage",
                            "Fibroblast",
                            "Neutrophil")

Idents(your_sc_data) <- your_sc_data$Subset

exp_data <- FetchData(
  your_sc_data,
  vars  = c(PMNs03_genes, "Subset"),
  layer = "data"
)

## ── 3.2 Calculate statistics per cell type ───────────────
gene_stats <- exp_data %>%
  filter(Subset %in% cell_types_of_interest) %>%
  pivot_longer(cols = -Subset, names_to = "Gene",
               values_to = "Expression") %>%
  group_by(Subset, Gene) %>%
  summarise(
    mean_exp = mean(Expression, na.rm = TRUE),
    pct_exp  = sum(Expression > 0) / n() * 100,
    .groups  = "drop"
  )

## ── 3.3 Define thresholds and filter genes ───────────────
mean_threshold <- 0.1  # minimum mean expression
pct_threshold  <- 1    # minimum % expressing cells

expressed_genes <- gene_stats %>%
  filter(mean_exp > mean_threshold & pct_exp > pct_threshold) %>%
  select(Subset, Gene) %>%
  distinct()

# Genes expressed in all three cell types
common_genes <- expressed_genes %>%
  group_by(Gene) %>%
  filter(n() == 3) %>%
  pull(Gene) %>%
  unique()

# Neutrophil-specific genes
neutrophil_specific <- expressed_genes %>%
  group_by(Gene) %>%
  filter(n() == 1 & Subset == "Neutrophil") %>%
  pull(Gene) %>%
  unique()

cat("\nCommon genes (all 3 types):", length(common_genes), "\n")
cat("Neutrophil-specific genes:", length(neutrophil_specific), "\n\n")

## ── 3.4 Visualize gene expression patterns ───────────────

if (length(common_genes) > 0) {
  p_common <- DotPlot(
    your_sc_data,
    features = common_genes,
    group.by = "Subset",
    idents   = cell_types_of_interest
  ) +
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title  = element_text(size = 16, face = "bold")
    ) +
    ggtitle("Common Genes Across 3 Cell Types")
  
  ggsave(p_common,
         file   = "../output/04_PMN_common_genes_dotplot.pdf",
         width  = max(8, length(common_genes) * 0.3 + 4),
         height = 8,
         dpi    = 300)
}

if (length(neutrophil_specific) > 0) {
  p_neutrophil <- DotPlot(
    your_sc_data,
    features = neutrophil_specific,
    group.by = "Subset",
    idents   = cell_types_of_interest
  ) +
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      plot.title  = element_text(size = 16, face = "bold")
    ) +
    ggtitle("Neutrophil-Specific PMNs_S03 Marker Genes")
  
  ggsave(p_neutrophil,
         file   = "../output/04_PMN_neutrophil_specific_dotplot.pdf",
         width  = max(8, length(neutrophil_specific) * 0.3 + 4),
         height = 8,
         dpi    = 300)
}

cat("── Done. Outputs saved to ../output/ ──────────\n")

# ============================================================
# END OF SCRIPT
# ============================================================