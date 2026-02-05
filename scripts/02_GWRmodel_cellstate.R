# ============================================================
# 02_GWRmodel_cellstate.R
#
# Purpose:
#   Perform Geographically Weighted Regression (GWR) analysis to identify
#   cell states with spatially varying associations with TLS presence.
#   Focus on macrophage, T cell, and plasma cell states that show
#   significant local correlations with TLS organization.
#
#   This script includes:
#   1) Variable screening via univariate/multivariate logistic regression
#   2) GWR model fitting with optimized bandwidth selection
#   3) Identification of significant cell states (Macrophages S01/S08, 
#      CD8+ T S01, CD4+ T S02, Plasma Cells S01)
#   4) Spatial distribution analysis of key cell states
#   5) Functional characterization via GO enrichment analysis
#
# Input:
#   (1) combined_obj (R object loaded in environment):
#       A named list of Visium samples with EcoTyper visium mode cell state output.
#       Each element must contain:
#         - $cellstate_raw / $cellstate_norm(normalized)
#             (ID, pixel_x, pixel_y, Label [TLS / NonTLS], cell state abundances)
#         - $hs
#             (sparse matrix of gene expression for each spot)
#
#   (2) all_marker_gene (data frame):
#       Marker genes for each cell state defined by EcoTyper
#       Columns: cell_state, Gene
#
# Output:
#   Five main figures:
#   
#   (1) 02_Impact_of_cell_states.pdf
#       - Barplot of cell state impact on TLS (median |β| and fraction significant spots)
#   
#   (2) 02_macrophage_spatial_associations.pdf
#       - Macrophages S01 vs S08 local regression coefficients across samples
#       - Shows consistent positive association of S01 with TLS
#   
#   (3) 02_Tcell_Plasmacell_analysis.pdf
#       - Violin plots of local coefficient distributions (CD8+ T S01, CD4+ T S02, PC S01)
#       - Statistical significance vs effect size scatter plots
#   
#   (4) 02_spatial_distribution.pdf
#       - Spatial visualization of CD4+ T S02, CD8+ T S01, PC S01, and cancer cells
#       - B26 sample with TLS regions outlined in red
#   
#   (5) 02_go_plot.pdf
#       - GO enrichment analysis (Biological Process) for high-abundance spot DEGs
#       - Functional characterization of CD8+ T S01, CD4+ T S02, and PC S01
#
# Dependencies:
#   dplyr, tidyr, ggplot2, broom, purrr, GWmodel, sp, scales,
#   ggrepel, patchwork, viridis, limma, clusterProfiler, org.Hs.eg.db
#
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(purrr)
library(GWmodel)
library(sp)
library(scales)
library(ggrepel)
library(patchwork)
library(viridis)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 1: VARIABLE SCREENING AND GWR MODEL FITTING
# ═══════════════════════════════════════════════════════════════════════════
load("../data/example_data/combined_obj.RData")
## 1.1 Prepare combined dataset from all samples ----
df_all <- purrr::map_dfr(names(combined_obj), function(samp) {
  cs <- combined_obj[[samp]]$cellstate_norm
  if (is.null(cs)) return(NULL)
  
  cs %>%
    dplyr::select(ID, pixel_x, pixel_y, Label, matches("_S[0-9]+$")) %>%
    mutate(
      Sample = samp,
      y = ifelse(Label == "TLS", 1, 0)
    )
})

## 1.2 Univariate logistic regression screening ----
uni_or <- df_all %>%
  pivot_longer(
    cols = matches("_S[0-9]+$"),
    names_to = "State",
    values_to = "Abundance"
  ) %>%
  group_by(State) %>%
  do(tidy(glm(y ~ Abundance, data = ., family = "binomial"))) %>%
  filter(term == "Abundance") %>%
  mutate(
    OR = exp(estimate),
    lower = exp(estimate - 1.96 * std.error),
    upper = exp(estimate + 1.96 * std.error),
    pval = p.value
  ) %>%
  arrange(desc(abs(estimate)))

# Select top candidates based on univariate effect size
top_candidates <- uni_or$State[1:30]

## 1.3 Multivariate logistic regression refinement ----
multi_fit <- glm(
  formula = as.formula(paste("y ~", paste(top_candidates, collapse = " + "))),
  data = df_all,
  family = "binomial"
)

multi_or <- tidy(multi_fit) %>%
  filter(term %in% top_candidates) %>%
  mutate(
    OR = exp(estimate),
    lower = exp(estimate - 1.96 * std.error),
    upper = exp(estimate + 1.96 * std.error),
    pval = p.value
  ) %>%
  arrange(desc(abs(estimate)))

# Select final states with p < 0.05 in multivariate model
final_states <- multi_or %>%
  filter(pval < 0.05) %>%
  pull(term)

cat("Selected states for GWR:", paste(final_states, collapse = ", "), "\n")

## 1.5 Fit GWR models with selected states ----

gwr_cand <- list()
###This process takes long time, you can load the prepared gwr_results for following analysis.
###load("../data/example_data/gwr_results.RData")
for (samp in names(combined_obj)) {
  cs <- combined_obj[[samp]]$cellstate_norm %>% na.omit()
  if (!"TLS" %in% cs$Label) next
  
  # Prepare response variable
  cs$y <- ifelse(cs$Label == "TLS", 1, 0)
  expl <- intersect(final_states, colnames(cs))
  if (length(expl) == 0) next
  
  # Create spatial data
  fml <- as.formula(paste("y ~", paste(expl, collapse = " + ")))
  coords <- cs[, c("pixel_x", "pixel_y")]
  rownames(cs) <- cs$ID
  spdf <- SpatialPointsDataFrame(coords, data = cs, proj4string = CRS(NA_character_))
  
  # Optimize bandwidth
  bw0 <- bw.gwr(
    formula = fml,
    data = spdf,
    approach = "AICc",
    kernel = "bisquare",
    adaptive = FALSE
  )
  
  # Fit GWR model
  gwr0 <- gwr.basic(
    formula = fml,
    data = spdf,
    bw = bw0,
    kernel = "bisquare",
    adaptive = FALSE
  )
  
  gwr_cand[[samp]] <- list(
    bandwidth = bw0,
    model = gwr0
  )
  
  message(sprintf("GWR fitted for %s (bandwidth = %.1f)", samp, bw0))
}


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 2: GWR RESULTS ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

## 2.1 Extract coefficients from GWR models ----
coef_df <- purrr::imap_dfr(gwr_cand, ~{
  tibble::as_tibble(.x$model$SDF@data, rownames = "ID") %>%
    mutate(Sample = .y)
})
total_spots <- nrow(coef_df)

## 2.2 Reshape to long format ----
coef_long <- coef_df %>%
  pivot_longer(
    cols = -c(ID, Sample, yhat, residual),
    names_to = c("Variable", "Metric"),
    names_pattern = "(.+?)(?:_(SE|TV))?$",
    values_to = "value"
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "SE" ~ "SE",
      Metric == "TV" ~ "TV",
      TRUE ~ "coef"
    )
  ) %>%
  pivot_wider(
    names_from = Metric,
    values_from = value
  )

## 2.3 Calculate SE Q1 threshold for each variable ----
se_q1 <- coef_long %>%
  group_by(Sample, Variable) %>%
  summarise(
    SE_Q1 = quantile(SE, 0.25, na.rm = TRUE),
    .groups = "drop"
  )

## 2.4 Filter significant and high-precision spots ----
sig_coef <- coef_long %>%
  filter(abs(TV) > 2) %>%                          # Significant spots
  inner_join(se_q1, by = c("Sample", "Variable")) %>%
  filter(SE <= SE_Q1)                              # High precision

## 2.5 Variable summary with directionality ----
var_summary <- sig_coef %>%
  group_by(Variable) %>%
  summarise(
    med_coef = median(coef, na.rm = TRUE),         # Median coefficient
    med_abs_coef = median(abs(coef), na.rm = TRUE), # Effect magnitude
    prop_pos = sum(coef > 0) / total_spots,         # Proportion positive
    prop_neg = sum(coef < 0) / total_spots,         # Proportion negative
    prop_sig = n() / total_spots,                   # Proportion significant
    .groups = "drop"
  ) %>%
  arrange(desc(med_abs_coef))

print(var_summary)

## 2.6 Visualization: Impact of cell states on TLS ----
plot_df_impact <- var_summary %>%
  filter(Variable != "Intercept") %>%
  mutate(
    pct_label = percent(prop_sig, accuracy = 1),
    label_x = med_coef + if_else(med_coef > 0, 0.02, -0.02),
    hjust_val = if_else(med_coef > 0, 0, 1)
  )

p_impact <- ggplot(plot_df_impact, aes(
  x = med_coef,
  y = reorder(Variable, med_coef),
  fill = prop_sig
)) +
  geom_col(width = 0.7, color = "gray30") +
  scale_fill_gradient(
    low = "lightblue",
    high = "darkblue",
    name = "Fraction significant spots"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_text(aes(
    x = label_x,
    label = pct_label,
    hjust = hjust_val
  ), size = 3) +
  labs(
    title = "Impact of Cell States on TLS Presence",
    x = "Impact degree (median |β|)",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )

suppressWarnings(
ggsave(plot = p_impact,file = "../output/02_Impact_of_cell_states.pdf"))


# ═══════════════════════════════════════════════════════════════════════════
# SECTION 3: MACROPHAGE CELL STATE ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

## 3.1 Macrophage S01 vs S08 coefficients across samples ----
sig_macro_plot <- sig_coef %>%
  filter(Variable %in% c("Monocytes.and.Macrophages_S01", 
                         "Monocytes.and.Macrophages_S08")) %>%
  ggplot(aes(x = Sample, y = coef, color = Variable, fill = Variable)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.2, color = "gray70") +
  stat_summary(fun = mean, geom = "col", 
               position = position_dodge(width = 0.8),
               alpha = 0.7, width = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 3, shape = 16,
               position = position_dodge(width = 0.8),
               color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_fill_manual(
    values = c("Monocytes.and.Macrophages_S01" = "lightblue",
               "Monocytes.and.Macrophages_S08" = "lightcoral"),
    labels = c("Macrophages S01", "Macrophages S08"),
    name = "Cell State"
  ) +
  scale_color_manual(
    values = c("Monocytes.and.Macrophages_S01" = "lightblue",
               "Monocytes.and.Macrophages_S08" = "lightcoral"),
    guide = "none"
  ) +
  labs(
    title = "Spatial Associations of Macrophage Cell States with TLS",
    subtitle = "Significant spots only (|t-value| > 2 & SE ≤ Q1)",
    x = "Sample ID",
    y = "Local Regression Coefficient",
    caption = "Bars indicate mean values per sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray60"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(sig_macro_plot)
suppressWarnings(
  ggsave(plot = sig_macro_plot, file = "../output/02_macrophage_spatial_associations.pdf"))




# ═══════════════════════════════════════════════════════════════════════════
# SECTION 4: T CELL AND PLASMA CELL ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════

## 4.1 Overall statistics for key cell states ----
target_cells <- c("PCs_S01", "CD8.T_S01", "CD4.T_S02")

overall_stats <- sig_coef %>%
  filter(Variable %in% target_cells) %>%
  group_by(Variable) %>%
  summarise(
    total_sig_spots = n(),
    median_coef = median(coef, na.rm = TRUE),
    median_abs_coef = median(abs(coef), na.rm = TRUE),
    q25_coef = quantile(coef, 0.25, na.rm = TRUE),
    q75_coef = quantile(coef, 0.75, na.rm = TRUE),
    prop_positive = sum(coef > 0) / n(),
    mean_tv = mean(abs(TV), na.rm = TRUE),
    prop_very_sig = sum(abs(TV) > 3) / n(),
    .groups = "drop"
  )

print("Overall Statistics for Key Cell States:")
print(overall_stats)

## 4.2 Coefficient distribution (violin plot) ----
plot1 <- sig_coef %>%
  filter(Variable %in% target_cells) %>%
  mutate(Variable_clean = str_replace_all(Variable, "[._]", " ")) %>%
  ggplot(aes(x = Variable_clean, y = coef, fill = Variable)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(
    values = c("PCs_S01" = "#E31A1C", 
               "CD8.T_S01" = "#1F78B4", 
               "CD4.T_S02" = "#33A02C")
  ) +
  labs(
    title = "Distribution of Local Coefficients for Top Cell States",
    x = "Cell State",
    y = expression(bold("Local Coefficient (" * beta * ")"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

## 4.3 Statistical significance vs effect size ----
plot2 <- sig_coef %>%
  filter(Variable %in% target_cells) %>%
  mutate(Variable_clean = str_replace_all(Variable, "[._]", " ")) %>%
  ggplot(aes(x = coef, y = abs(TV), color = Variable)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "gray60", alpha = 0.8) +
  geom_hline(yintercept = 2.5, linetype = "dotted", color = "gray60", alpha = 0.8) +
  geom_hline(yintercept = 3, linetype = "solid", color = "gray60", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("PCs_S01" = "#E31A1C", 
               "CD8.T_S01" = "#1F78B4", 
               "CD4.T_S02" = "#33A02C"),
    name = "Cell State"
  ) +
  facet_wrap(~Variable_clean, scales = "free") +
  labs(
    title = "Statistical Significance vs. Effect Size",
    x = expression(bold("Local Coefficient (" * beta * ")")),
    y = "|t-value| (significance)",
    caption = "Horizontal lines: |t| = 2 (dashed), 2.5 (dotted), 3 (solid)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

## 4.4 Combined plot ----
combined_tcell_plot <- plot1 / plot2 +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Impact Analysis of Top Cell States on TLS",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.margin = margin(20, 20, 20, 20)
    )
  )

print(combined_tcell_plot)
suppressWarnings(
  ggsave(plot = combined_tcell_plot, file = "../output/02_Tcell_Plasmacell_analysis.pdf")
)

# ═══════════════════════════════════════════════════════════════════════════
# SECTION 5: SPATIAL DISTRIBUTION IN B26 SAMPLE (CD4S02,CD8S01,PCS01,cancercell)
# ═══════════════════════════════════════════════════════════════════════════

## 5.1 Prepare B26 data ----
# Extract cell states
df26 <- combined_obj$B26$cellstate_raw %>%
  dplyr::select(ID, pixel_x, pixel_y,
         PCs_S01, CD8.T_S01, CD4.T_S02,
         Monocytes.and.Macrophages_S01,
         Monocytes.and.Macrophages_S08,
         Label) %>%
  dplyr::rename(
    PC_S01 = PCs_S01,
    CD8_S01 = CD8.T_S01,
    CD4_S02 = CD4.T_S02,
    Macro_S01 = Monocytes.and.Macrophages_S01,
    Macro_S08 = Monocytes.and.Macrophages_S08
  ) %>%
  mutate(isTLS = (Label == "TLS"))

# Calculate cancer cell score
epithelial_markers <- c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1")
available_markers <- epithelial_markers[epithelial_markers %in% rownames(combined_obj$B26$hs)]
cat("Available epithelial markers:", paste(available_markers, collapse = ", "), "\n")

if (length(available_markers) > 0) {
  expr_matrix <- combined_obj$B26$hs[available_markers, , drop = FALSE]
  cancer_scores <- Matrix::colMeans(expr_matrix)
  df26$Cancer_score <- cancer_scores[df26$ID]
}

## 5.2 Spatial plotting function ----
plot_spatial <- function(var, title, color_low = "lightgrey", 
                         color_high = "darkblue", legend_name = NULL) {
  if (is.null(legend_name)) legend_name <- title
  
  ggplot(df26, aes(x = pixel_x, y = -pixel_y)) +
    geom_point(aes(color = .data[[var]]), size = 1.2, alpha = 0.9) +
    geom_point(data = filter(df26, isTLS),
               shape = 21, fill = NA, color = "red", 
               size = 1.2, stroke = 0.6) +
    scale_color_gradient(
      low = color_low, 
      high = color_high, 
      name = legend_name,
      guide = guide_colorbar(
        title.theme = element_text(size = 14, face = "bold"),
        label.theme = element_text(size = 12),
        barwidth = 1.5,
        barheight = 8
      )
    ) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
}

## 5.3 Create spatial panels ----
p_cd4 <- plot_spatial("CD4_S02", "CD4_S02", "lightgrey", "forestgreen")
p_cd8 <- plot_spatial("CD8_S01", "CD8_S01", "lightgrey", "blue")
p_pc <- plot_spatial("PC_S01", "PC_S01", "lightgrey", "red")
p_cancer <- plot_spatial("Cancer_score", "Cancer Cells", 
                         "lightgrey", "darkorange",
                         legend_name = "Cancer Cell Score")

# Combined spatial plot
spatial_combined <- (p_cd4 | p_cd8 | p_pc | p_cancer) +
  plot_annotation(
    title = "Spatial Distribution of Key Cell States and Cancer Cells",
    subtitle = "TLS regions outlined in red",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray30")
    )
  ) &
  theme(plot.margin = margin(5, 5, 5, 5))

print(spatial_combined)
suppressWarnings(
  ggsave(plot = spatial_combined, file = "../output/02_spatial_distribution.pdf")
)

# ═══════════════════════════════════════════════════════════════════════════
# SECTION 6: FUNCTIONAL ANALYSIS (GO ENRICHMENT)
# ═══════════════════════════════════════════════════════════════════════════

## 6.1 Prepare expression data ----
expr <- as.matrix(combined_obj$B26$hs)

# Identify high abundance spots (top 25%)
cs <- combined_obj$B26$cellstate_raw %>%
  dplyr::select(ID, CD8.T_S01, PCs_S01, CD4.T_S02) %>%
  mutate(
    high_CD8 = CD8.T_S01 >= quantile(CD8.T_S01, 0.75, na.rm = TRUE),
    high_PC = PCs_S01 >= quantile(PCs_S01, 0.75, na.rm = TRUE),
    high_CD4 = CD4.T_S02 >= quantile(CD4.T_S02, 0.75, na.rm = TRUE)
  )

# Match expression matrix to cell state data
expr <- expr[, cs$ID]

## 6.2 Convert marker genes to ENTREZID ----
#marker gene defined by EcoTyper
all_marker_gene <- readRDS("~/Result3/gene_info/all_marker_gene.rds")
get_markers <- function(cellstate_name) {
  markers <- all_marker_gene %>%
    dplyr::filter(cell_state == cellstate_name) %>%
    dplyr::pull(Gene)
  
  converted <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = markers,
    keytype  = "SYMBOL",
    columns  = "ENTREZID"
  )
  converted <- converted[!is.na(converted$ENTREZID), ]
  
  list(
    symbols     = markers,
    entrez      = unique(converted$ENTREZID),
    n_original  = length(markers),
    n_converted = nrow(converted)
  )
}
CD8_markers <- get_markers("CD8.T_S01")
PC_markers <- get_markers("PCs_S01")
CD4_markers <- get_markers("CD4.T_S02")

## 6.3 GO enrichment for marker genes ----
perform_go_analysis <- function(entrez_ids, description) {
  if (length(entrez_ids) < 5) {
    warning(paste("Too few genes for", description))
    return(NULL)
  }
  
  result <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
  
  if (is.null(result) || nrow(result@result) == 0) {
    warning(paste("No significant GO terms for", description))
    return(NULL)
  }
  
  return(result)
}

go_CD8_markers <- perform_go_analysis(CD8_markers$entrez, "CD8_S01 markers")
go_PC_markers <- perform_go_analysis(PC_markers$entrez, "PCs_S01 markers")
go_CD4_markers <- perform_go_analysis(CD4_markers$entrez, "CD4_S02 markers")

## 6.4 Differential expression in high abundance spots ----
perform_de_analysis <- function(group_var, group_name = NULL,
                                logFC_cutoff = 1,
                                padj_cutoff = 0.05) {
  design <- model.matrix(as.formula(paste0("~ ", group_var)), data = cs)
  fit <- lmFit(expr, design) %>% eBayes()
  
  coef_name <- paste0(group_var, "TRUE")
  top_table <- topTable(fit, coef = coef_name, number = Inf)
  
  # Extract significant DEGs (SYMBOL)
  de_genes <- top_table %>%
    tibble::rownames_to_column("Gene") %>%
    dplyr::filter(
      logFC > logFC_cutoff,
      adj.P.Val < padj_cutoff
    ) %>%
    dplyr::arrange(desc(logFC))
  
  list(
    top_table = top_table,
    de_genes  = de_genes,
    symbols   = de_genes$Gene,
    n_de      = nrow(de_genes),
    group     = group_name %||% group_var
  )
}

de_CD8 <- perform_de_analysis("high_CD8", "CD8_S01")
de_PC <- perform_de_analysis("high_PC", "PCs_S01")
de_CD4 <- perform_de_analysis("high_CD4", "CD4_S02")

cat("DE analysis summary:\n")
cat(sprintf("CD8_S01 high spots: %d DEGs\n", de_CD8$n_de))
cat(sprintf("PCs_S01 high spots: %d DEGs\n", de_PC$n_de))
cat(sprintf("CD4_S02 high spots: %d DEGs\n", de_CD4$n_de))

## 6.5 GO enrichment for DEGs ----
go_CD8_high <- perform_go_analysis(de_CD8$entrez_ids, "CD8_S01 high spots")
go_PC_high <- perform_go_analysis(de_PC$entrez_ids, "PCs_S01 high spots")
go_CD4_high <- perform_go_analysis(de_CD4$entrez_ids, "CD4_S02 high spots")

## 6.6 Visualize GO results ----
p4 <- if (!is.null(go_CD8_high)) {
  dotplot(go_CD8_high, showCategory = 8) + 
    ggtitle("CD8T S01 High Spots GO BP") + 
    theme(axis.text.y = element_text(size = 10))
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5, 
                      label = "No significant GO terms", size = 6)
}

p5 <- if (!is.null(go_PC_high)) {
  dotplot(go_PC_high, showCategory = 8) + 
    ggtitle("Plasma Cells S01 High Spots GO BP") + 
    theme(axis.text.y = element_text(size = 10))
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5, 
                      label = "No significant GO terms", size = 6)
}

p6 <- if (!is.null(go_CD4_high)) {
  dotplot(go_CD4_high, showCategory = 8) + 
    ggtitle("CD4T S02 High Spots GO BP") + 
    theme(axis.text.y = element_text(size = 10))
} else {
  ggplot() + annotate("text", x = 0.5, y = 0.5, 
                      label = "No significant GO terms", size = 6)
}

combined_go_plot <- (p4 | p5 | p6) +
  plot_annotation(
    title = "Functional Analysis of Key Cell States in TLS Formation",
    subtitle = "GO analysis of high-abundance spot DEGs",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30"),
      plot.margin = margin(20, 20, 20, 20)
    )
  )

print(combined_go_plot)
suppressWarnings(
  ggsave(plot = combined_go_plot, file = "../output/02_go_plot.pdf")
)

# ═══════════════════════════════════════════════════════════════════════════
# END OF SCRIPT
# ═══════════════════════════════════════════════════════════════════════════