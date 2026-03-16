# ============================================================
# 03_Abundance_vs_Centrality.R
#
# Purpose:
#   Test whether chemokine network hub status is simply
#   driven by cell state abundance in TLS, or reflects
#   genuine interaction specificity.
#
# Approach:
#   x-axis: mean proportional abundance in TLS spots
#   y-axis: median weighted degree across samples
#   Key expectation: if hub = abundance artifact, points
#   fall along diagonal. PMNs_S03 (low abundance, high
#   centrality) should be a clear outlier.
#
# Input:
#   combined_obj.RData: Visium samples
#   LR_analysis.RData: chemokine network results
#     (output of analyze_LR_network, one list per axis)
#
# Output:
#   03b_abundance_vs_centrality.pdf
#   03b_abundance_vs_centrality_per_axis.pdf
#
# Dependencies:
#   dplyr, tidyr, ggplot2, ggrepel, purrr
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(purrr)

load("../data/example_data/combined_obj.RData")
load("../data/example_data/LR_analysis.RData")

# ════════════════════════════════════════════════════════════
# SECTION 1: MEAN ABUNDANCE IN TLS SPOTS
# ════════════════════════════════════════════════════════════

# Collect all TLS spots across samples
df_tls <- purrr::map_dfr(names(combined_obj), function(samp) {
  meta <- combined_obj[[samp]]$cellstate_norm
  if (is.null(meta)) return(NULL)
  
  meta %>%
    dplyr::filter(Label == "TLS") %>%
    dplyr::select(ID, matches("_S[0-9]+$")) %>%
    mutate(Sample = samp)
})

# Mean abundance per cell state across all TLS spots
cell_states <- df_tls %>%
  dplyr::select(matches("_S[0-9]+$")) %>%
  colnames()

abundance_summary <- df_tls %>%
  dplyr::select(all_of(cell_states)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(),
               names_to  = "cell_state",
               values_to = "mean_tls_abundance")

# ════════════════════════════════════════════════════════════
# SECTION 2: WEIGHTED DEGREE FROM NETWORK SCORE MATRICES
# ════════════════════════════════════════════════════════════

axis_results <- list(
  CXCL13   = LR_analysis$CXCL13_CXCR5,
  CXCL12   = LR_analysis$CXCL12_CXCR4,
  CCL19_21 = LR_analysis$CCL19_21_CCR7
)

# Compute weighted degree per cell state per sample per axis
# Weighted degree = row sum + col sum of score matrix
#   (total send + receive interaction strength)
compute_weighted_degree <- function(results_list) {
  purrr::map_dfr(names(results_list), function(sid) {
    res <- results_list[[sid]]
    if (is.null(res) || is.null(res$score)) return(NULL)
    
    score <- res$score
    cs <- rownames(score)
    
    w_out <- rowSums(score)   # send
    w_in  <- colSums(score)   # receive
    w_deg <- w_out + w_in
    
    data.frame(
      Sample     = sid,
      cell_state = cs,
      w_degree   = w_deg,
      stringsAsFactors = FALSE
    )
  })
}

# Process each axis
centrality_all <- purrr::map_dfr(names(axis_results), function(axis) {
  compute_weighted_degree(axis_results[[axis]]) %>%
    mutate(axis = axis)
})

# Median weighted degree across samples, per axis
centrality_median <- centrality_all %>%
  group_by(cell_state, axis) %>%
  summarise(
    median_w_degree = median(w_degree, na.rm = TRUE),
    n_samples       = n(),
    .groups         = "drop"
  )

# Average median weighted degree across the three axes
centrality_overall <- centrality_median %>%
  group_by(cell_state) %>%
  summarise(
    mean_median_w_degree = mean(median_w_degree, na.rm = TRUE),
    .groups = "drop"
  )

# ════════════════════════════════════════════════════════════
# SECTION 3: MERGE AND COMPUTE RANKS
# ════════════════════════════════════════════════════════════

df_plot <- abundance_summary %>%
  inner_join(centrality_overall, by = "cell_state") %>%
  mutate(
    rank_abundance  = rank(-mean_tls_abundance),
    rank_centrality = rank(-mean_median_w_degree),
    highlight = case_when(
      cell_state == "PMNs_S03" ~ "PMNs_S03",
      cell_state %in% c("B_S01", "PCs_S01", "CD4.T_S02",
                        "Dendritic_S01", "CD8.T_S01") ~ "C1 Hub",
      TRUE ~ "Other"
    )
  )

# Spearman correlation between abundance and centrality
cor_test <- cor.test(df_plot$mean_tls_abundance,
                     df_plot$mean_median_w_degree,
                     method = "spearman")

cat("── Abundance vs Centrality ──────────────────────\n")
cat("Spearman rho:", round(cor_test$estimate, 3), "\n")
cat("p-value     :", signif(cor_test$p.value, 3), "\n\n")

# Flag PMNs_S03 specifically
pmn <- df_plot %>% filter(cell_state == "PMNs_S03")
cat("── PMNs_S03 ────────────────────────────────────\n")
cat("Abundance rank :", pmn$rank_abundance, "/", nrow(df_plot), "\n")
cat("Centrality rank:", pmn$rank_centrality, "/", nrow(df_plot), "\n\n")

# ════════════════════════════════════════════════════════════
# SECTION 4: SCATTER PLOT
# ════════════════════════════════════════════════════════════

df_plot <- df_plot %>%
  mutate(
    abundance_z  = as.numeric(scale(mean_tls_abundance)),
    centrality_z = as.numeric(scale(mean_median_w_degree))
  )

color_map <- c("PMNs_S03" = "#E64B35",
               "C1 Hub"   = "#4DBBD5",
               "Other"    = "gray60")

p_scatter <- ggplot(df_plot,
                    aes(x = abundance_z, y = centrality_z)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_point(aes(color = highlight, size = highlight), alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = df_plot %>% filter(highlight != "Other"),
    aes(label = cell_state),
    size          = 3.5,
    box.padding   = 0.6,
    point.padding = 0.3,
    max.overlaps  = 20,
    segment.color = "gray60"
  ) +
  scale_color_manual(values = color_map) +
  scale_size_manual(values = c("PMNs_S03" = 4.5,
                               "C1 Hub"   = 3.5,
                               "Other"    = 2),
                    guide = "none") +
  annotate("text",
           x = max(df_plot$abundance_z) * 0.6,
           y = min(df_plot$centrality_z) * 0.8,
           label = paste0("Spearman rho = ", round(cor_test$estimate, 2),
                          "\np = ", signif(cor_test$p.value, 2)),
           color = "gray40", size = 3.5, hjust = 0) +
  annotate("text",
           x = min(df_plot$abundance_z) * 0.8,
           y = max(df_plot$centrality_z) * 0.9,
           label = "Low abundance\nHigh centrality",
           color = "#E64B35", size = 3, fontface = "italic",
           alpha = 0.7) +
  labs(
    title    = "TLS Abundance vs. Network Centrality",
    subtitle = "Testing whether hub status is driven by abundance",
    x        = "Mean TLS Abundance (z-score)",
    y        = "Median Weighted Degree (z-score)",
    color    = NULL,
    caption  = paste0("Weighted degree averaged across CXCL13, CXCL12, ",
                      "CCL19/21 axes\n",
                      "Dashed line: y = x (abundance-driven expectation)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = c(0.85, 0.15),
    legend.background = element_rect(fill = "white", color = NA),
    plot.caption     = element_text(color = "gray50", size = 9)
  )

ggsave(plot   = p_scatter,
       file   = "../output/03b_abundance_vs_centrality.pdf",
       width  = 8,
       height = 7,
       dpi    = 300)

# ════════════════════════════════════════════════════════════
# SECTION 5: PER-AXIS SCATTER (SUPPLEMENTARY)
# ════════════════════════════════════════════════════════════

df_per_axis <- abundance_summary %>%
  inner_join(centrality_median, by = "cell_state") %>%
  mutate(
    highlight = case_when(
      cell_state == "PMNs_S03" ~ "PMNs_S03",
      cell_state %in% c("B_S01", "PCs_S01", "CD4.T_S02",
                        "Dendritic_S01", "CD8.T_S01") ~ "C1 Hub",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(axis) %>%
  mutate(
    abundance_z  = as.numeric(scale(mean_tls_abundance)),
    centrality_z = as.numeric(scale(median_w_degree))
  ) %>%
  ungroup()

p_per_axis <- ggplot(df_per_axis,
                     aes(x = abundance_z, y = centrality_z)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray70", linewidth = 0.5) +
  geom_point(aes(color = highlight, size = highlight), alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = df_per_axis %>% filter(highlight != "Other"),
    aes(label = cell_state),
    size = 3, box.padding = 0.5, max.overlaps = 15,
    segment.color = "gray60"
  ) +
  scale_color_manual(values = color_map) +
  scale_size_manual(values = c("PMNs_S03" = 4, "C1 Hub" = 3, "Other" = 1.8),
                    guide = "none") +
  facet_wrap(~ axis, scales = "free") +
  labs(
    title = "Abundance vs. Centrality by Chemokine Axis",
    x     = "Mean TLS Abundance (z-score)",
    y     = "Median Weighted Degree (z-score)",
    color = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    strip.text       = element_text(face = "bold", size = 11)
  )

ggsave(plot   = p_per_axis,
       file   = "../output/03b_abundance_vs_centrality_per_axis.pdf",
       width  = 13,
       height = 5.5,
       dpi    = 300)

cat("── Done ─────────────────────────────────────────\n")