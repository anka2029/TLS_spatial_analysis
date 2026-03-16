# ============================================================
# 02_CellState_DistanceGradient.R
#
# Purpose:
#   Analyze spatial distribution of TLS-associated cell states
#   relative to TLS centroids, complementing the logistic
#   regression analysis (02_LogisticModel_cellstate.R).
#
# Approach:
#   1) Identify discrete TLS regions per sample via connected
#      component analysis on TLS spot adjacency graph
#   2) Compute geometric centroid of each TLS region
#   3) For each spot (TLS + non-TLS), calculate Euclidean
#      distance to nearest TLS centroid (array units × 100 μm)
#   4) Bin distances using pooled quartile thresholds
#      (consistent with CE distance analysis, Section S2)
#   5) Visualize abundance gradient for top associated states
#   6) Cell state composition by cell type
#   7) Hub community (C1/C2) spatial gradient from TLS centroid (SFig.3c)
#
# Input:
#   combined_obj: named list of Visium samples
#   multi_results: from 02_Logistic.R
#
# Output:
#   02b_cellstate_distance_gradient.pdf
#   02b_cellstate_composition.pdf
#   SFig.3c_hub_community_gradient.pdf
#
# Dependencies:
#   igraph, dplyr, tidyr, ggplot2, patchwork, stringr, purrr
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(igraph)
library(stringr)
library(purrr)
library(conflicted)

set.seed(42)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::setdiff)
conflicts_prefer(base::intersect)

load("../data/example_data/combined_obj.RData")
load("../data/example_data/multi_regression.RData")

# ════════════════════════════════════════════════════════════
# SECTION 0: DATA PREPARATION
# ════════════════════════════════════════════════════════════

# --- Build spot-level data frame (needed for Section 5) ---
df_all <- purrr::map_dfr(names(combined_obj), function(samp) {
  cs <- combined_obj[[samp]]$cellstate_norm
  if (is.null(cs)) return(NULL)
  cs %>%
    dplyr::select(ID, array_x, array_y, Label,
                  matches("_S[0-9]+$")) %>%
    mutate(Sample = samp,
           y      = as.integer(Label == "TLS"))
})

tls_count <- tapply(df_all$y, df_all$Sample, sum)
samples_with_TLS <- names(tls_count[tls_count > 0])
df_all <- df_all[df_all$Sample %in% samples_with_TLS, ]

all_states <- df_all %>%
  dplyr::select(matches("_S[0-9]+$")) %>%
  colnames() %>%
  setdiff("Fibroblasts_S07")  # remove doublet

# --- Select cell states to analyze ---
# Top TLS-enriched and TLS-depleted states by |statistic|
# Exclude Fibroblasts_S07 (identified as doublet)
multi_filtered <- multi_results %>%
  dplyr::filter(term != "Fibroblasts_S07",
                !grepl("Intercept", term)) %>%
  mutate(direction = ifelse(estimate > 1, "Enriched", "Depleted"))

n_top <- 7  # top N from each direction

top_enriched <- multi_filtered %>%
  dplyr::filter(direction == "Enriched", is_interaction == "FALSE") %>%
  arrange(desc(abs(statistic))) %>%
  head(n_top) %>%
  pull(term)

top_depleted <- multi_filtered %>%
  dplyr::filter(direction == "Depleted") %>%
  arrange(desc(abs(statistic))) %>%
  head(n_top) %>%
  pull(term)

target_states <- c(top_enriched, top_depleted)

# ════════════════════════════════════════════════════════════
# SECTION 1: TLS CLUSTER IDENTIFICATION
# ════════════════════════════════════════════════════════════

# Adjacency threshold for connected component analysis
# Visium hex grid: nearest neighbors are sqrt(1) or sqrt(3) apart
# Use threshold slightly above sqrt(3) ~ 1.73 to capture hex neighbors
adj_threshold <- 2.0

compute_tls_centroids <- function(meta, adj_thr = adj_threshold) {
  
  tls_spots <- meta %>%
    filter(Label == "TLS") %>%
    filter(!is.na(array_x), !is.na(array_y))
  
  if (nrow(tls_spots) == 0) return(NULL)
  
  coords <- as.matrix(tls_spots[, c("array_x", "array_y")])
  
  dmat <- as.matrix(stats::dist(coords))
  adj_mat <- (dmat > 0 & dmat <= adj_thr) * 1
  
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  comp <- igraph::components(g)
  
  tls_spots$tls_region <- comp$membership
  
  centroids <- tls_spots %>%
    group_by(tls_region) %>%
    summarise(
      cx = mean(array_x),
      cy = mean(array_y),
      n_spots = n(),
      .groups = "drop"
    )
  
  list(centroids = centroids, tls_spots = tls_spots)
}

# ════════════════════════════════════════════════════════════
# SECTION 2: DISTANCE COMPUTATION
# ════════════════════════════════════════════════════════════

# Helper function to compute distance to nearest TLS centroid
compute_spot_distances <- function(samp, state_cols) {
  meta <- combined_obj[[samp]]$cellstate_norm
  if (is.null(meta)) return(NULL)
  if (sum(meta$Label == "TLS") == 0) return(NULL)
  
  tls_info <- compute_tls_centroids(meta)
  if (is.null(tls_info)) return(NULL)
  
  centroids <- tls_info$centroids
  
  spot_coords <- as.matrix(meta[, c("array_x", "array_y")])
  centroid_coords <- as.matrix(centroids[, c("cx", "cy")])
  
  # Distance matrix: spots x centroids
  dist_to_centroids <- sapply(seq_len(nrow(centroid_coords)), function(k) {
    sqrt((spot_coords[, 1] - centroid_coords[k, 1])^2 +
           (spot_coords[, 2] - centroid_coords[k, 2])^2)
  })
  
  if (is.null(dim(dist_to_centroids))) {
    min_dist <- dist_to_centroids
  } else {
    min_dist <- apply(dist_to_centroids, 1, min)
  }
  
  # Convert to micrometers (array unit x 100)
  dist_um <- min_dist * 100
  
  available_cols <- base::intersect(state_cols, colnames(meta))
  
  meta %>%
    dplyr::select(ID, array_x, array_y, Label,
                  all_of(available_cols)) %>%
    mutate(
      Sample = samp,
      dist_to_centroid_um = dist_um,
      is_TLS = Label == "TLS"
    )
}

# Compute distances for target cell states
df_gradient <- purrr::map_dfr(names(combined_obj), function(samp) {
  compute_spot_distances(samp, target_states)
})

cat("── Distance summary ─────────────────────────────\n")
cat("Total spots  :", nrow(df_gradient), "\n")
cat("TLS spots    :", sum(df_gradient$is_TLS), "\n")
cat("Samples      :", n_distinct(df_gradient$Sample), "\n\n")

# ════════════════════════════════════════════════════════════
# SECTION 3: DISTANCE BINNING (POOLED QUARTILE APPROACH)
# ════════════════════════════════════════════════════════════

# Compute quartiles from pooled non-TLS distance distribution
# (TLS spots form their own bin; quartiles are for non-TLS spots)
non_tls_dist <- df_gradient %>%
  filter(!is_TLS) %>%
  pull(dist_to_centroid_um)

q_vals <- quantile(non_tls_dist, probs = c(0.25, 0.5, 0.75))

cat("── Pooled distance quartiles (non-TLS, um) ─────\n")
cat("Q1 (25%):", round(q_vals[1]), "um\n")
cat("Q2 (50%):", round(q_vals[2]), "um\n")
cat("Q3 (75%):", round(q_vals[3]), "um\n\n")

# For TLS spots: subdivide by within-TLS distance to centroid
# to reveal core vs. edge structure
tls_dist <- df_gradient %>%
  filter(is_TLS) %>%
  pull(dist_to_centroid_um)

tls_median <- median(tls_dist)

assign_distance_bins <- function(df, q_vals, tls_median) {
  df %>%
    mutate(
      dist_bin = case_when(
        is_TLS & dist_to_centroid_um <= tls_median ~ "TLS Core",
        is_TLS & dist_to_centroid_um > tls_median  ~ "TLS Edge",
        !is_TLS & dist_to_centroid_um <= q_vals[1] ~ "Near",
        !is_TLS & dist_to_centroid_um <= q_vals[2] ~ "Medium",
        !is_TLS & dist_to_centroid_um <= q_vals[3] ~ "Far",
        TRUE                                        ~ "Very Far"
      ),
      dist_bin = factor(dist_bin,
                        levels = c("TLS Core", "TLS Edge",
                                   "Near", "Medium", "Far", "Very Far"))
    )
}

df_gradient <- assign_distance_bins(df_gradient, q_vals, tls_median)

cat("── Spots per bin ────────────────────────────────\n")
print(table(df_gradient$dist_bin))
cat("\n")

# ════════════════════════════════════════════════════════════
# SECTION 4: COMPUTE GRADIENT STATISTICS & VISUALIZE
# ════════════════════════════════════════════════════════════

# For each distance bin x cell state:
#   median abundance per sample -> then mean +/- SE across samples
# (Same approach as CE distance analysis)

states_present <- base::intersect(target_states, colnames(df_gradient))

gradient_stats <- df_gradient %>%
  pivot_longer(
    cols = all_of(states_present),
    names_to = "cell_state",
    values_to = "abundance"
  ) %>%
  group_by(cell_state, dist_bin, Sample) %>%
  summarise(median_abund = median(abundance, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(cell_state, dist_bin) %>%
  summarise(
    mean_abund = mean(median_abund, na.rm = TRUE),
    se_abund   = sd(median_abund, na.rm = TRUE) / sqrt(n()),
    n_samples  = n(),
    .groups    = "drop"
  )

gradient_stats <- gradient_stats %>%
  mutate(direction = ifelse(cell_state %in% top_enriched,
                            "TLS-Enriched", "TLS-Depleted"))

# Color palettes
enriched_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#7E6148", "#B09C85"
)
depleted_colors <- c(
  "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
  "#E15759", "#76B7B2"
)

# --- Plot function for one direction ---
plot_gradient <- function(data, direction_label, colors) {
  states_in_data <- unique(data$cell_state)
  color_map <- setNames(colors[seq_along(states_in_data)], states_in_data)
  
  ggplot(data, aes(x = dist_bin, y = mean_abund,
                   color = cell_state, group = cell_state)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_abund - se_abund,
                      ymax = mean_abund + se_abund),
                  width = 0.15, linewidth = 0.5) +
    geom_vline(xintercept = 2.5, linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    annotate("text", x = 1.5,
             y = max(data$mean_abund + data$se_abund) * 0.95,
             label = "TLS", color = "gray40", size = 3,
             fontface = "italic") +
    annotate("text", x = 4.5,
             y = max(data$mean_abund + data$se_abund) * 0.95,
             label = "non-TLS", color = "gray40", size = 3,
             fontface = "italic") +
    scale_color_manual(values = color_map) +
    labs(
      title = direction_label,
      x     = "Distance to TLS centroid",
      y     = "Mean proportional abundance (+/- SE)",
      color = "Cell State"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.x      = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}

p_enriched <- gradient_stats %>%
  filter(direction == "TLS-Enriched") %>%
  plot_gradient("TLS-Enriched Cell States", enriched_colors)

p_depleted <- gradient_stats %>%
  filter(direction == "TLS-Depleted") %>%
  plot_gradient("TLS-Depleted Cell States", depleted_colors)

p_combined <- p_enriched / p_depleted +
  plot_annotation(
    title   = "Cell State Abundance Gradient from TLS Centroid",
    caption = paste0("Bins: TLS Core/Edge (median split) | ",
                     "Non-TLS: pooled quartile thresholds\n",
                     "Error bars: SE across samples | ",
                     "Dashed line: TLS boundary"),
    theme   = theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.caption = element_text(color = "gray50", size = 9)
    )
  )

ggsave(plot   = p_combined,
       file   = "../output/02b_cellstate_distance_gradient.pdf",
       width  = 11,
       height = 10,
       dpi    = 300)

cat("── Gradient plot saved ─────────────────────────\n\n")


# ════════════════════════════════════════════════════════════
# SECTION 5: CELL STATE COMPOSITION BY CELL TYPE
# ════════════════════════════════════════════════════════════

abundance_summary <- df_all %>%
  dplyr::select(all_of(all_states)) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "state",
               values_to = "mean_abundance") %>%
  mutate(
    cell_type = str_replace(state, "_S[0-9]+$", "")
  ) %>%
  group_by(cell_type) %>%
  mutate(proportion = mean_abundance / sum(mean_abundance)) %>%
  ungroup()

state_classification <- abundance_summary %>%
  inner_join(
    multi_results %>%
      filter(!is_interaction) %>%
      select(term, estimate, FDR) %>%
      rename(state = term),
    by = "state"
  ) %>%
  mutate(
    direction = ifelse(log2(estimate) > 0, "TLS", "Non-TLS"),
    significance = ifelse(FDR < 0.05, "Sign.", "NS"),
    category = factor(paste(direction, significance, sep = "_"),
                      levels = c("TLS_Sign.", "TLS_NS",
                                 "Non-TLS_NS", "Non-TLS_Sign."))
  ) %>%
  group_by(cell_type) %>%
  mutate(proportion = mean_abundance / sum(mean_abundance)) %>%
  ungroup() %>%
  mutate(state_label = str_extract(state, "S[0-9]+$")) %>%
  arrange(cell_type, category)

composition_colors <- c(
  "TLS_Sign."     = "#1a5276",
  "TLS_NS"        = "#85c1e9",
  "Non-TLS_NS"    = "#f0b27a",
  "Non-TLS_Sign." = "#ba4a00"
)

p_composition <- ggplot(state_classification,
                        aes(x = cell_type, y = proportion,
                            fill = category)) +
  geom_bar(stat = "identity", position = "stack",
           color = "white", linewidth = 0.3) +
  geom_text(aes(label = state_label),
            position = position_stack(vjust = 0.5), size = 2.5) +
  scale_fill_manual(values = composition_colors) +
  coord_flip() +
  labs(title = "Cell State Composition by Cell Type",
       subtitle = "TLS association from multivariable logistic regression",
       x = NULL, y = "Relative Proportion", fill = NULL) +
  theme_minimal(base_size = 12)

ggsave(plot   = p_composition,
       file   = "../output/02b_cellstate_composition.pdf",
       width  = 10,
       height = 7,
       dpi    = 300)

cat("── Composition plot saved ──────────────────────\n\n")


# ════════════════════════════════════════════════════════════
# SECTION 6: HUB COMMUNITY (C1/C2) SPATIAL GRADIENT (SFig.3c used the same distant gradient analysis)
# ════════════════════════════════════════════════════════════

C1_members <- c("B_S01", "CD4.T_S02", "Dendritic_S01",
                "PMNs_S03", "PCs_S01")
C2_members <- c("Endothelial_S02", "Epithelial_S04",
                "Monocytes.and.Macrophages_S01", "NK_S01", "NK_S03")

all_hub_states <- c(C1_members, C2_members)

# Compute distances for hub states (reuse helper function)
df_hub_gradient <- purrr::map_dfr(names(combined_obj), function(samp) {
  compute_spot_distances(samp, all_hub_states)
})

# Recompute quartiles for hub gradient data
non_tls_dist_hub <- df_hub_gradient %>%
  filter(!is_TLS) %>%
  pull(dist_to_centroid_um)

q_vals_hub <- quantile(non_tls_dist_hub, probs = c(0.25, 0.5, 0.75))

tls_dist_hub <- df_hub_gradient %>%
  filter(is_TLS) %>%
  pull(dist_to_centroid_um)

tls_median_hub <- median(tls_dist_hub)

df_hub_gradient <- assign_distance_bins(df_hub_gradient,
                                        q_vals_hub, tls_median_hub)

# Compute gradient statistics
available_hub_states <- base::intersect(all_hub_states,
                                        colnames(df_hub_gradient))

hub_gradient_stats <- df_hub_gradient %>%
  pivot_longer(
    cols = all_of(available_hub_states),
    names_to = "cell_state",
    values_to = "abundance"
  ) %>%
  group_by(cell_state, dist_bin, Sample) %>%
  summarise(median_abund = median(abundance, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(cell_state, dist_bin) %>%
  summarise(
    mean_abund = mean(median_abund, na.rm = TRUE),
    se_abund   = sd(median_abund, na.rm = TRUE) / sqrt(n()),
    n_samples  = n(),
    .groups    = "drop"
  )

# Add community annotation
hub_gradient_stats <- hub_gradient_stats %>%
  mutate(
    community = case_when(
      cell_state %in% C1_members ~ "C1",
      cell_state %in% C2_members ~ "C2",
      TRUE ~ "Other"
    )
  )

# Normalize each cell state to its peak abundance = 1
gradient_norm <- hub_gradient_stats %>%
  group_by(cell_state) %>%
  mutate(
    peak_abund = max(mean_abund, na.rm = TRUE),
    norm_abund = mean_abund / peak_abund,
    norm_se    = se_abund / peak_abund
  ) %>%
  ungroup()

cat("── Normalized gradient per cell state ──────────\n")
gradient_norm %>%
  dplyr::select(community, cell_state, dist_bin,
                mean_abund, norm_abund) %>%
  pivot_wider(names_from = dist_bin,
              values_from = c(mean_abund, norm_abund)) %>%
  arrange(community, cell_state) %>%
  print(n = 20, width = 120)

# Aggregate per community
community_agg <- gradient_norm %>%
  group_by(community, dist_bin) %>%
  summarise(
    agg_norm  = mean(norm_abund, na.rm = TRUE),
    agg_se    = sd(norm_abund, na.rm = TRUE) / sqrt(n()),
    n_members = n(),
    .groups   = "drop"
  )

cat("\n── Aggregate community gradients ────────────────\n")
community_agg %>%
  pivot_wider(names_from = dist_bin,
              values_from = c(agg_norm, agg_se)) %>%
  print(width = 120)

# Plot: aggregate C1 vs C2
community_colors <- c("C1" = "#E64B35", "C2" = "#4DBBD5")

p_agg <- ggplot(community_agg,
                aes(x = dist_bin, y = agg_norm,
                    color = community, group = community)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = agg_norm - agg_se,
                  ymax = agg_norm + agg_se,
                  fill = community),
              alpha = 0.15, color = NA) +
  geom_vline(xintercept = 2.5, linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  annotate("text", x = 1.5, y = 1.05,
           label = "TLS", color = "gray40", size = 3.5,
           fontface = "italic") +
  annotate("text", x = 4.5, y = 1.05,
           label = "non-TLS", color = "gray40", size = 3.5,
           fontface = "italic") +
  scale_color_manual(values = community_colors) +
  scale_fill_manual(values = community_colors) +
  labs(
    title = "Hub Community Spatial Gradient from TLS Centroid",
    subtitle = "Normalized abundance (peak = 1) | Mean across community members",
    x = "Distance to TLS centroid",
    y = "Normalized abundance",
    color = "Community", fill = "Community"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 10),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(plot   = p_agg,
       file   = "../output/SFig.3c_hub_community_gradient.pdf",
       width  = 10,
       height = 7,
       dpi    = 300)

cat("── Hub community gradient saved ────────────────\n")
cat("── Done. All outputs saved to ../output/ ───────\n")

# ============================================================
# END OF SCRIPT
# ============================================================