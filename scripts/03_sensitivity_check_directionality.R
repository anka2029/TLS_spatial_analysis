# ============================================================
# 03_Sensitivity_and_Directionality.R
#
# Purpose:
#   1) Sensitivity analysis: test robustness of hub membership
#      across different edge weight thresholds (85th/90th/95th)
#   2) Hub threshold justification: examine degree distribution
#      drop-off to support top-15 cutoff
#   3) Directionality analysis: decompose hub cell state
#      weighted degree into sender (ligand) vs receiver
#      (receptor) components
#
# Input:
#   combined_obj.RData: Visium samples
#   LR_analysis.RData: pre-computed chemokine network results
#
# Output:
#   03c_hub_directionality.pdf
#   03c_hub_network_role.pdf
#
# Dependencies:
#   dplyr, tidyr, ggplot2, purrr, stringr, patchwork, igraph,
#   scales, conflicted
#
# NOTE:
#   This script requires analyze_LR_network() and
#   get_hub_modules() defined in 03_LR_Network_Analysis.R.
#   Either source that script first or run it in the same
#   session.
# ============================================================

library(conflicted)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(igraph)
library(scales)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::lag)
conflicts_prefer(base::intersect)
conflicts_prefer(base::setdiff)
conflicts_prefer(igraph::V)
conflicts_prefer(igraph::E)
conflicts_prefer(igraph::strength)
conflicts_prefer(igraph::degree)

load("../data/example_data/combined_obj.RData")
load("../data/example_data/LR_analysis.RData")

# ════════════════════════════════════════════════════════════
# SECTION 1: SENSITIVITY ANALYSIS — EDGE WEIGHT THRESHOLD
# ════════════════════════════════════════════════════════════
# Test whether hub membership is robust across different
# percentile thresholds for edge weight filtering.

sensitivity_check <- function(combined_obj, ligand, receptor,
                              thresholds = c(0.85, 0.90, 0.95)) {
  
  results <- list()
  
  for (thr in thresholds) {
    cat("\n=== Testing threshold:", thr, "===\n")
    data <- analyze_LR_network(combined_obj, ligand, receptor,
                               quant_thr = thr)
    modules <- get_hub_modules(data, top_n = 15)
    results[[as.character(thr)]] <- modules$hub_state
  }
  
  cat("\n=== Hub overlap across thresholds ===\n")
  for (i in 1:(length(thresholds) - 1)) {
    a <- results[[i]]
    b <- results[[i + 1]]
    jacc <- length(intersect(a, b)) / length(union(a, b))
    cat(thresholds[i], "vs", thresholds[i + 1],
        "- Jaccard:", round(jacc, 3), "\n")
  }
  
  return(results)
}

# Run sensitivity analysis for each axis
cat("\n══ CXCL13-CXCR5 Sensitivity ═════════════════════\n")
sensitivity_CXCL13 <- sensitivity_check(combined_obj, "CXCL13", "CXCR5")

cat("\n══ CXCL12-CXCR4 Sensitivity ═════════════════════\n")
sensitivity_CXCL12 <- sensitivity_check(combined_obj, "CXCL12", "CXCR4")

cat("\n══ CCL19/21-CCR7 Sensitivity ════════════════════\n")
sensitivity_CCL19 <- sensitivity_check(combined_obj,
                                       c("CCL19", "CCL21"), "CCR7")


# ════════════════════════════════════════════════════════════
# SECTION 2: HUB THRESHOLD JUSTIFICATION
# ════════════════════════════════════════════════════════════
# Examine the degree distribution drop-off across ranks
# to support the top-15 hub cutoff.

check_hub_threshold <- function(axis_data, axis_name, top_range = 20) {
  
  deg_summary <- axis_data %>%
    purrr::map_dfr("node_metrics") %>%
    group_by(cell_state) %>%
    summarise(med_deg = median(degree), .groups = "drop") %>%
    arrange(desc(med_deg)) %>%
    mutate(rank = row_number()) %>%
    filter(rank <= top_range) %>%
    mutate(
      drop     = lag(med_deg) - med_deg,
      drop_pct = drop / lag(med_deg) * 100
    )
  
  cat("\n=== Degree distribution (top", top_range, ") for",
      axis_name, "===\n")
  print(deg_summary[, c("rank", "cell_state", "med_deg", "drop_pct")])
  
  # Find largest drop within top range
  elbow_rank <- deg_summary %>%
    filter(!is.na(drop_pct)) %>%
    arrange(desc(drop_pct)) %>%
    slice(1) %>%
    pull(rank)
  
  cat("Suggested elbow within top", top_range, ": rank", elbow_rank, "\n")
  
  return(deg_summary)
}

cat("\n══ Hub Threshold Check ══════════════════════════\n")
deg_CXCL13 <- check_hub_threshold(LR_analysis$CXCL13_CXCR5, "CXCL13-CXCR5")
deg_CXCL12 <- check_hub_threshold(LR_analysis$CXCL12_CXCR4, "CXCL12-CXCR4")
deg_CCL19  <- check_hub_threshold(LR_analysis$CCL19_21_CCR7, "CCL19/21-CCR7")


# ════════════════════════════════════════════════════════════
# SECTION 3: DEFINE HUB STATES
# ════════════════════════════════════════════════════════════

C1_members <- c("B_S01", "CD4.T_S02", "Dendritic_S01",
                "PMNs_S03", "PCs_S01")
C2_members <- c("Endothelial_S02", "Epithelial_S04",
                "Monocytes.and.Macrophages_S01", "NK_S01", "NK_S03")
all_hubs <- c(C1_members, C2_members)


# ════════════════════════════════════════════════════════════
# SECTION 4: DIRECTIONAL DEGREE COMPUTATION
# ════════════════════════════════════════════════════════════

compute_directional_degree <- function(axis_data, axis_name) {
  
  map_dfr(names(axis_data), function(sid) {
    g <- axis_data[[sid]]$graph
    nodes <- V(g)$name
    
    hubs_present <- intersect(all_hubs, nodes)
    if (length(hubs_present) == 0) return(NULL)
    
    tibble(
      sample     = sid,
      cell_state = hubs_present,
      # Weighted degree (interaction strength)
      wdeg_out   = igraph::strength(g, v = hubs_present, mode = "out",
                                    weights = E(g)$score),
      wdeg_in    = igraph::strength(g, v = hubs_present, mode = "in",
                                    weights = E(g)$score),
      wdeg_total = igraph::strength(g, v = hubs_present, mode = "all",
                                    weights = E(g)$score),
      # Unweighted degree (number of partners)
      deg_out    = igraph::degree(g, v = hubs_present, mode = "out"),
      deg_in     = igraph::degree(g, v = hubs_present, mode = "in"),
      deg_total  = igraph::degree(g, v = hubs_present, mode = "all")
    )
  }) %>%
    mutate(axis = axis_name)
}

dir_degree <- bind_rows(
  compute_directional_degree(LR_analysis$CXCL13_CXCR5, "CXCL13-CXCR5"),
  compute_directional_degree(LR_analysis$CXCL12_CXCR4, "CXCL12-CXCR4"),
  compute_directional_degree(LR_analysis$CCL19_21_CCR7, "CCL19/21-CCR7")
)


# ════════════════════════════════════════════════════════════
# SECTION 5: AGGREGATE — MEDIAN ACROSS SAMPLES
# ════════════════════════════════════════════════════════════

dir_summary <- dir_degree %>%
  group_by(cell_state, axis) %>%
  summarise(
    median_wdeg_out   = median(wdeg_out, na.rm = TRUE),
    median_wdeg_in    = median(wdeg_in, na.rm = TRUE),
    median_wdeg_total = median(wdeg_total, na.rm = TRUE),
    median_deg_out    = median(deg_out, na.rm = TRUE),
    median_deg_in     = median(deg_in, na.rm = TRUE),
    median_deg_total  = median(deg_total, na.rm = TRUE),
    n_samples         = n(),
    .groups           = "drop"
  ) %>%
  mutate(
    # Out-degree proportion (0 = pure receiver, 1 = pure sender)
    out_proportion = median_wdeg_out /
      (median_wdeg_out + median_wdeg_in + 1e-10),
    community = case_when(
      cell_state %in% C1_members ~ "C1",
      cell_state %in% C2_members ~ "C2",
      TRUE ~ "Other"
    )
  )

cat("\n══ Hub Degree Directionality Summary ═════════════\n\n")
dir_summary %>%
  arrange(community, cell_state, axis) %>%
  select(community, cell_state, axis,
         median_wdeg_out, median_wdeg_in, out_proportion,
         median_deg_out, median_deg_in) %>%
  print(n = 50, width = 120)


# ════════════════════════════════════════════════════════════
# SECTION 6: CROSS-AXIS SUMMARY
# ════════════════════════════════════════════════════════════

cross_axis_summary <- dir_summary %>%
  group_by(cell_state, community) %>%
  summarise(
    mean_out_prop   = mean(out_proportion, na.rm = TRUE),
    sd_out_prop     = sd(out_proportion, na.rm = TRUE),
    mean_wdeg_out   = mean(median_wdeg_out, na.rm = TRUE),
    mean_wdeg_in    = mean(median_wdeg_in, na.rm = TRUE),
    mean_wdeg_total = mean(median_wdeg_total, na.rm = TRUE),
    .groups         = "drop"
  ) %>%
  mutate(
    role = case_when(
      mean_out_prop > 0.6 ~ "Predominantly sender",
      mean_out_prop < 0.4 ~ "Predominantly receiver",
      TRUE                ~ "Balanced (sender/receiver)"
    )
  ) %>%
  arrange(community, desc(mean_out_prop))

cat("\n══ Cross-Axis Summary (averaged over 3 axes) ════\n\n")
print(cross_axis_summary, width = 120)


# ════════════════════════════════════════════════════════════
# SECTION 7: VISUALIZATION
# ════════════════════════════════════════════════════════════

# Shorten names for display
dir_summary_plot <- dir_summary %>%
  mutate(cell_state_label = str_replace(cell_state,
                                        "Monocytes\\.and\\.Macrophages",
                                        "Macrophages"))

# --- 7a: Stacked bar — out vs in proportion per hub per axis ---
p_direction <- ggplot(
  dir_summary_plot %>%
    select(cell_state_label, community, axis,
           median_wdeg_out, median_wdeg_in) %>%
    pivot_longer(cols = c(median_wdeg_out, median_wdeg_in),
                 names_to  = "direction",
                 values_to = "wdeg") %>%
    mutate(direction = ifelse(direction == "median_wdeg_out",
                              "Out-degree (sender)",
                              "In-degree (receiver)")),
  aes(x = reorder(cell_state_label, -wdeg),
      y = wdeg, fill = direction)
) +
  geom_bar(stat = "identity", position = "stack",
           color = "white", linewidth = 0.3) +
  facet_grid(axis ~ community, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("Out-degree (sender)"   = "#E64B35",
                               "In-degree (receiver)" = "#4DBBD5"),
                    name = "Degree Direction") +
  labs(
    title    = "Hub Cell State Degree Directionality",
    subtitle = paste0("Weighted degree decomposed into sender (ligand) ",
                      "and receiver (receptor) components"),
    x = NULL,
    y = "Median Weighted Degree"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 10),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

# --- 7b: Out-proportion dot plot (averaged across axes) ---
cross_axis_plot <- cross_axis_summary %>%
  mutate(cell_state_label = str_replace(cell_state,
                                        "Monocytes\\.and\\.Macrophages",
                                        "Macrophages"))

p_balance <- ggplot(cross_axis_plot,
                    aes(x = mean_out_prop,
                        y = reorder(cell_state_label, mean_out_prop),
                        color = community)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray50") +
  geom_point(size = 4) +
  geom_errorbarh(aes(xmin = mean_out_prop - sd_out_prop,
                     xmax = mean_out_prop + sd_out_prop),
                 height = 0.2) +
  scale_color_manual(values = c("C1" = "#E64B35", "C2" = "#4DBBD5")) +
  scale_x_continuous(limits = c(0, 1),
                     labels = percent,
                     breaks = seq(0, 1, 0.25)) +
  annotate("text", x = 0.15, y = 1,
           label = "<- Receiver\n(integrator)",
           color = "gray40", size = 3, fontface = "italic") +
  annotate("text", x = 0.85, y = 1,
           label = "Sender ->\n(broadcaster)",
           color = "gray40", size = 3, fontface = "italic") +
  labs(
    title    = "Hub Cell State Network Role",
    subtitle = "Out-degree proportion (mean +/- SD across 3 chemokine axes)",
    x        = "Out-degree proportion\n(0 = pure receiver, 1 = pure sender)",
    y        = NULL,
    color    = "Community"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 10),
    legend.position  = "bottom"
  )

# --- 7c: Combined output ---
p_combined <- p_direction / p_balance +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = "Directionality Analysis of Hub Cell State Network Interactions",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    )
  )

ggsave(p_combined,
       file   = "../output/03c_hub_directionality_combined.pdf",
       width  = 12,
       height = 14,
       dpi    = 300)

# Also save individual panels
ggsave(p_direction,
       file   = "../output/03c_hub_directionality.pdf",
       width  = 11,
       height = 9,
       dpi    = 300)

ggsave(p_balance,
       file   = "../output/03c_hub_network_role.pdf",
       width  = 8,
       height = 6,
       dpi    = 300)

cat("\n── Done. Outputs saved to ../output/ ────────────\n")

# ============================================================
# END OF SCRIPT
# ============================================================