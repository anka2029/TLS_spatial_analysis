# ============================================================
# 01_TLSvsNonTLS_EcoTypern.R
#
# Purpose:
#   Analyze the association between EcoTyper Carcinoma Ecotypes (CEs)
#   and TLS presence, with a focus on CE9 and CE10.
#
#   This script includes:
#   1) TLS vs NonTLS comparison of CE abundance (boxplots, Wilcoxon test)
#   2) CE9 / CE10 association with TLS at spot level
#   3) Spatial gradient analysis of CE abundance relative to TLS distance
#
# Input:
#   (1) EcoTyper bulk mode output:
#       exampledata/EcoTyper_Output/
#
#   (2) combined_obj (R object loaded in environment):
#       A named list of Visium samples.With output of EcoTyper Visium mode.
#       Each element must contain:
#         - $cellstate_raw or $cellstate_norm(cell state abundance normalized)
#             (ID, pixel_x, pixel_y, Label [TLS / NonTLS])
#         - $ecotype
#             (spot-level CE1–CE10 abundance)
#         - $hs
#             (sparse matrix of gene expression in each barcode)

#
# Output:
#   - Figures comparing CE abundance in TLS vs NonTLS (01_Ecotype_TLS_comparison_with_pvalues.pdf)
#   - CE9 / CE10 TLS association statistics (01_CE9_CE10_MannWhitney_boxplot.pdf)
#   - Spatial CE abundance gradient plots relative to TLS distance (01_CE_abundance_from_TLS.pdf)
#
# ============================================================
## ── Required packages ──────────────────────────────────────────
library(Matrix)      # For sparse matrix rowMeans()
library(dplyr)       # Data manipulation
library(purrr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(rstatix)
library(readr)


# ============================================================================
# PART 1: Overall Ecotype Abundance Comparison (TLS vs NonTLS)
# ============================================================================
## ── Load example data ──────────────────────────────────────────
load("../data/example_data/combined_obj.RData")

## ── 1.1 Load ecotype abundance data ──────────────────────────────────────
eco_file <- "../data/example_data/EcoTyper_Output/Ecotypes/ecotype_abundance.txt"
eco_df <- t(read.csv(eco_file, sep = "\t")) %>% as.data.frame()

## ── 1.2 Assign TLS/NonTLS labels ─────────────────────────────────────────
eco_df <- eco_df %>% 
  mutate(Histology = if_else(grepl("_TLS$", rownames(eco_df)), "TLS", "NonTLS"))

## ── 1.3 Convert to long format ───────────────────────────────────────────
long_df <- eco_df %>% 
  pivot_longer(
    cols      = starts_with("CE"),        # CE1~CE10
    names_to  = "CE",
    values_to = "Abundance"
  ) %>% 
  mutate(
    CE        = factor(CE, levels = paste0("CE", 1:10)),
    Histology = factor(Histology, levels = c("NonTLS", "TLS"))
  )

## ── 1.4 Show p-values for CE9 and CE10 only ─────────────────
# Calculate p-values for CE9 and CE10
p_values <- long_df %>%
  filter(CE %in% c("CE9", "CE10")) %>%
  group_by(CE) %>%
  summarise(
    p_value = wilcox.test(Abundance ~ Histology)$p.value,
    y_pos   = max(Abundance) * 1.1
  )

# Plot with asterisks for all CEs and exact p-values for CE9/CE10
p_ecotype_comparison <- ggplot(long_df, aes(Histology, Abundance, fill = Histology)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.15, size = 0.8, alpha = 0.6) +
  facet_wrap(~ CE, ncol = 3) +
  # Add p-values for CE9 and CE10 only (in scientific notation)
  geom_text(
    data = p_values,
    aes(x = 1.5, y = y_pos, 
        label = paste("p =", formatC(p_value, format = "e", digits = 2))),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(TLS    = "#1E6091",   
               NonTLS = "#A8DADC")   
  ) +
  # Show asterisks for all CEs
  geom_signif(
    comparisons      = list(c("TLS", "NonTLS")),
    test             = "wilcox.test",
    map_signif_level = TRUE,
    margin_top       = 0.01,   
    step_increase    = 0.01 
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face = "bold")
  ) +
  labs(
    title = "Ecotype Abundance in TLS vs NonTLS Spots",
    x     = NULL,
    y     = "EcoType abundance"
  )

# Save the plot
ggsave(
  filename = "../output/01_Ecotype_TLS_comparison_with_pvalues.pdf",
  plot     = p_ecotype_comparison,
  width    = 10,
  height   = 12,
  dpi      = 300
)

# ============================================================================
# PART 2: CE9/CE10 Association with TLS at Spot Level
# ============================================================================

## ── 2.1 Extract CE9/CE10 abundance and TLS labels per spot ───────────────
res_list <- list()
k <- 1

for (samp_name in names(combined_obj)) {
  
  samp_obj <- combined_obj[[samp_name]]
  
  # Check if cellstate_raw exists
  if (is.null(samp_obj$cellstate_raw)) {
    message("[SKIP] ", samp_name, ": cellstate_raw not found")
    next
  }
  cs <- samp_obj$cellstate_raw
  
  # Check if Label column exists
  if (!"Label" %in% colnames(cs)) {
    message("[SKIP] ", samp_name, ": Label column not found in cellstate_raw")
    next
  }
  
  # Check if ecotype data exists
  if (is.null(samp_obj$ecotype)) {
    message("[SKIP] ", samp_name, ": ecotype information not found")
    next
  }
  eco <- samp_obj$ecotype
  
  # Find common barcodes
  common_bc <- intersect(rownames(eco), cs$ID)
  if (length(common_bc) == 0) {
    message("[SKIP] ", samp_name, ": No common barcodes found")
    next
  }
  
  # Extract CE9/CE10 and assign TLS flag
  df <- eco[common_bc, c("CE9", "CE10")] %>%
    as.data.frame() %>%
    mutate(
      Barcode  = common_bc,
      TLS_flag = cs[cs$ID %in% common_bc, ]$Label,
      Sample   = samp_name
    )
  rownames(df) <- NULL
  
  # Add to list
  res_list[[k]] <- df
  k <- k + 1
}

## ── 2.2 Combine all samples ──────────────────────────────────────────────
ce_tls_list <- bind_rows(res_list)

## ── 2.3 Mann-Whitney U test ────────────────────────
# Prepare data in long format
ce_long <- ce_tls_list %>%
  filter(!is.na(CE9) & !is.na(CE10)) %>%
  pivot_longer(
    cols      = c(CE9, CE10),
    names_to  = "CE_type",
    values_to = "Abundance"
  )

# Perform Wilcoxon rank-sum tests
ce9_test <- wilcox.test(
  CE9 ~ TLS_flag, 
  data  = ce_tls_list %>% filter(!is.na(CE9)),
  exact = FALSE
)

ce10_test <- wilcox.test(
  CE10 ~ TLS_flag,
  data  = ce_tls_list %>% filter(!is.na(CE10)),
  exact = FALSE
)

# Prepare p-value labels
p_labels <- data.frame(
  CE_type = c("CE9", "CE10"),
  p_label = c("p < 2.2e-16", "p < 2.2e-16"),
  y_pos   = c(1.15, 1.15)
)

# Create boxplot with Mann-Whitney test results
p <- ggplot(ce_long, aes(x = TLS_flag, y = Abundance, fill = TLS_flag)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  facet_wrap(~ CE_type, scales = "free_y") +
  geom_text(
    data        = p_labels,
    aes(x = 1.5, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size        = 4
  ) +
  scale_fill_manual(values = c("NonTLS" = "#E8E8E8", "TLS" = "#1E6091")) +
  scale_y_continuous(limits = c(0, 1.25)) +
  labs(
    x     = NULL,
    y     = "CE Abundance",
    title = "CE9 and CE10 Abundance in TLS vs Non-TLS Regions"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 12, face = "bold"),
    axis.text       = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, size = 14)
  )



# Save figure
ggsave(plot = p,filename = "../output/01_CE9_CE10_MannWhitney_boxplot.pdf", width = 8, height = 5)


# ============================================================================
# PART 3: Distance-Based Analysis (Unified Physical Distance)
# ============================================================================

## ── 3.1 Calculate distances (same as Part 3.1) ───────────────────────────
df_list <- lapply(names(combined_obj), function(samp) {
  samp_obj <- combined_obj[[samp]]
  if (is.null(samp_obj$cellstate_norm) || is.null(samp_obj$ecotype)) return(NULL)
  
  cs  <- samp_obj$cellstate_norm %>% select(ID, pixel_x, pixel_y, Label)
  eco <- samp_obj$ecotype %>% 
    as.data.frame() %>% 
    select(ID, paste0("CE", 1:10))
  
  df <- inner_join(cs, eco, by = "ID") %>% mutate(Sample = samp)
  
  # Get TLS coordinates
  tls_coords <- df %>% filter(Label == "TLS") %>% select(pixel_x, pixel_y)
  
  # Handle samples without TLS
  if (nrow(tls_coords) == 0) {
    df$dist_to_TLS <- NA
    return(df)
  }
  
  # Calculate distance
  df$dist_to_TLS <- apply(
    df[, c("pixel_x", "pixel_y")], 1,
    function(pt) {
      d <- sqrt((tls_coords$pixel_x - pt[1])^2 + 
                  (tls_coords$pixel_y - pt[2])^2)
      min(d, na.rm = TRUE)
    }
  )
  df
})
df_all <- bind_rows(df_list)

## ── 3.2 Filter valid samples ─────────────────────────────────────────────
valid_samples <- df_all %>%
  group_by(Sample) %>%
  summarise(
    n_tls   = sum(Label == "TLS", na.rm = TRUE),
    n_total = n(),
    .groups = "drop"
  ) %>%
  filter(n_tls > 0) %>%
  pull(Sample)

cat("Valid samples with TLS:", length(valid_samples), "\n")

## ── 3.3 Check distance distribution ──────────────────────────────────────
dist_summary <- df_all %>%
  filter(Sample %in% valid_samples, dist_to_TLS > 0) %>%
  summarise(
    min    = min(dist_to_TLS, na.rm = TRUE),
    q10    = quantile(dist_to_TLS, 0.10, na.rm = TRUE),
    q25    = quantile(dist_to_TLS, 0.25, na.rm = TRUE),
    median = median(dist_to_TLS, na.rm = TRUE),
    q75    = quantile(dist_to_TLS, 0.75, na.rm = TRUE),
    q90    = quantile(dist_to_TLS, 0.90, na.rm = TRUE),
    max    = max(dist_to_TLS, na.rm = TRUE),
    mean   = mean(dist_to_TLS, na.rm = TRUE),
    sd     = sd(dist_to_TLS, na.rm = TRUE)
  )

print("Distance distribution (in pixels):")
print(dist_summary)

## ── 3.4 Bin by unified physical distance ─────────────────────────────────
# Option 1: Quartile-based thresholds (rounded values)
df_binned <- df_all %>%
  filter(Sample %in% valid_samples) %>%
  mutate(
    dist_category = case_when(
      dist_to_TLS == 0            ~ "TLS (0)",
      dist_to_TLS <= 750          ~ "Near (<750)",         # ~Q25
      dist_to_TLS <= 1750         ~ "Medium (750-1750)",   # ~Median
      dist_to_TLS <= 3000         ~ "Far (1750-3000)",     # ~Q75
      dist_to_TLS > 3000          ~ "Very Far (>3000)",
      TRUE                        ~ NA_character_
    ),
    dist_category = factor(dist_category, 
                           levels = c("TLS (0)", "Near (<750)", "Medium (750-1750)", 
                                      "Far (1750-3000)", "Very Far (>3000)"))
  ) %>%
  filter(!is.na(dist_category))

# Check distribution
cat("\nDistribution of spots by distance category:\n")
category_dist <- df_binned %>%
  group_by(dist_category) %>%
  summarise(
    n_spots   = n(),
    pct       = n() / nrow(.) * 100,
    n_samples = n_distinct(Sample),
    .groups   = "drop"
  )
print(category_dist)

cat("\nDistance category distribution per sample:\n")
sample_category_dist <- df_binned %>%
  group_by(Sample, dist_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = dist_category, values_from = n, values_fill = 0)
print(sample_category_dist)

## ── 3.5 Calculate statistics per bin ─────────────────────────────────────
df_bin_summary <- df_binned %>%
  select(Sample, dist_category, paste0("CE", 1:10)) %>%
  pivot_longer(cols = paste0("CE", 1:10), names_to = "CE", values_to = "Abundance") %>%
  group_by(Sample, dist_category, CE) %>%
  summarise(
    median_abundance = median(Abundance, na.rm = TRUE),
    mean_abundance   = mean(Abundance, na.rm = TRUE),
    n                = n(),
    .groups          = "drop"
  )

# Average across samples
df_bin_avg <- df_bin_summary %>%
  group_by(dist_category, CE) %>%
  summarise(
    mean      = mean(median_abundance, na.rm = TRUE),
    se        = sd(median_abundance, na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups   = "drop"
  )

## ── 3.6 Visualize gradient ───────────────────────────────────────────────
df_bin_avg_styled <- df_bin_avg %>%
  mutate(
    CE           = factor(CE, levels = paste0("CE", 1:10)),
    is_highlight = CE %in% c("CE9", "CE10")
  )

# Prepare labels
df_labels_connected <- df_bin_avg_styled %>%
  group_by(CE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(mean) %>%
  mutate(
    point_y = mean,
    label_y = seq(from = min(mean) - 0.02, 
                  to = max(mean) + 0.08, 
                  length.out = 10),
    point_x = 1,
    label_x = 0.4
  )

df_value_labels <- df_bin_avg_styled %>%
  filter(is_highlight) %>%
  mutate(label_text = sprintf("%.2f", mean))

# Create plot
p1_with_value_labels <- ggplot(df_bin_avg_styled, 
                               aes(x = dist_category, y = mean, group = CE)) +
  geom_line(data = df_bin_avg_styled %>% filter(!is_highlight),
            aes(color = CE), linewidth = 0.8, alpha = 0.5) +
  geom_point(data = df_bin_avg_styled %>% filter(!is_highlight),
             aes(color = CE), size = 2, alpha = 0.5) +
  geom_errorbar(data = df_bin_avg_styled %>% filter(!is_highlight),
                aes(ymin = mean - se, ymax = mean + se, color = CE), 
                width = 0.2, alpha = 0.3) +
  geom_line(data = df_bin_avg_styled %>% filter(is_highlight),
            aes(color = CE), linewidth = 1.8) +
  geom_point(data = df_bin_avg_styled %>% filter(is_highlight),
             aes(color = CE), size = 4) +
  geom_errorbar(data = df_bin_avg_styled %>% filter(is_highlight),
                aes(ymin = mean - se, ymax = mean + se, color = CE), 
                width = 0.3, alpha = 0.8, linewidth = 1) +
  geom_label(
    data = df_value_labels,
    aes(label = label_text, fill = CE),
    color         = "white",
    vjust         = -0.8,
    size          = 3,
    fontface      = "bold",
    label.padding = unit(0.15, "lines"),
    label.size    = 0.2,
    show.legend   = FALSE
  ) +
  geom_segment(
    data = df_labels_connected,
    aes(x = point_x, y = point_y, xend = label_x, yend = label_y, color = CE),
    linewidth = 0.4,
    alpha     = 0.6,
    linetype  = "solid"
  ) +
  geom_text(
    data = df_labels_connected,
    aes(x = label_x, y = label_y, label = CE, color = CE),
    hjust    = 1,
    size     = ifelse(df_labels_connected$is_highlight, 5, 3.5),
    fontface = ifelse(df_labels_connected$is_highlight, "bold", "plain")
  ) +
  scale_color_manual(values = ce_colors) +
  scale_fill_manual(values = ce_colors) +
  scale_x_discrete(expand = expansion(mult = c(0.25, 0.05))) +
  coord_cartesian(clip = "off") +
  labs(
    title    = "CE Abundance Gradient from TLS",
    subtitle = sprintf("Physical distance categories (n=%d samples with TLS) | CE9 & CE10 highlighted", 
                       length(valid_samples)),
    x        = "Distance from TLS (pixels)",
    y        = "Mean CE Abundance ± SE"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y   = element_text(size = 12),
    axis.title    = element_text(size = 13),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.margin   = margin(10, 10, 10, 60)
  )

ggsave(plot = p1_with_value_labels, file = "../output/01_CE_abundance_from_TLS.pdf")


# ============================================================================
# END OF SCRIPT
# ============================================================================