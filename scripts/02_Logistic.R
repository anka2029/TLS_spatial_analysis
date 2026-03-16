# ============================================================
# 02_LogisticModel_cellstate.R
#
# Purpose:
#   Identify cell states associated with TLS presence:
#   1) Random Forest (ranger) 
#   2) H-statistic (hstats) for interaction detection
#   3) Multivariable mixed-effect logistic regression
#      — p-value cutoff for interaction significance
#   4) Visualization
#
# Input:
#   combined_obj: named list of Visium samples
#   multi_regression: pre-computed regression result
#
# Dependencies:
#   ranger, hstats, lme4, broom.mixed, car,
#   dplyr, tidyr, ggplot2, ggrepel, purrr, stringr, patchwork
# ============================================================

library(conflicted)
library(ranger)
library(hstats)
library(lme4)
library(broom.mixed)
library(car)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(purrr)
library(stringr)
library(patchwork)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::setdiff)
conflicts_prefer(base::intersect)
set.seed(42)

load("../data/example_data/combined_obj.RData")
load("../data/example_data/multi_regression.RData")
# ════════════════════════════════════════════════════════════
# SECTION 1: DATA PREPARATION
# ════════════════════════════════════════════════════════════

df_all <- purrr::map_dfr(names(combined_obj), function(samp) {
  cs <- combined_obj[[samp]]$cellstate_norm
  if (is.null(cs)) return(NULL)
  cs %>%
    dplyr::select(ID, array_x, array_y, Label,
                  matches("_S[0-9]+$")) %>%
    mutate(Sample = samp,
           y      = as.integer(Label == "TLS"))
})

# Select samples with TLS
tls_count <- tapply(df_all$y, df_all$Sample, sum)
samples_with_TLS <- names(tls_count[tls_count > 0])
df_all <- df_all[df_all$Sample %in% samples_with_TLS, ]

all_states <- df_all %>%
  dplyr::select(matches("_S[0-9]+$")) %>%
  colnames() %>% setdiff("Fibroblasts_S07") # remove doublet

cat("── Data summary ──────────────────────────────────\n")
cat("Total spots   :", nrow(df_all), "\n")
cat("TLS spots     :", sum(df_all$y), "\n")
cat("Non-TLS spots :", sum(df_all$y == 0), "\n")
cat("TLS prevalence:", round(mean(df_all$y) * 100, 1), "%\n")
cat("Samples       :", n_distinct(df_all$Sample), "\n")
cat("Cell states   :", length(all_states), "\n\n")

# ════════════════════════════════════════════════════════════
# SECTION 2: RANDOM FOREST — VARIABLE IMPORTANCE
# ════════════════════════════════════════════════════════════

# --- 2.1 Prepare data for RF ---
df_rf <- df_all %>%
  dplyr::select(y, all_of(all_states)) %>%
  mutate(y = factor(y, levels = c(0, 1), labels = c("nonTLS", "TLS"))) %>%
  na.omit()

rf_fit <- ranger(
  y ~ .,
  data        = df_rf,
  num.trees   = 1000,
  importance  = "permutation",
  probability = TRUE
)

print(rf_fit)

# --- 2.4 Variable Importance ---
vi <- sort(rf_fit$variable.importance, decreasing = TRUE)

# Select top 30
n_top_main <- 30
top_vars <- names(vi)[1:n_top_main]

cat("\nSelected main effects (top", n_top_main, "):\n")
cat(paste(seq_along(top_vars), top_vars, sep = ". ", collapse = "\n"), "\n\n")

# Variable importance plot
vi_df <- data.frame(
  variable   = names(vi),
  importance = vi,
  rank       = seq_along(vi),
  selected   = names(vi) %in% top_vars
) %>%
  head(30)

p_vi <- ggplot(vi_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", aes(fill = selected)) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "gray60"),
                    guide  = "none") +
  coord_flip() +
  labs(title = "Random Forest Variable Importance (Permutation)",
       subtitle = "Top 30 selected",
       x = NULL, y = "Importance") +
  theme_minimal(base_size = 12)



# ════════════════════════════════════════════════════════════
# SECTION 3: INTERACTION DETECTION (H-STATISTIC)
# ════════════════════════════════════════════════════════════

cat("Computing H-statistics for top", n_top_main, "variables ...\n")

# hstats internally calls predict(), which needs ALL predictor columns
X_sub <- df_rf %>% dplyr::select(all_of(all_states)) %>% as.data.frame()

H <- hstats(rf_fit, v = top_vars, X = X_sub)

# Overall interaction strength
h2_overall_vals <- h2_overall(H, squared = FALSE, zero = FALSE)

# Pairwise interaction strength
h2_pair <- h2_pairwise(H, squared = FALSE, zero = FALSE)

cat("\n── Pairwise Interactions (all) ───────────────────\n")
print(h2_pair)

# Extract pairwise H values as data frame
h2_df <- as.data.frame(h2_pair$M)

if (ncol(h2_df) == 2) {
  h2_df <- data.frame(
    interaction_pair = rownames(h2_df),
    H_stat           = h2_df[, 1]
  ) %>%
    arrange(desc(H_stat))
}

cat("\n── All pairwise H-statistics ─────────────────────\n")
cat("Total pairs:", nrow(h2_df), "\n")
cat("H-stat range:", round(min(h2_df$H_stat), 4), "-",
    round(max(h2_df$H_stat), 4), "\n\n")

# ════════════════════════════════════════════════════════════
# SECTION 4: MULTIVARIABLE MIXED-EFFECT LOGISTIC REGRESSION
# ════════════════════════════════════════════════════════════

# --- 4.1 Prepare interaction term names ---
interaction_terms <- h2_df$interaction_pair

# Collect all unique variables involved in interactions
vars_in_interactions <- unique(unlist(strsplit(interaction_terms, ":")))

# Final set of main effects: union of top_vars and any vars in interactions
selected_main <- unique(c(top_vars, vars_in_interactions))
selected_main <- intersect(selected_main, all_states)

cat("\n── Final Model ──────────────────────────────────\n")
cat("Main effects     :", length(selected_main), "\n")
cat("Interaction terms:", length(interaction_terms), "\n")

# --- 4.2 Prepare data (full dataset, standardized) ---
df_multi <- df_all %>%
  dplyr::select(y, Sample, all_of(selected_main)) %>%
  na.omit() %>%
  mutate(across(all_of(selected_main), ~ as.numeric(scale(.x))))

cat("Regression data  :", nrow(df_multi), "spots\n\n")

# --- 4.3 Build formula ---
formula_str <- paste(
  "y ~",
  paste(selected_main, collapse = " + "),
  "+",
  paste(interaction_terms, collapse = " + "),
  "+ (1 | Sample)"
)

cat("Formula:\n", formula_str, "\n\n")
formula_multi <- as.formula(formula_str)

# --- 4.4 Fit model ---
cat("Fitting multivariable mixed logistic regression ...\n")
fit_multi <- glmer(
  formula_multi,
  family  = binomial,
  data    = df_multi,
  control = glmerControl(optimizer = "bobyqa",
                         optCtrl   = list(maxfun = 2e5))
)

summary(fit_multi)

# --- 4.5 VIF check (main effects only) ---
# Per reviewer: VIF check is not necessary for interaction terms
cat("\n── VIF Check (main effects only) ─────────────────\n")

formula_main_only <- paste("y ~", paste(selected_main, collapse = " + "))
fit_main_only <- glm(
  as.formula(formula_main_only),
  data   = df_multi %>% dplyr::select(-Sample),
  family = binomial
)
vif_values <- car::vif(fit_main_only)

cat("VIF > 5 (problematic):\n")
high_vif <- vif_values[vif_values > 5]
if (length(high_vif) > 0) {
  print(high_vif)
} else {
  cat("  None — all VIF values acceptable\n")
}
cat("VIF range:", round(min(vif_values), 2), "-",
    round(max(vif_values), 2),
    "(median:", round(median(vif_values), 2), ")\n\n")

# --- 4.6 Extract results ---
multi_results <- broom.mixed::tidy(
  fit_multi,
  conf.int     = TRUE,
  exponentiate = TRUE
) %>%
  dplyr::filter(term != "(Intercept)") %>%
  mutate(
    FDR            = p.adjust(p.value, method = "fdr"),
    is_interaction = grepl(":", term)
  ) %>%
  arrange(p.value)

cat("── Results (all terms) ──────────────────────────\n")
print(multi_results, n = 50)

# Separate significant interactions
sig_interactions <- multi_results %>%
  filter(is_interaction, p.value < 0.05)

sig_interactions_fdr <- multi_results %>% filter(is_interaction, FDR < 0.05)

cat("\n── Summary ──────────────────────────────────────\n")
cat("Total interaction terms tested:", sum(multi_results$is_interaction), "\n")
cat("Significant at p < 0.05      :", nrow(sig_interactions), "\n")
cat("Significant at FDR < 0.05    :", nrow(sig_interactions_fdr), "\n")

save(multi_results, rf_fit, vi, H, h2_pair, top_interactions,
     h2_df, h_cutoff,
     file = "../data/example_data/multi_regression.RData")


# ════════════════════════════════════════════════════════════
# SECTION 5: VISUALIZATION
# ════════════════════════════════════════════════════════════

plot_df <- multi_results %>%
  mutate(
    log2OR        = log2(estimate),
    neg_log10_FDR = -log10(FDR),
    sig           = FDR < 0.05,
    term_clean    = str_replace_all(term, "[._]", " ")
  )

# --- 5.1 Volcano plot (main effects) ---
plot_main <- plot_df %>% dplyr::filter(!is_interaction)

p_volcano <- ggplot(plot_main, aes(x = log2OR, y = neg_log10_FDR)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", linewidth = 0.6) +
  geom_point(aes(size = sig, alpha = sig, color = sig)) +
  ggrepel::geom_text_repel(
    data          = plot_main %>% filter(sig),
    aes(label = term_clean),
    size          = 3.5, color = "black",
    box.padding   = 0.5, point.padding = 0.3,
    max.overlaps  = 20, segment.color = "gray70"
  ) +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "gray60"), guide = "none") +
  scale_size_manual(values  = c("TRUE" = 3.5, "FALSE" = 2), guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.5), guide = "none") +
  annotate("text",
           x = max(plot_main$log2OR, na.rm = TRUE) * 0.8,
           y = max(plot_main$neg_log10_FDR, na.rm = TRUE) * 0.95,
           label = "TLS-enriched →", color = "gray40", size = 3.5) +
  annotate("text",
           x = min(plot_main$log2OR, na.rm = TRUE) * 0.8,
           y = max(plot_main$neg_log10_FDR, na.rm = TRUE) * 0.95,
           label = "← TLS-depleted", color = "gray40", size = 3.5) +
  labs(
    title    = "Cell State Associations with TLS (Main Effects)",
    subtitle = "Multivariable mixed-effect logistic regression | FDR < 0.05",
    x = "log\u2082(Odds Ratio)", y = "-log\u2081\u2080(FDR)",
    caption  = "(1|Sample) random intercept | Predictors standardized"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 15),
    plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 10),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    plot.caption     = element_text(color = "gray50", size = 9)
  )

ggsave(plot = p_volcano, file = "../output/02_volcano_main_effects.pdf",
       width = 9, height = 7, dpi = 300)

# --- 5.2 Forest plot (interaction effects) ---
plot_int <- plot_df %>%
  filter(is_interaction) %>%
  arrange(log2OR) %>%
  mutate(term_clean = factor(term_clean, levels = term_clean))

if (nrow(plot_int) > 0) {
  p_forest <- ggplot(plot_int, aes(x = log2OR, y = term_clean)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = log2(conf.low), xmax = log2(conf.high)),
                   height = 0.25, color = "gray40") +
    geom_point(aes(color = sig), size = 3) +
    scale_color_manual(
      values = c("TRUE" = "#E64B35", "FALSE" = "gray60"),
      labels = c("TRUE" = "FDR < 0.05", "FALSE" = "NS"),
      name   = ""
    ) +
    labs(
      title    = "Interaction Effects on TLS Presence",
      x = "log\u2082(Odds Ratio)", y = NULL,
      caption  = "Error bars: 95% CI | Mixed-effect logistic regression"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle    = element_text(hjust = 0.5, color = "gray40", size = 10),
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(color = "gray50", size = 9)
    )
  
  ggsave(plot = p_forest, file = "../output/02_forest_interactions.pdf",
         width = 10, height = max(4, nrow(plot_int) * 0.4 + 2), dpi = 300)
}

cat("\n── Done ──────────────────────────────────────────\n")
cat("Outputs saved to ../output/\n")