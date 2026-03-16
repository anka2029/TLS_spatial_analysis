# ============================================================
# 02_CellState_Functional_Characterization.R
#
# Figure 3: Functional differences among cell states with
#           divergent TLS associations within the same cell type
#
#
# Cell types analyzed:
#   Macrophage: S01 (TLS+), S07 (TLS−), S08 (TLS−)
#   CD4.T:     S01–S04, S06, S07
#
# ============================================================
library(conflicted)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ComplexHeatmap)
library(circlize)
library(purrr)
library(stringr)
library(Matrix)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::setdiff)
conflicts_prefer(base::intersect)

set.seed(42)

load("../data/example_data/combined_obj.RData")
load("../data/example_data/multi_regression.RData")
all_marker_gene <- readRDS("../data/example_data/gene_info/all_marker_gene.rds")

# ════════════════════════════════════════════════════════════
# SECTION 0: SETUP
# ════════════════════════════════════════════════════════════

analysis_targets <- list(
  Macrophage = list(
    prefix = "Monocytes.and.Macrophages",
    states = c("S01", "S07", "S08")
  ),
  CD4T = list(
    prefix = "CD4.T",
    states = c("S01", "S02", "S03", "S04", "S06", "S07")
  )
)

# Check marker gene naming convention
# NOTE: all_marker_gene$cell_state may use spaces or dots.
#       Adjust this line if needed after checking:
#       unique(all_marker_gene$cell_state)
cat("── Marker gene naming examples ──────────────────\n")
print(head(unique(all_marker_gene$cell_state), 10))
cat("\n")

# Helper: convert between dot and space naming if needed
to_marker_name <- function(prefix, state) {
  # Try dot version first; if no match, try space version
  dot_name   <- paste0(prefix, "_", state)
  space_name <- gsub("\\.", " ", dot_name)
  
  if (dot_name %in% all_marker_gene$cell_state) return(dot_name)
  if (space_name %in% all_marker_gene$cell_state) return(space_name)
  
  # Try other variants
  under_name <- gsub("\\.", "_", dot_name)
  if (under_name %in% all_marker_gene$cell_state) return(under_name)
  
  warning("No match for: ", dot_name)
  return(dot_name)
}

# Classify TLS association from regression results
get_tls_group <- function(prefix, states) {
  state_terms <- paste0(prefix, "_", states)
  multi_results %>%
    filter(term %in% state_terms, !is_interaction) %>%
    mutate(
      tls_group = case_when(
        FDR < 0.05 & estimate > 1 ~ "TLS-Enriched",
        FDR < 0.05 & estimate < 1 ~ "TLS-Depleted",
        TRUE                       ~ "NS"
      ),
      state_id = str_extract(term, "S[0-9]+$")
    ) %>%
    arrange(desc(estimate))
}

state_class <- map_dfr(names(analysis_targets), function(ct) {
  info <- analysis_targets[[ct]]
  get_tls_group(info$prefix, info$states) %>% mutate(cell_type = ct)
})

cat("── TLS classification ───────────────────────────\n")
print(state_class %>% select(cell_type, term, estimate, FDR, tls_group))
cat("\n")


# ════════════════════════════════════════════════════════════
# PART 1: DEG HEATMAP
# ════════════════════════════════════════════════════════════

cat("═══ PART 1: DEG Heatmaps ═══════════════════════\n\n")

make_deg_heatmap <- function(ct_name, info, n_genes = 15) {
  
  # --- Get marker genes per state ---
  markers <- map_dfr(info$states, function(s) {
    state_name <- to_marker_name(info$prefix, s)
    all_marker_gene %>%
      filter(cell_state == state_name) %>%
      slice_head(n = n_genes) %>%
      mutate(state = s)
  })
  
  if (nrow(markers) == 0) {
    cat("  WARNING: no markers found for", ct_name, "\n")
    return(NULL)
  }
  
  all_genes <- unique(markers$Gene)
  cat("  ", ct_name, ":", length(all_genes), "unique genes across",
      length(info$states), "states\n")
  
  # --- Compute mean expression per state across samples ---
  # For each sample: assign spots to dominant state, compute mean expr
  state_cols <- paste0(info$prefix, "_", info$states)
  
  expr_list <- map(names(combined_obj), function(samp) {
    obj <- combined_obj[[samp]]
    hs  <- obj$hs               # dgCMatrix: genes × spots
    cs  <- obj$cellstate_norm
    
    # Check gene and state column availability
    genes_present <- intersect(all_genes, rownames(hs))
    cols_present  <- intersect(state_cols, colnames(cs))
    if (length(genes_present) == 0 || length(cols_present) == 0) return(NULL)
    
    # Match spots between hs and cellstate_norm
    common_spots <- intersect(colnames(hs), cs$ID)
    if (length(common_spots) == 0) return(NULL)
    
    # Get dominant state per spot
    state_mat <- cs[match(common_spots, cs$ID), cols_present, drop = FALSE]
    state_mat <- as.matrix(state_mat)
    storage.mode(state_mat) <- "numeric"
    
    dominant_idx <- apply(state_mat, 1, function(x) {
      if (all(is.na(x))) return(NA)
      which.max(x)
    })
    
    dominant <- cols_present[dominant_idx]
    dominant_id <- str_extract(dominant, "S[0-9]+$")
    # Extract expression
    expr_sub <- as.matrix(hs[genes_present, common_spots, drop = FALSE])
    
    # Mean per state
    map_dfr(info$states, function(s) {
      idx <- which(dominant_id == s)
      if (length(idx) == 0) return(NULL)
      tibble(
        state = s,
        gene  = genes_present,
        mean_expr = rowMeans(expr_sub[, idx, drop = FALSE])
      )
    })
  })
  
  expr_df <- bind_rows(expr_list)
  if (nrow(expr_df) == 0) return(NULL)
  
  # Average across samples
  expr_avg <- expr_df %>%
    group_by(state, gene) %>%
    summarise(expr = mean(mean_expr, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = state, values_from = expr) %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  
  # Z-score per gene (across states)
  expr_z <- t(scale(t(expr_avg)))
  expr_z[is.nan(expr_z)] <- 0
  
  # --- Draw heatmap ---
  # Column annotation: TLS group
  ct_class <- state_class %>% filter(cell_type == ct_name)
  col_groups <- ct_class$tls_group[match(colnames(expr_z), ct_class$state_id)]
  
  col_ha <- HeatmapAnnotation(
    TLS = col_groups,
    col = list(TLS = c("TLS-Enriched" = "#E64B35",
                       "TLS-Depleted" = "#4DBBD5",
                       "NS"           = "gray70")),
    show_annotation_name = TRUE
  )
  
  # Row split: which state the gene is a marker for
  gene_state <- markers %>%
    distinct(Gene, state) %>%
    group_by(Gene) %>%
    summarise(origin = first(state), .groups = "drop")
  row_split <- gene_state$origin[match(rownames(expr_z), gene_state$Gene)]
  
  # Column labels
  col_labels <- paste0(info$prefix, "_", colnames(expr_z))
  col_labels <- gsub("Monocytes.and.Macrophages", "M/M", col_labels)
  col_labels <- gsub("CD4.T", "CD4.T", col_labels)
  
  ht <- Heatmap(
    expr_z,
    name = "Z-score",
    col  = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
    
    top_annotation  = col_ha,
    cluster_columns = FALSE,
    cluster_rows    = TRUE,
    show_row_dend   = FALSE,
    
    row_split    = factor(row_split, levels = info$states),
    row_gap      = unit(2, "mm"),
    row_title_gp = gpar(fontsize = 9, fontface = "bold"),
    
    column_labels   = col_labels,
    row_names_gp    = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    
    column_title    = paste0(ct_name, " — Marker Gene Expression"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    
    border = TRUE,
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  pdf(paste0("../output/03_fig3a_DEG_heatmap_", ct_name, ".pdf"),
      width  = 4 + length(info$states) * 1.2,
      height = max(6, nrow(expr_z) * 0.2 + 3))
  draw(ht, heatmap_legend_side = "bottom",
       padding = unit(c(10, 10, 10, 10), "mm"))
  dev.off()
  
  cat("  Saved:", ct_name, "heatmap\n")
  invisible(ht)
}

# Run
walk(names(analysis_targets), function(ct) {
  make_deg_heatmap(ct, analysis_targets[[ct]])
})

cat("\n")


# ════════════════════════════════════════════════════════════
# PART 2: GO ENRICHMENT DOT PLOTS
# ════════════════════════════════════════════════════════════

cat("═══ PART 2: GO Enrichment ══════════════════════\n\n")

all_go_tables <- list()

for (ct in names(analysis_targets)) {
  info <- analysis_targets[[ct]]
  cat("── GO:", ct, "────────────────────────────────\n")
  
  # Build named gene list (ENTREZID) per state
  gene_lists <- setNames(
    lapply(info$states, function(s) {
      state_name <- to_marker_name(info$prefix, s)
      symbols <- all_marker_gene %>%
        filter(cell_state == state_name) %>%
        pull(Gene)
      
      entrez <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = symbols, keytype = "SYMBOL", columns = "ENTREZID"
      ) %>%
        filter(!is.na(ENTREZID)) %>%
        pull(ENTREZID) %>%
        unique()
      
      cat("  ", paste0(info$prefix, "_", s), ":",
          length(symbols), "→", length(entrez), "mapped\n")
      entrez
    }),
    paste0(info$prefix, "_", info$states)
  )
  
  gene_lists <- gene_lists[sapply(gene_lists, length) > 0]
  
  if (length(gene_lists) < 2) {
    cat("  Skipping (< 2 valid states)\n")
    next
  }
  
  # compareCluster
  cc_go <- compareCluster(
    geneClusters  = gene_lists,
    fun           = "enrichGO",
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  go_df <- as.data.frame(cc_go)
  go_df$cell_type <- ct
  all_go_tables[[ct]] <- go_df
  cat("  Enriched terms:", nrow(go_df), "\n")
  
  # Dot plot
  n_show <- 5
  p_go <- dotplot(cc_go, showCategory = n_show) +
    labs(
      title   = paste0("GO BP: ", ct, " States"),
      caption = paste0("Top ", n_show,
                       " per state | BH p.adj < 0.05")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y  = element_text(size = 8),
      plot.caption = element_text(color = "gray50", size = 8)
    )
  
  ggsave(p_go,
         file   = paste0("../output/03_fig3b_GO_dotplot_", ct, ".pdf"),
         width  = 12,
         height = max(6, n_show * length(info$states) * 0.4 + 3),
         dpi    = 300)
  cat("  Saved dotplot\n\n")
}


# ════════════════════════════════════════════════════════════
# PART 3: SPATIAL CO-OCCURRENCE (TLS-ENRICHED STATES)
# ════════════════════════════════════════════════════════════

cat("═══ PART 3: Spatial Co-occurrence ══════════════\n\n")

# --- 3a: Get TLS-enriched states ---
tls_enriched <- multi_results %>%
  filter(!is_interaction, FDR < 0.1, estimate > 1) %>%
  pull(term)

cat("TLS-enriched states:\n")
cat(paste(" ", tls_enriched, collapse = "\n"), "\n\n")

# --- 3b: Build spot-level abundance matrix ---
df_abund <- map_dfr(names(combined_obj), function(samp) {
  cs <- combined_obj[[samp]]$cellstate_norm
  if (is.null(cs)) return(NULL)
  
  cols <- intersect(tls_enriched, colnames(cs))
  if (length(cols) == 0) return(NULL)
  
  cs %>%
    select(ID, Label, all_of(cols)) %>%
    mutate(Sample = samp)
})

cat("Spots:", nrow(df_abund), "| Samples:", n_distinct(df_abund$Sample), "\n\n")

# --- 3c: METHOD 1 — Spearman correlation (per sample, then average) ---
cat("── Method 1: Correlation ────────────────────────\n")

cols_cor <- intersect(tls_enriched, colnames(df_abund))

cor_list <- df_abund %>%
  group_by(Sample) %>%
  group_map(~ {
    mat <- .x %>% select(all_of(cols_cor)) %>% as.matrix()
    cor(mat, use = "pairwise.complete.obs", method = "spearman")
  })

cor_avg <- Reduce("+", cor_list) / length(cor_list)
rownames(cor_avg) <- gsub("\\.", " ", rownames(cor_avg))
colnames(cor_avg) <- gsub("\\.", " ", colnames(cor_avg))

ht_cor <- Heatmap(
  cor_avg,
  name = "Spearman ρ",
  col  = colorRamp2(c(-0.3, 0, 0.3, 0.6),
                    c("#4DBBD5", "white", "#F39B7F", "#E64B35")),
  
  cluster_rows = TRUE, cluster_columns = TRUE,
  
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", cor_avg[i, j]), x, y,
              gp = gpar(fontsize = 8))
  },
  
  column_title    = "Spot-level Correlation: TLS-Enriched States",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_names_gp    = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  border = TRUE,
  heatmap_legend_param = list(direction = "horizontal")
)

pdf("../output/03_fig3c_cooccurrence_correlation.pdf", width = 8, height = 7)
draw(ht_cor, heatmap_legend_side = "right",
     padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()
cat("  Saved correlation heatmap\n\n")

# ════════════════════════════════════════════════════════════
# PART 4: CD4.T S01 vs S02 — Functional Programs via
#         Correlation in TLS Core vs Edge
# ════════════════════════════════════════════════════════════
# --- 4.0: Build TLS Core/Edge distance bins ---
# (Self-contained; does not depend on 02b environment)

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
    summarise(cx = mean(array_x), cy = mean(array_y),
              n_spots = n(), .groups = "drop")
  list(centroids = centroids, tls_spots = tls_spots)
}

cd4_states <- c("CD4.T_S01", "CD4.T_S02")

df_gradient <- purrr::map_dfr(names(combined_obj), function(samp) {
  meta <- combined_obj[[samp]]$cellstate_norm
  if (is.null(meta)) return(NULL)
  if (sum(meta$Label == "TLS") == 0) return(NULL)
  
  tls_info <- compute_tls_centroids(meta)
  if (is.null(tls_info)) return(NULL)
  centroids <- tls_info$centroids
  
  spot_coords <- as.matrix(meta[, c("array_x", "array_y")])
  centroid_coords <- as.matrix(centroids[, c("cx", "cy")])
  
  dist_to_centroids <- sapply(seq_len(nrow(centroid_coords)), function(k) {
    sqrt((spot_coords[, 1] - centroid_coords[k, 1])^2 +
           (spot_coords[, 2] - centroid_coords[k, 2])^2)
  })
  
  if (is.null(dim(dist_to_centroids))) {
    min_dist <- dist_to_centroids
  } else {
    min_dist <- apply(dist_to_centroids, 1, min)
  }
  
  dist_um <- min_dist * 100
  available_cols <- base::intersect(cd4_states, colnames(meta))
  
  meta %>%
    dplyr::select(ID, array_x, array_y, Label,
                  all_of(available_cols)) %>%
    mutate(Sample = samp,
           dist_to_centroid_um = dist_um,
           is_TLS = Label == "TLS")
})

# Assign TLS Core / Edge bins
tls_dist <- df_gradient %>% filter(is_TLS) %>% pull(dist_to_centroid_um)
tls_median <- median(tls_dist)

non_tls_dist <- df_gradient %>% filter(!is_TLS) %>% pull(dist_to_centroid_um)
q_vals <- quantile(non_tls_dist, probs = c(0.25, 0.5, 0.75))

df_gradient <- df_gradient %>%
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

cat("── df_gradient for PART 4 ──────────────────────\n")
cat("Total spots:", nrow(df_gradient), "\n")
cat("TLS Core   :", sum(df_gradient$dist_bin == "TLS Core"), "\n")
cat("TLS Edge   :", sum(df_gradient$dist_bin == "TLS Edge"), "\n\n")
# --- 4a: Extract marker genes ---
s01_name <- to_marker_name("CD4.T", "S01")
s02_name <- to_marker_name("CD4.T", "S02")

s01_markers <- all_marker_gene %>% filter(cell_state == s01_name)
s02_markers <- all_marker_gene %>% filter(cell_state == s02_name)

s01_genes <- s01_markers$Gene
s02_genes <- s02_markers$Gene

cat("  CD4.T_S01 markers:", length(s01_genes), "\n")
cat("  CD4.T_S02 markers:", length(s02_genes), "\n")
cat("  Overlap:", length(intersect(s01_genes, s02_genes)), "\n\n")

# --- 4b: Define functional categories ---
func_categories <- list(
  "Treg signature"                    = c("CTLA4", "IL2RA", "FOXP3", "IKZF2", "TNFRSF18"),
  "Immune checkpoint / co-regulation" = c("CD80", "CD86", "CD28", "CD226", "TIGIT", "ICOS"),
  "Tfh signature"                     = c("CD40LG", "ICOS", "SH2D1A", "PDCD1", "BCL6", "IL21"),
  "TLS organizing"                    = c("LTA", "LTB", "CCL19", "CCR7", "CXCL13", "CCL21"),
  "TCR complex / signaling"           = c("CD3G", "CD247", "TRAT1", "TESPA1", "THEMIS", "CD5"),
  "T cell activation"                 = c("CD69", "IL2", "IL7R", "IL2RB", "IL12RB1", "CD2"),
  "Chemokine / Migration"             = c("CXCR6", "CCR7", "CCR4", "CCR2", "GPR25", "SELPLG")
)

# --- 4c: Build gene annotation ---
all_func_genes <- unique(unlist(func_categories))
# Keep only genes that are actually S01 or S02 markers
all_s01s02_genes <- union(s01_genes, s02_genes)
target_genes <- intersect(all_func_genes, all_s01s02_genes)

gene_anno <- tibble(Gene = target_genes) %>%
  mutate(
    Marker_of = case_when(
      Gene %in% s01_genes & Gene %in% s02_genes ~ "Both",
      Gene %in% s01_genes ~ "S01",
      TRUE ~ "S02"
    )
  )

gene_anno$Category <- sapply(gene_anno$Gene, function(g) {
  for (cat_name in names(func_categories)) {
    if (g %in% func_categories[[cat_name]]) return(cat_name)
  }
  return("Other")
})

gene_anno <- gene_anno %>% filter(Category != "Other")
cat("  Target genes:", nrow(gene_anno), "\n\n")

# --- 4d: Get TLS Core/Edge spots with S01/S02 abundance ---
core_edge_spots <- df_gradient %>%
  filter(dist_bin %in% c("TLS Core", "TLS Edge")) %>%
  select(ID, Sample, dist_bin, CD4.T_S01, CD4.T_S02)

cat("  TLS Core spots:", sum(core_edge_spots$dist_bin == "TLS Core"), "\n")
cat("  TLS Edge spots:", sum(core_edge_spots$dist_bin == "TLS Edge"), "\n\n")

# --- 4e: Compute per-sample Spearman ρ for each gene × state × zone ---

cor_results <- map_dfr(unique(core_edge_spots$Sample), function(samp) {
  obj <- combined_obj[[samp]]
  if (is.null(obj)) return(NULL)
  hs <- obj$hs
  
  samp_spots <- core_edge_spots %>% filter(Sample == samp)
  spot_ids   <- intersect(samp_spots$ID, colnames(hs))
  if (length(spot_ids) < 10) return(NULL)  # need enough spots
  
  genes_present <- intersect(target_genes, rownames(hs))
  if (length(genes_present) == 0) return(NULL)
  
  # Expression matrix
  expr_mat <- as.matrix(hs[genes_present, spot_ids, drop = FALSE])
  
  # Match spots
  samp_spots_matched <- samp_spots[match(spot_ids, samp_spots$ID), ]
  
  # For each zone × state × gene: Spearman ρ
  map_dfr(c("TLS Core", "TLS Edge"), function(zone) {
    zone_idx <- which(samp_spots_matched$dist_bin == zone)
    if (length(zone_idx) < 5) return(NULL)  # minimum spots
    
    s01_abund <- samp_spots_matched$CD4.T_S01[zone_idx]
    s02_abund <- samp_spots_matched$CD4.T_S02[zone_idx]
    
    map_dfr(genes_present, function(g) {
      gene_expr <- expr_mat[g, zone_idx]
      
      # Spearman with S01
      rho_s01 <- tryCatch(
        cor(gene_expr, s01_abund, method = "spearman", use = "complete.obs"),
        error = function(e) NA_real_
      )
      # Spearman with S02
      rho_s02 <- tryCatch(
        cor(gene_expr, s02_abund, method = "spearman", use = "complete.obs"),
        error = function(e) NA_real_
      )
      
      tibble(
        Gene   = g,
        zone   = zone,
        rho_S01 = rho_s01,
        rho_S02 = rho_s02,
        Sample = samp
      )
    })
  })
})

cat("  Correlation records:", nrow(cor_results), "\n\n")

# --- 4f: Aggregate by category ---
# Add category annotation
cor_annotated <- cor_results %>%
  inner_join(gene_anno %>% select(Gene, Category), by = "Gene")

# Average ρ: first across genes within category per sample,
# then across samples
cat_cor <- cor_annotated %>%
  # Per sample, per category, per zone: mean ρ across genes
  group_by(Category, zone, Sample) %>%
  summarise(
    mean_rho_S01 = mean(rho_S01, na.rm = TRUE),
    mean_rho_S02 = mean(rho_S02, na.rm = TRUE),
    n_genes      = n_distinct(Gene),
    .groups      = "drop"
  ) %>%
  # Across samples: mean and SE
  group_by(Category, zone) %>%
  summarise(
    rho_S01     = mean(mean_rho_S01, na.rm = TRUE),
    rho_S01_se  = sd(mean_rho_S01, na.rm = TRUE) / sqrt(n()),
    rho_S02     = mean(mean_rho_S02, na.rm = TRUE),
    rho_S02_se  = sd(mean_rho_S02, na.rm = TRUE) / sqrt(n()),
    n_genes     = first(n_genes),
    n_samples   = n(),
    .groups     = "drop"
  )

cat("── Category-level correlations ─────────────────\n")
print(cat_cor %>% arrange(Category, zone), n = 30)
cat("\n")

# --- 4g: Build heatmap matrix ---
# Pivot to: rows = Category, cols = S01_Core, S01_Edge, S02_Core, S02_Edge

heat_long <- cat_cor %>%
  pivot_longer(
    cols = c(rho_S01, rho_S02),
    names_to = "state", values_to = "rho"
  ) %>%
  mutate(
    state = gsub("rho_", "", state),
    col_id = paste0(state, "_", gsub(" ", "_", zone))
  ) %>%
  select(Category, col_id, rho)

heat_mat <- heat_long %>%
  pivot_wider(names_from = col_id, values_from = rho) %>%
  tibble::column_to_rownames("Category") %>%
  as.matrix()

# Column and row order
desired_cols <- c("S01_TLS_Core", "S01_TLS_Edge", "S02_TLS_Core", "S02_TLS_Edge")
desired_cols <- intersect(desired_cols, colnames(heat_mat))

cat_order <- c(
  "Chemokine / Migration",
  "T cell activation",
  "Immune checkpoint / co-regulation",
  "TCR complex / signaling",
  "TLS organizing",
  "Tfh signature",
  "Treg signature"
)
cat_order <- intersect(cat_order, rownames(heat_mat))

heat_mat <- heat_mat[rev(cat_order), desired_cols, drop = FALSE]
heat_mat[is.na(heat_mat)] <- 0

cat("  Heatmap dimensions:", nrow(heat_mat), "x", ncol(heat_mat), "\n\n")

# --- 4h: Draw heatmap (raw ρ values, no z-score) ---

# Column annotation
col_state <- gsub("_TLS.*", "", desired_cols)
col_zone  <- gsub(".*TLS_", "TLS ", desired_cols)

col_ha_hm <- HeatmapAnnotation(
  State = col_state,
  Zone  = col_zone,
  col = list(
    State = c("S01" = "#7B68EE", "S02" = "#2E8B57"),
    Zone  = c("TLS Core" = "#C0392B", "TLS Edge" = "#F39B7F")
  ),
  show_annotation_name = TRUE,
  annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
  gap = unit(1, "mm")
)

# Bottom annotation for column labels
col_label_text <- c("S01\nCore", "S01\nEdge", "S02\nCore", "S02\nEdge")

bot_ha <- HeatmapAnnotation(
  label = anno_text(col_label_text,
                    gp = gpar(fontsize = 10, fontface = "bold"),
                    just = "center",
                    location = unit(0.5, "npc")),
  which = "column",
  annotation_name_side = "left",
  show_annotation_name = FALSE
)

# Color: white to dark red, scaled to data range
rho_min <- min(heat_mat, na.rm = TRUE)
rho_max <- max(heat_mat, na.rm = TRUE)

cat("  ρ range:", round(rho_min, 3), "to", round(rho_max, 3), "\n")

ht_func <- Heatmap(
  heat_mat,
  name = "Spearman ρ",
  col  = colorRamp2(c(rho_min, (rho_min + rho_max) / 2, rho_max),
                    c("white", "#F39B7F", "#B2182B")),
  
  top_annotation    = col_ha_hm,
  bottom_annotation = bot_ha,
  
  cluster_columns = FALSE,
  cluster_rows    = FALSE,
  
  # Show ρ values in cells
  cell_fun = function(j, i, x, y, w, h, fill) {
    rho_val <- heat_mat[i, j]
    grid.text(sprintf("%.2f", rho_val), x, y,
              gp = gpar(fontsize = 8,
                        col = ifelse(rho_val > (rho_min + rho_max * 2) / 3,
                                     "white", "black")))
  },
  
  show_column_names = FALSE,
  row_names_gp      = gpar(fontsize = 9),
  
  column_title     = "CD4.T S01 vs S02: Correlation with Functional Programs\nin TLS Core vs Edge",
  column_title_gp  = gpar(fontsize = 12, fontface = "bold"),
  
  border = TRUE,
  heatmap_legend_param = list(direction = "horizontal"),
  
  width  = unit(8, "cm"),
  height = unit(nrow(heat_mat) * 0.8, "cm")
)

pdf("../output/03_fig3d_CD4T_S01S02_Core_Edge_correlation.pdf",
    width = 9, height = 7)
draw(ht_func,
     heatmap_legend_side    = "bottom",
     annotation_legend_side = "right",
     padding = unit(c(10, 15, 25, 10), "mm"))
dev.off()
cat("  Saved: correlation heatmap\n\n")

# ════════════════════════════════════════════════════════════
# DONE
# ════════════════════════════════════════════════════════════
