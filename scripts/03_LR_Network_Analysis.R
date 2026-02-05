# ============================================================
# 03_LR_Network_Analysis.R
#
# Purpose:
#   Analyze ligand-receptor interaction networks across three key 
#   chemokine axes in TLS: CXCL13-CXCR5, CXCL12-CXCR4, and CCL19/21-CCR7.
#   Identify hub cell states, characterize network organization through
#   community detection, and assess functional enrichment patterns.
#
#   This script includes:
#   1) Cross-sample hubness analysis (network prominence bubble plots)
#   2) Hub community identification via Louvain clustering
#   3) Functional enrichment analysis (GO BP terms) for hub communities
#   4) Directional interaction heatmaps (sender-receiver relationships)
#
# Input:
#   (1) LR_analysis.RData(optional,If the computing takes too long time):
#       ../data/example_data/LR_analysis.RData
#       Contains pre-computed ligand-receptor analysis results for 
#       three chemokine axes. Each axis element includes:
#         - $graph: igraph object with interaction network
#         - $score: interaction score matrix (sender × receiver)
#         - $node_metrics: degree and weighted-degree metrics per cell state
#
#   (2) combined_obj.RData:
#       ../data/example_data/combined_obj.RData
#       Named list of Visium samples with TLS annotations and Cell state, Ecotype information
#
#   (3) all_marker_gene.rds:
#       ../data/example_data/all_marker_gene.rds
#       Marker genes defined by EcoTyper for each cell state (for GO enrichment)
#
# Output:
#   - Network prominence plots (03_*_Network_Prominence.pdf, 3 files)
#   - Hub community chord diagrams (03_*_Hub_community.pdf, 3 files)
#   - GO enrichment bubble plot (displayed, not saved in example)
#   - Significant interactions heatmaps (returned as list objects)
#
# =============================================================================
# 1. PACKAGE LOADING
# =============================================================================

library(conflicted)
library(clusterProfiler)
library(dplyr)
library(forcats)
library(ggplot2)
library(ggraph)
library(igraph)
library(org.Hs.eg.db)
library(purrr)
library(pheatmap)
library(stringr)
library(tidyr)
library(tibble)
library(viridis)
library(circlize)
library(RColorBrewer)

# Resolve function conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("setdiff", "base")
conflict_prefer("intersect", "base")
conflict_prefer("V", "igraph")
conflict_prefer("E", "igraph")
conflict_prefer("strength", "igraph")
conflict_prefer("degree", "igraph")
# =============================================================================
# 2. DATA LOADING AND REORGANIZATION
# =============================================================================

# Load individual LR analysis results
load("../data/example_data/combined_obj.RData")
load("../data/example_data/LR_analysis.RData")

# =============================================================================
# 3. CROSS-SAMPLE ANALYSIS FUNCTIONS
# =============================================================================

# -----------------------------------------------------------------------------
# 3.1 Hubness bubble plot (improved version)
# -----------------------------------------------------------------------------

plot_hubness_analysis <- function(axis_data, axis_name, 
                                  percentile_filter = 0.6) {
  
  # Shorten cell state names
  axis_data_modified <- map(axis_data, function(x) {
    V(x$graph)$name <- gsub("Monocytes\\.and\\.Macrophages", 
                            "Macrophages", 
                            V(x$graph)$name)
    return(x)
  })
  
  # Collect node metrics across samples
  node_tbl <- map_dfr(names(axis_data_modified), function(sid) {
    g <- axis_data_modified[[sid]]$graph
    tibble::tibble(
      sample = sid,
      cell_state = V(g)$name,
      wdeg = igraph::strength(g, mode = "all", weights = E(g)$score),
      edge_n = igraph::degree(g, mode = "all")
    )
  })
  
  # Check if we have data
  if (nrow(node_tbl) == 0) {
    stop("No node data found for ", axis_name)
  }
  
  # Aggregate by cell state
  summary_tbl <- node_tbl %>% 
    dplyr::group_by(cell_state) %>% 
    dplyr::summarise(
      median_wdeg = median(wdeg, na.rm = TRUE),
      mean_wdeg = mean(wdeg, na.rm = TRUE),
      edge_freq = sum(edge_n, na.rm = TRUE),
      median_edge = median(edge_n, na.rm = TRUE),
      n_sample = n(),
      .groups = "drop"
    )
  
  # Filter low-activity cell states
  wdeg_threshold <- quantile(summary_tbl$median_wdeg, percentile_filter, na.rm = TRUE)
  edge_threshold <- quantile(summary_tbl$median_edge, percentile_filter, na.rm = TRUE)
  
  filtered_tbl <- summary_tbl %>%
    dplyr::filter(median_wdeg > wdeg_threshold | median_edge > edge_threshold) %>%
    dplyr::filter(median_wdeg > 0.01 & median_edge > 0.5)
  
  cat("\n=== Hubness Analysis:", axis_name, "===\n")
  cat("Original cell states:", nrow(summary_tbl), "\n")
  cat("Filtered cell states:", nrow(filtered_tbl), "\n")
  cat("Removed low-activity states:", nrow(summary_tbl) - nrow(filtered_tbl), "\n")
  
  # Create bubble plot
  p <- ggplot(filtered_tbl,
              aes(x = median_wdeg,
                  y = median_edge,
                  size = edge_freq,
                  colour = n_sample,
                  label = cell_state)) +
    geom_point(alpha = 0.8) +
    scale_size(
      range = c(4, 20), 
      guide = guide_legend(
        title = "Total Edge\nFrequency",
        title.theme = element_text(size = 14, face = "bold"),
        label.theme = element_text(size = 12)
      )
    ) +
    scale_colour_gradient(
      low = "#c6dbef",
      high = "#1f78b4",
      name = "# Samples",
      guide = guide_colorbar(
        title.theme = element_text(size = 14, face = "bold"),
        label.theme = element_text(size = 12)
      )
    ) +
    ggrepel::geom_text_repel(
      max.overlaps = Inf, 
      size = 4.5,
      color = "black",
      fontface = "bold",
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      segment.size = 0.3
    ) +
    labs(
      title = paste(axis_name, "Network Prominence Across All Samples"),
      subtitle = paste("Showing", nrow(filtered_tbl), 
                       "active cell states (low-activity states excluded)"),
      x = "Median Weighted-Degree (Total Interaction Score)",
      y = "Median Edge Frequency (Σ Send + Receive)"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.3),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  # Show top 10 hubs
  cat("\nTop 10 hub cell states:\n")
  top_10 <- filtered_tbl %>%
    dplyr::arrange(desc(median_wdeg)) %>%
    head(10)
  print(top_10[, c("cell_state", "median_wdeg", "median_edge", "n_sample")])
  
  return(list(
    plot = p,
    filtered_data = filtered_tbl,
    summary_data = summary_tbl
  ))
}

# =============================================================================
# 4. HUB COMMUNITY IDENTIFICATION
# =============================================================================

# -----------------------------------------------------------------------------
# 4.1 Identify hub modules
# -----------------------------------------------------------------------------

get_hub_modules <- function(axis_data, top_n = 15) {
  
  # Identify top hub cell states
  hubs <- axis_data %>%
    purrr::map_dfr("node_metrics") %>%
    dplyr::group_by(cell_state) %>%
    dplyr::summarise(med_deg = median(degree), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(med_deg)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(cell_state)
  
  if (length(hubs) < 2) {
    return(tibble::tibble(hub_state = hubs, module = 1))
  }
  
  # Calculate co-membership counts
  co <- matrix(0, length(hubs), length(hubs), dimnames = list(hubs, hubs))
  
  for (x in axis_data) {
    tb <- tibble::tibble(
      cell = igraph::V(x$graph)$name,
      com = igraph::V(x$graph)$community
    ) %>%
      dplyr::filter(cell %in% hubs)
    
    same <- dplyr::inner_join(tb, tb, by = "com", suffix = c("_A", "_B")) %>%
      dplyr::filter(cell_A < cell_B)
    
    for (i in seq_len(nrow(same))) {
      a <- same$cell_A[i]
      b <- same$cell_B[i]
      co[a, b] <- co[a, b] + 1
    }
  }
  
  co <- co + t(co)
  diag(co) <- 0
  
  # Module detection
  if (all(co == 0)) {
    mod_vec <- rep(1, length(hubs))
    names(mod_vec) <- hubs
  } else {
    # Jaccard similarity
    rs <- rowSums(co)
    jacc <- co / (outer(rs, rs, "+") - co + 1e-9)
    
    g <- igraph::graph_from_adjacency_matrix(jacc, "undirected",
                                             weighted = TRUE, diag = FALSE)
    mod_vec <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)$membership
    names(mod_vec) <- hubs
  }
  
  tibble::tibble(
    hub_state = names(mod_vec),
    module = unname(mod_vec)
  ) %>%
    dplyr::arrange(module, hub_state)
}

# -----------------------------------------------------------------------------
# 4.2 Chord diagram visualization
# -----------------------------------------------------------------------------

draw_hub_chord <- function(tbl, title,
                           start_deg = 130,
                           label_cex = 1.2,
                           label_rad = 1.4,
                           title_cex = 1.5) {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(oma = c(2, 2, 3, 2))
  
  # Ensure column names
  if (!"cell_state" %in% names(tbl)) names(tbl)[1L] <- "cell_state"
  if (!"module" %in% names(tbl)) names(tbl)[2L] <- "module"
  
  # Shorten cell state names
  tbl <- tbl %>%
    dplyr::mutate(cell_state = stringr::str_replace(
      cell_state, 
      "Monocytes\\.and\\.Macrophages", 
      "Macrophages"
    ))
  
  # Create edge list (Community -> cell state)
  edges <- tbl %>%
    dplyr::transmute(
      from = paste0("C", module),
      to = cell_state,
      value = 1
    )
  
  # Color setup
  mod_ids <- sort(unique(edges$from))
  n_mod <- length(mod_ids)
  
  if (n_mod <= 3) {
    mod_col <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A")[seq_len(n_mod)], mod_ids)
  } else {
    mod_col <- setNames(
      RColorBrewer::brewer.pal(max(3, n_mod), "Set2")[seq_len(n_mod)],
      mod_ids
    )
  }
  
  cs_ids <- setdiff(unique(edges$to), mod_ids)
  grid_col <- c(mod_col, setNames(rep("#CCCCCC", length(cs_ids)), cs_ids))
  
  # Circos parameters
  sectors <- c(mod_ids, cs_ids)
  gap.after <- c(rep(3, length(sectors) - 1), 10)
  
  circos.clear()
  circos.par(
    start.degree = start_deg,
    gap.after = gap.after,
    track.margin = c(0.02, 0.02),
    canvas.xlim = c(-1.3, 1.3),
    canvas.ylim = c(-1.3, 1.3)
  )
  
  # Draw chord diagram
  circlize::chordDiagram(
    edges,
    order = sectors,
    grid.col = grid_col,
    transparency = 0.2,
    directional = 0,
    annotationTrack = "grid",
    annotationTrackHeight = 0.04
  )
  
  # Add labels
  circos.trackPlotRegion(
    track.index = 1,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      
      # Bold for community labels
      font_weight <- ifelse(grepl("^C\\d+$", sector.index), 2, 1)
      
      circos.text(mean(xlim), label_rad, sector.index,
                  facing = "clockwise", 
                  niceFacing = TRUE,
                  cex = label_cex,
                  font = font_weight,
                  adj = c(0.5, 0))
    }
  )
  
  # Add title
  mtext(title, side = 3, line = 0.5, cex = title_cex, font = 2)
}

# =============================================================================
# 5. FUNCTIONAL ENRICHMENT ANALYSIS
# =============================================================================

# Load marker genes
gene_tbl <- readRDS("../data/example_data/all_marker_gene.rds")

# Function to create gene lists for each module
make_gene_lists <- function(hub_tbl, network_tag) {
  gene_tbl %>%
    dplyr::inner_join(hub_tbl, by = c(cell_state = "hub_state")) %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(gene_vec = list(unique(Gene)), .groups = "drop") %>%
    dplyr::mutate(module_tag = paste0(network_tag, "_M", module)) %>%
    dplyr::select(module_tag, gene_vec)
}

# Function to run GO enrichment
run_go_enrichment <- function(hub_modules_list) {
  
  # Combine gene lists from all axes
  gene_lists <- bind_rows(
    make_gene_lists(hub_modules_list$CXCL13_CXCR5, "CXCL13_CXCR5"),
    make_gene_lists(hub_modules_list$CXCL12_CXCR4, "CXCL12_CXCR4"),
    make_gene_lists(hub_modules_list$CCL19_21_CCR7, "CCL19_21_CCR7")
  )
  
  # Run GO enrichment
  go_results <- purrr::pmap_dfr(gene_lists, function(module_tag, gene_vec) {
    ego <- enrichGO(
      gene = gene_vec,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    if (is.null(ego) || nrow(ego) == 0) {
      return(tibble::tibble())
    }
    
    as_tibble(ego) %>% dplyr::mutate(module = module_tag)
  })
  
  return(go_results)
}

# Function to plot GO enrichment
plot_go_enrichment <- function(go_results, topn = 10) {
  
  # Ordered community levels
  ordered_levels <- c(
    "CXCL13_CXCR5_C1", "CXCL12_CXCR4_C1", "CCL19_21_CCR7_C1",
    "CXCL13_CXCR5_C2", "CXCL12_CXCR4_C2", "CCL19_21_CCR7_C2"
  )
  
  # Prepare plot data
  plot_df <- go_results %>% 
    dplyr::mutate(module = str_replace(module, "_M(\\d+)$", "_C\\1")) %>% 
    dplyr::filter(!grepl("adaptive immune response based on somatic recombination", 
                         Description, ignore.case = TRUE)) %>%
    dplyr::group_by(module) %>% 
    dplyr::slice_min(p.adjust, n = topn, with_ties = FALSE) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      module = factor(module, levels = rev(ordered_levels)),
      Description = str_wrap(Description, width = 40),
      community_type = ifelse(grepl("_C1$", module), "C1", "C2")
    )
  
  # Create bubble plot
  p <- ggplot(plot_df,
              aes(x = Description,
                  y = module,
                  size = Count,
                  colour = -log10(p.adjust))) +
    geom_point() +
    scale_colour_viridis_c(name = "–log10 FDR") +
    scale_size(range = c(3, 10), name = "Gene count") +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      x = "GO BP term",
      y = "Hub community",
      title = "Functional enrichment of hub communities"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, 
                                 size = 14, margin = margin(t = 5)),
      axis.text.y = element_text(
        size = 13,
        colour = ifelse(grepl("_C1$", rev(ordered_levels)), "#E41A1C", "#377EB8")
      ),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin = margin(10, 10, 30, 10),
      panel.grid.major.x = element_line(linewidth = 0.3),
      panel.background = element_rect(fill = "white", colour = NA)
    )
  
  return(p)
}

# =============================================================================
# 6. DIRECTIONAL INTERACTION HEATMAP
# =============================================================================

create_significant_interactions_heatmap <- function(axis_data, 
                                                    axis_name,
                                                    score_threshold_quantile = 0.8,
                                                    min_samples = 2,
                                                    top_n_pairs = 50) {
  
  cat("=== Creating Significant Interactions Heatmap:", axis_name, "===\n")
  cat("Parameters:\n")
  cat("  Score threshold quantile:", score_threshold_quantile, "\n")
  cat("  Minimum samples required:", min_samples, "\n")
  cat("  Top N pairs to show:", top_n_pairs, "\n\n")
  
  # Shorten cell state names
  axis_data <- map(axis_data, function(x) {
    rownames(x$score) <- gsub("Monocytes\\.and\\.Macrophages", 
                              "Macrophages", rownames(x$score))
    colnames(x$score) <- gsub("Monocytes\\.and\\.Macrophages", 
                              "Macrophages", colnames(x$score))
    return(x)
  })
  
  # Collect all cell state names
  cell_types <- unique(unlist(map(axis_data, ~ rownames(.x$score))))
  cat("Total cell states found:", length(cell_types), "\n")
  
  # Initialize matrices
  n <- length(cell_types)
  sum_mat <- matrix(0, n, n, dimnames = list(cell_types, cell_types))
  cnt_mat <- sum_mat
  all_values <- list()
  
  # Aggregate scores across samples
  cat("Processing", length(axis_data), "samples...\n")
  
  for (sample_id in names(axis_data)) {
    s <- axis_data[[sample_id]]$score
    
    sum_mat[rownames(s), colnames(s)] <- sum_mat[rownames(s), colnames(s)] + s
    cnt_mat[rownames(s), colnames(s)] <- cnt_mat[rownames(s), colnames(s)] + 1
    
    # Store individual values for variance calculation
    for (i in rownames(s)) {
      for (j in colnames(s)) {
        pair_name <- paste(i, j, sep = "->")
        if (is.null(all_values[[pair_name]])) {
          all_values[[pair_name]] <- c()
        }
        all_values[[pair_name]] <- c(all_values[[pair_name]], s[i, j])
      }
    }
  }
  
  # Calculate mean scores
  mean_mat <- sum_mat / pmax(cnt_mat, 1)
  mean_mat[cnt_mat == 0] <- NA
  
  # Statistical filtering
  sufficient_samples <- cnt_mat >= min_samples
  
  all_nonzero_scores <- as.vector(mean_mat[!is.na(mean_mat) & mean_mat > 0])
  score_threshold <- quantile(all_nonzero_scores, score_threshold_quantile)
  cat("Score threshold (", score_threshold_quantile*100, "th percentile):", 
      round(score_threshold, 3), "\n")
  
  high_score <- mean_mat >= score_threshold
  high_score[is.na(high_score)] <- FALSE
  
  # Statistical significance test
  significant_mat <- matrix(FALSE, n, n, dimnames = list(cell_types, cell_types))
  
  for (i in 1:n) {
    for (j in 1:n) {
      pair_name <- paste(cell_types[i], cell_types[j], sep = "->")
      if (pair_name %in% names(all_values) && 
          length(all_values[[pair_name]]) >= min_samples) {
        values <- all_values[[pair_name]]
        if (length(values) > 1 && var(values) > 0) {
          t_result <- t.test(values, mu = 0)
          significant_mat[i, j] <- t_result$p.value < 0.05
        }
      }
    }
  }
  
  # Combined filter
  final_filter <- sufficient_samples & high_score & significant_mat
  cat("Interactions passing all filters:", sum(final_filter, na.rm = TRUE), "\n")
  
  # Select top pairs
  significant_scores <- mean_mat * final_filter
  significant_scores[!final_filter] <- 0
  
  interaction_df <- data.frame(
    sender = rep(cell_types, each = n),
    receiver = rep(cell_types, times = n),
    score = as.vector(significant_scores),
    n_samples = as.vector(cnt_mat),
    significant = as.vector(final_filter),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(significant & score > 0) %>%
    dplyr::arrange(desc(score)) %>%
    head(top_n_pairs)
  
  cat("Top", nrow(interaction_df), "significant interactions selected\n")
  
  # Create filtered matrix
  involved_cells <- unique(c(interaction_df$sender, interaction_df$receiver))
  cat("Cell states involved in top interactions:", length(involved_cells), "\n")
  
  filtered_mat <- matrix(0, length(involved_cells), length(involved_cells),
                         dimnames = list(involved_cells, involved_cells))
  
  for (k in 1:nrow(interaction_df)) {
    sender <- interaction_df$sender[k]
    receiver <- interaction_df$receiver[k]
    score <- interaction_df$score[k]
    filtered_mat[sender, receiver] <- score
  }
  
  filtered_mat[filtered_mat == 0] <- NA
  filtered_mat[is.infinite(filtered_mat)] <- NA
  filtered_mat[!is.finite(filtered_mat) & !is.na(filtered_mat)] <- NA
  
  row_has_data <- apply(filtered_mat, 1, function(x) any(!is.na(x) & x > 0))
  col_has_data <- apply(filtered_mat, 2, function(x) any(!is.na(x) & x > 0))
  
  if (sum(row_has_data) == 0 || sum(col_has_data) == 0) {
    stop("No valid data remaining after filtering. Try relaxing thresholds.")
  }
  
  filtered_mat <- filtered_mat[row_has_data, col_has_data, drop = FALSE]
  cat("Matrix size after cleaning:", nrow(filtered_mat), "x", ncol(filtered_mat), "\n")
  
  # Create annotations
  sender_scores <- rowSums(filtered_mat, na.rm = TRUE)
  sender_scores[is.na(sender_scores) | is.infinite(sender_scores)] <- 0
  
  row_annotation <- data.frame(
    Total_Sent = sender_scores,
    row.names = names(sender_scores)
  )
  
  receiver_scores <- colSums(filtered_mat, na.rm = TRUE)
  receiver_scores[is.na(receiver_scores) | is.infinite(receiver_scores)] <- 0
  
  col_annotation <- data.frame(
    Total_Received = receiver_scores,
    row.names = names(receiver_scores)
  )
  
  annotation_colors <- list(
    Total_Sent = colorRampPalette(c("white", "coral", "red"))(100),
    Total_Received = colorRampPalette(c("white", "lightblue", "darkblue"))(100)
  )
  
  # Check clustering feasibility
  can_cluster_rows <- nrow(filtered_mat) > 1 && 
    sum(!is.na(filtered_mat)) > nrow(filtered_mat)
  can_cluster_cols <- ncol(filtered_mat) > 1 && 
    sum(!is.na(filtered_mat)) > ncol(filtered_mat)
  
  cat("Clustering options: rows =", can_cluster_rows, ", cols =", can_cluster_cols, "\n")
  
  # Create labels
  lab_row <- paste0("S: ", rownames(filtered_mat))
  lab_col <- paste0("R: ", colnames(filtered_mat))
  dimnames(filtered_mat) <- list(lab_row, lab_col)
  
  # Draw heatmap
  tryCatch({
    p <- pheatmap(
      filtered_mat,
      color = colorRampPalette(c("#f7fbff", "#6baed6", "#2171b5", "#08306b"))(100),
      na_col = "grey95",
      cluster_rows = can_cluster_rows,
      cluster_cols = can_cluster_cols,
      annotation_row = if(can_cluster_rows) row_annotation else NULL,
      annotation_col = if(can_cluster_cols) col_annotation else NULL,
      annotation_colors = annotation_colors,
      legend = FALSE,
      main = paste0("Significant ", axis_name, " Interactions\n",
                    "(p<0.05, >", score_threshold_quantile*100, 
                    "th percentile)\nS: sender / R: receiver"),
      fontsize_row = 8,
      fontsize_col = 8,
      angle_col = 45,
      border_color = "white",
      cellwidth = 15,
      cellheight = 15,
      silent = TRUE
    )
  }, error = function(e) {
    cat("Clustering failed, creating simple heatmap...\n")
    p <<- pheatmap(
      filtered_mat,
      color = colorRampPalette(c("#f7fbff", "#6baed6", "#2171b5", "#08306b"))(100),
      na_col = "grey95",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      legend = FALSE,
      main = paste0("Significant ", axis_name, " Interactions\n",
                    "(p<0.05, >", score_threshold_quantile*100, 
                    "th percentile)\nS: sender / R: receiver"),
      fontsize_row = 8,
      fontsize_col = 8,
      angle_col = 45,
      border_color = "white",
      cellwidth = 15,
      cellheight = 15,
      silent = TRUE
    )
  })
  
  # Print summary
  cat("\n=== Results Summary ===\n")
  cat("Final selected pairs:", nrow(interaction_df), "\n")
  cat("Final matrix size:", nrow(filtered_mat), "x", ncol(filtered_mat), "\n")
  
  cat("\nTop 10 interactions:\n")
  print(head(interaction_df[, c("sender", "receiver", "score", "n_samples")], 10))
  
  return(list(
    heatmap = p,
    interaction_df = interaction_df,
    filtered_matrix = filtered_mat,
    involved_cells = rownames(filtered_mat)
  ))
}

# =============================================================================
# 7. EXAMPLE USAGE
# =============================================================================

# Run hubness analysis for all axes
hubness_results <- list(
  CXCL13_CXCR5 = plot_hubness_analysis(LR_analysis$CXCL13_CXCR5, "CXCL13–CXCR5"),
  CXCL12_CXCR4 = plot_hubness_analysis(LR_analysis$CXCL12_CXCR4, "CXCL12–CXCR4"),
  CCL19_21_CCR7 = plot_hubness_analysis(LR_analysis$CCL19_21_CCR7, "CCL19/21–CCR7")
)

# Display plots
suppressWarnings({
  ggsave(hubness_results$CXCL13_CXCR5$plot, file = "../output/03_CXCL13-CXCR5_Network_Prominence.pdf",
         width      = 10,
         height     = 8,
         units      = "in",
         dpi        = 300,
         limitsize  = FALSE)  
  ggsave(hubness_results$CXCL12_CXCR4$plot, file = "../output/03_CXCL12_CXCR4_Network_Prominence.pdf",
         width      = 10,
         height     = 8,
         units      = "in",
         dpi        = 300,
         limitsize  = FALSE)  
  ggsave(hubness_results$CCL19_21_CCR7$plot, file = "../output/03_CCL19_21_CCR7_Network_Prominence.pdf",
         width      = 10,
         height     = 8,
         units      = "in",
         dpi        = 300,
         limitsize  = FALSE)  
})


# Identify hub modules
hub_modules <- list(
  CXCL13_CXCR5 = get_hub_modules(LR_analysis$CXCL13_CXCR5, top_n = 20),
  CXCL12_CXCR4 = get_hub_modules(LR_analysis$CXCL12_CXCR4, top_n = 15),
  CCL19_21_CCR7 = get_hub_modules(LR_analysis$CCL19_21_CCR7, top_n = 20)
)

# Draw chord diagrams
p1 <- suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CXCL13_CXCR5, title = "CXCL13–CXCR5: Hub Community to Cell State")))
p2 <- suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CXCL12_CXCR4, title = "CXCL12–CXCR4: Hub Community to Cell State")))
p3 <- suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CCL19_21_CCR7, title = "CCL19/21–CCR7: Hub Community to Cell State")))

# Display plots
# CXCL13–CXCR5
pdf(file   = "../output/03_CXCL13-CXCR5_Hub_community.pdf", width  = 10,height = 8)
suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CXCL13_CXCR5,title = "CXCL13–CXCR5: Hub Community to Cell State")))
dev.off()


# CXCL12–CXCR4
pdf(file   = "../output/03_CXCL12_CXCR4_Hub_community.pdf",width  = 10,height = 8)
suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CXCL12_CXCR4,title = "CXCL12–CXCR4: Hub Community to Cell State")))
dev.off()

# CCL19/21–CCR7
pdf(file   = "../output/03_CCL19_21_CCR7_Hub_community.pdf",width  = 10,height = 8)
suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CCL19_21_CCR7,title = "CCL19/21–CCR7: Hub Community to Cell State")))
dev.off()


# Run GO enrichment
go_results <- run_go_enrichment(hub_modules)

# Plot GO enrichment
go_plot <- plot_go_enrichment(go_results, topn = 10)
print(go_plot)

# Generate directional interaction heatmaps
interaction_results <- list(
  CXCL13_CXCR5 = create_significant_interactions_heatmap(
    LR_analysis$CXCL13_CXCR5, "CXCL13–CXCR5"
  ),
  CXCL12_CXCR4 = create_significant_interactions_heatmap(
    LR_analysis$CXCL12_CXCR4, "CXCL12–CXCR4"
  ),
  CCL19_21_CCR7 = create_significant_interactions_heatmap(
    LR_analysis$CCL19_21_CCR7, "CCL19/21–CCR7"
  )
)

# =============================================================================
# END OF SCRIPT
# =============================================================================