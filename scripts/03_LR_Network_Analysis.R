# =============================================================================
# 03_LR_Network_Analysis.R
#
# Purpose:
#   Analyze ligand-receptor interaction networks across three key
#   chemokine axes in TLS: CXCL13-CXCR5, CXCL12-CXCR4, and CCL19/21-CCR7.
#
#   This script includes:
#   1) Ligand-receptor network computation (spatially weighted)
#   2) Cross-sample hubness analysis (network prominence bubble plots)
#   3) Hub community identification via Louvain clustering
#   4) Functional enrichment analysis (GO BP terms) for hub communities
#   5) Directional interaction heatmaps (sender-receiver relationships)
#   6) Hub cell state interaction partner lineage diversity
#   7) C2 hub community members: TLS vs NonTLS comparison
#
# Input:
#   (1) LR_analysis.RData (optional, if computation takes too long):
#       ../data/example_data/LR_analysis.RData
#       Contains pre-computed ligand-receptor analysis results for
#       three chemokine axes. Each axis element includes:
#         - $graph: igraph object with interaction network
#         - $score: interaction score matrix (sender x receiver)
#         - $node_metrics: degree and weighted-degree metrics per cell state
#
#   (2) combined_obj.RData:
#       ../data/example_data/combined_obj.RData
#       Named list of Visium samples with TLS annotations and
#       cell state / ecotype information
#
#   (3) all_marker_gene.rds:
#       ../data/example_data/all_marker_gene.rds
#       Marker genes defined by EcoTyper for each cell state
#
# Output:
#   - 03_*_Network_Prominence.pdf (3 files, one per axis)
#   - 03_*_Hub_community.pdf (3 files, chord diagrams)
#   - 03_GO_enrichment_hub_communities.pdf
#   - 03_Hub_Lineage_Diversity_C1.pdf
#   - 03_Hub_Lineage_Diversity_C2.pdf
#   - 03_C2_TLS_vs_NonTLS.pdf
#
# Dependencies:
#   conflicted, clusterProfiler, dplyr, forcats, ggplot2, ggraph,
#   igraph, org.Hs.eg.db, purrr, pheatmap, stringr, tidyr, tibble,
#   viridis, circlize, RColorBrewer, ggrepel
# =============================================================================

# =============================================================================
# SECTION 1: PACKAGE LOADING
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
library(ggrepel)

# Resolve function conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(igraph::V)
conflicts_prefer(igraph::E)
conflicts_prefer(igraph::strength)
conflicts_prefer(igraph::degree)

# =============================================================================
# SECTION 2: DATA LOADING
# =============================================================================

load("../data/example_data/combined_obj.RData")
load("../data/example_data/LR_analysis.RData")

# =============================================================================
# SECTION 3: LIGAND-RECEPTOR NETWORK COMPUTATION
# =============================================================================
# This section defines the core analysis function and runs it for all axes.
# Pre-computed results are available in LR_analysis.RData;
# re-running is only needed if parameters change.

analyze_LR_network <- function(combined_obj,
                               ligand_gene,
                               receptor_gene,
                               lambda = 2,
                               quant_thr = 0.9,
                               sample_id = NULL) {
  
  # If sample_id provided, analyze single sample; otherwise analyze all
  samples <- if (!is.null(sample_id)) sample_id else names(combined_obj)
  
  results <- list()
  
  for (sid in samples) {
    message("[Processing] ", sid)
    
    # Extract metadata and filter for TLS regions
    meta <- combined_obj[[sid]]$cellstate_norm
    rownames(meta) <- meta$ID
    tls_ids <- meta$ID[meta$Label == "TLS"]
    
    if (length(tls_ids) == 0) {
      message("  No TLS spots, skipping")
      next
    }
    
    cell_states <- grep("_S\\d{2}$", colnames(meta), value = TRUE)
    
    # Extract expression matrix and coordinates
    expr <- as.matrix(combined_obj[[sid]]$hs)[, tls_ids]
    meta_t <- meta[tls_ids, ]
    cellmat <- as.matrix(meta_t[, cell_states])
    coords <- as.matrix(meta_t[, c("array_x", "array_y")])
    
    # Check gene presence
    all_genes <- c(ligand_gene, receptor_gene)
    if (!all(all_genes %in% rownames(expr))) {
      missing <- setdiff(all_genes, rownames(expr))
      message("  ", paste(missing, collapse = ", "), " absent, skipping")
      next
    }
    
    # Handle single or multiple ligands
    if (length(ligand_gene) == 1) {
      ligand <- as.numeric(expr[ligand_gene, ])
    } else {
      ligand <- colMeans(expr[ligand_gene, , drop = FALSE])
    }
    receptor <- as.numeric(expr[receptor_gene, ])
    
    # Skip if no expression
    if (sum(ligand) == 0 || sum(receptor) == 0) {
      message("  Ligand or receptor expression all zero, skipping")
      next
    }
    
    # Calculate spatially-weighted interaction scores
    wmat <- exp(-as.matrix(dist(coords)) / lambda)
    ncs <- length(cell_states)
    score <- matrix(0, ncs, ncs, dimnames = list(cell_states, cell_states))
    
    for (i in seq_len(ncs)) {
      for (j in seq_len(ncs)) {
        send <- ligand * cellmat[, i]
        recv <- receptor * cellmat[, j]
        score[i, j] <- sum(outer(send, recv) * wmat)
      }
    }
    
    # Build network graph
    edges <- as.data.frame(as.table(score))
    names(edges) <- c("sender", "receiver", "score")
    edges <- edges[!is.na(edges$score), ]
    if (nrow(edges) == 0) {
      message("  All scores are NA, skipping")
      next
    }
    edges <- edges[edges$score >= quantile(edges$score, quant_thr, na.rm = TRUE), ]
    
    g <- igraph::graph_from_data_frame(edges, directed = TRUE)
    
    # Community detection
    if (vcount(g) > 2) {
      mem <- igraph::cluster_louvain(igraph::as.undirected(g))$membership
    } else {
      mem <- rep(1, vcount(g))
    }
    V(g)$community <- mem
    
    # Calculate centrality metrics
    V(g)$degree <- degree(g, mode = "all")
    V(g)$betweenness <- betweenness(g, directed = TRUE)
    V(g)$pagerank <- page.rank(g)$vector
    
    node_metrics <- data.frame(
      cell_state  = V(g)$name,
      degree      = degree(g, mode = "all"),
      betweenness = betweenness(g, directed = TRUE),
      pagerank    = page.rank(g)$vector,
      community   = as.factor(mem)
    )
    
    results[[sid]] <- list(
      score        = score,
      edge_list    = edges,
      node_metrics = node_metrics,
      graph        = g
    )
  }
  
  return(results)
}


# VISUALIZATION FUNCTION FOR SINGLE SAMPLE
plot_LR_network <- function(graph_obj, ligand_gene, receptor_gene,
                            sample_id = NULL, enhanced = FALSE) {
  
  g <- graph_obj$graph
  ligand_name <- paste(ligand_gene, collapse = "/")
  
  title <- if (!is.null(sample_id)) {
    paste0("Network of ", ligand_name, "-", receptor_gene,
           " Communication in ", sample_id)
  } else {
    paste0("Network of ", ligand_name, "-", receptor_gene, " Communication")
  }
  
  if (enhanced) {
    p <- ggraph(g, layout = "kk") +
      geom_edge_link(aes(edge_alpha = score, edge_width = score,
                         edge_color = score),
                     arrow = arrow(length = unit(6, "mm")),
                     end_cap = circle(5, "mm"),
                     edge_curved = 0.2) +
      geom_node_point(aes(size = degree, color = as.factor(community)),
                      show.legend = TRUE) +
      geom_node_text(aes(label = paste0(name, "(", degree, ")")),
                     repel = TRUE, size = 6, fontface = "bold") +
      scale_size_continuous(range = c(6, 18)) +
      scale_edge_color_gradient(low = "lightblue", high = "darkblue",
                                name = "Interaction score") +
      scale_edge_width(range = c(0.5, 4), name = "Interaction score") +
      scale_edge_alpha(range = c(0.4, 1), name = "Interaction score") +
      theme_graph(base_size = 16) +
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 12)
      ) +
      ggtitle(title)
  } else {
    p <- ggraph(g, layout = "kk") +
      geom_edge_link(aes(edge_alpha = score, edge_width = score,
                         edge_color = score),
                     arrow = arrow(length = unit(4, "mm")),
                     end_cap = circle(3, "mm"),
                     edge_curved = 0.2) +
      geom_node_point(aes(size = degree, color = as.factor(community)),
                      show.legend = TRUE) +
      geom_node_text(aes(label = paste0(name, "(", degree, ")")),
                     repel = TRUE, size = 4) +
      scale_size_continuous(range = c(4, 12)) +
      scale_edge_color_gradient(low = "lightblue", high = "darkblue",
                                name = "Interaction score") +
      scale_edge_width(range = c(0.2, 3), name = "Interaction score") +
      scale_edge_alpha(range = c(0.3, 1), name = "Interaction score") +
      theme_graph() +
      ggtitle(title)
  }
  
  return(p)
}


# --- Single sample examples (optional, uses first available sample) ---
# Uncomment to run on a specific sample:
# example_sample <- names(combined_obj)[1]
# result_example <- analyze_LR_network(
#   combined_obj, "CXCL13", "CXCR5", sample_id = example_sample
# )
# print(plot_LR_network(result_example[[example_sample]],
#                       "CXCL13", "CXCR5", sample_id = example_sample))


# --- Batch processing for all three chemokine axes ---
# (Skip if using pre-computed LR_analysis.RData)

# message("\n========== Analyzing CXCL13-CXCR5 axis ==========")
# LR_analysis$CXCL13_CXCR5 <- analyze_LR_network(
#   combined_obj, ligand_gene = "CXCL13", receptor_gene = "CXCR5"
# )
#
# message("\n========== Analyzing CCL19/21-CCR7 axis ==========")
# LR_analysis$CCL19_21_CCR7 <- analyze_LR_network(
#   combined_obj, ligand_gene = c("CCL19", "CCL21"), receptor_gene = "CCR7"
# )
#
# message("\n========== Analyzing CXCL12-CXCR4 axis ==========")
# LR_analysis$CXCL12_CXCR4 <- analyze_LR_network(
#   combined_obj, ligand_gene = "CXCL12", receptor_gene = "CXCR4"
# )


# =============================================================================
# SECTION 4: CROSS-SAMPLE HUBNESS ANALYSIS
# =============================================================================

plot_hubness_analysis <- function(axis_data, axis_name,
                                  percentile_filter = 0.6) {
  
  # Shorten cell state names for display
  axis_data_modified <- map(axis_data, function(x) {
    V(x$graph)$name <- gsub("Monocytes\\.and\\.Macrophages",
                            "Macrophages",
                            V(x$graph)$name)
    return(x)
  })
  
  # Collect node metrics across samples
  node_tbl <- map_dfr(names(axis_data_modified), function(sid) {
    g <- axis_data_modified[[sid]]$graph
    tibble(
      sample     = sid,
      cell_state = V(g)$name,
      wdeg       = igraph::strength(g, mode = "all", weights = E(g)$score),
      edge_n     = igraph::degree(g, mode = "all")
    )
  })
  
  if (nrow(node_tbl) == 0) {
    stop("No node data found for ", axis_name)
  }
  
  # Aggregate by cell state
  summary_tbl <- node_tbl %>%
    group_by(cell_state) %>%
    summarise(
      median_wdeg = median(wdeg, na.rm = TRUE),
      mean_wdeg   = mean(wdeg, na.rm = TRUE),
      edge_freq   = sum(edge_n, na.rm = TRUE),
      median_edge = median(edge_n, na.rm = TRUE),
      n_sample    = n(),
      .groups     = "drop"
    )
  
  # Filter low-activity cell states
  wdeg_threshold <- quantile(summary_tbl$median_wdeg,
                             percentile_filter, na.rm = TRUE)
  edge_threshold <- quantile(summary_tbl$median_edge,
                             percentile_filter, na.rm = TRUE)
  
  filtered_tbl <- summary_tbl %>%
    filter(median_wdeg > wdeg_threshold | median_edge > edge_threshold) %>%
    filter(median_wdeg > 0.01 & median_edge > 0.5)
  
  cat("\n=== Hubness Analysis:", axis_name, "===\n")
  cat("Original cell states:", nrow(summary_tbl), "\n")
  cat("Filtered cell states:", nrow(filtered_tbl), "\n")
  
  # Bubble plot
  p <- ggplot(filtered_tbl,
              aes(x = median_wdeg, y = median_edge,
                  size = edge_freq, colour = n_sample,
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
      low = "#c6dbef", high = "#1f78b4",
      name = "# Samples",
      guide = guide_colorbar(
        title.theme = element_text(size = 14, face = "bold"),
        label.theme = element_text(size = 12)
      )
    ) +
    ggrepel::geom_text_repel(
      max.overlaps = Inf, size = 4.5, color = "black",
      fontface = "bold", box.padding = 0.5,
      point.padding = 0.3, segment.color = "gray50",
      segment.size = 0.3
    ) +
    labs(
      title    = paste(axis_name, "Network Prominence Across All Samples"),
      subtitle = paste("Showing", nrow(filtered_tbl),
                       "active cell states (low-activity states excluded)"),
      x = "Median Weighted-Degree (Total Interaction Score)",
      y = "Median Edge Frequency (Send + Receive)"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title    = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray40"),
      axis.title.x  = element_text(size = 16, face = "bold"),
      axis.title.y  = element_text(size = 16, face = "bold"),
      axis.text.x    = element_text(size = 13),
      axis.text.y    = element_text(size = 13),
      legend.position = "right",
      legend.box      = "vertical",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.3),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  cat("\nTop 10 hub cell states:\n")
  top_10 <- filtered_tbl %>% arrange(desc(median_wdeg)) %>% head(10)
  print(top_10[, c("cell_state", "median_wdeg", "median_edge", "n_sample")])
  
  return(list(
    plot          = p,
    filtered_data = filtered_tbl,
    summary_data  = summary_tbl
  ))
}

# =============================================================================
# SECTION 5: HUB COMMUNITY IDENTIFICATION
# =============================================================================

get_hub_modules <- function(axis_data, top_n = 15) {
  
  # Identify top hub cell states by median degree
  hubs <- axis_data %>%
    purrr::map_dfr("node_metrics") %>%
    group_by(cell_state) %>%
    summarise(med_deg = median(degree), .groups = "drop") %>%
    arrange(desc(med_deg)) %>%
    slice_head(n = top_n) %>%
    pull(cell_state)
  
  if (length(hubs) < 2) {
    return(tibble(hub_state = hubs, module = 1))
  }
  
  # Calculate co-membership counts
  co <- matrix(0, length(hubs), length(hubs),
               dimnames = list(hubs, hubs))
  
  for (x in axis_data) {
    tb <- tibble(
      cell = V(x$graph)$name,
      com  = V(x$graph)$community
    ) %>%
      filter(cell %in% hubs)
    
    same <- inner_join(tb, tb, by = "com", suffix = c("_A", "_B")) %>%
      filter(cell_A < cell_B)
    
    for (i in seq_len(nrow(same))) {
      a <- same$cell_A[i]
      b <- same$cell_B[i]
      co[a, b] <- co[a, b] + 1
    }
  }
  
  co <- co + t(co)
  diag(co) <- 0
  
  # Module detection via Jaccard-weighted Louvain
  if (all(co == 0)) {
    mod_vec <- rep(1, length(hubs))
    names(mod_vec) <- hubs
  } else {
    rs <- rowSums(co)
    jacc <- co / (outer(rs, rs, "+") - co + 1e-9)
    g <- igraph::graph_from_adjacency_matrix(jacc, "undirected",
                                             weighted = TRUE, diag = FALSE)
    mod_vec <- igraph::cluster_louvain(g, weights = E(g)$weight)$membership
    names(mod_vec) <- hubs
  }
  
  tibble(
    hub_state = names(mod_vec),
    module    = unname(mod_vec)
  ) %>%
    arrange(module, hub_state)
}


# --- Chord diagram visualization ---

draw_hub_chord <- function(tbl, title,
                           start_deg = 130,
                           label_cex = 1.2,
                           label_rad = 1.4,
                           title_cex = 1.5) {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(oma = c(2, 2, 3, 2))
  
  if (!"cell_state" %in% names(tbl)) names(tbl)[1L] <- "cell_state"
  if (!"module" %in% names(tbl))     names(tbl)[2L] <- "module"
  
  # Shorten cell state names
  tbl <- tbl %>%
    mutate(cell_state = str_replace(cell_state,
                                    "Monocytes\\.and\\.Macrophages",
                                    "Macrophages"))
  
  # Create edge list (Community -> cell state)
  edges <- tbl %>%
    transmute(
      from  = paste0("C", module),
      to    = cell_state,
      value = 1
    )
  
  # Color setup
  mod_ids <- sort(unique(edges$from))
  n_mod <- length(mod_ids)
  
  if (n_mod <= 3) {
    mod_col <- setNames(c("#E41A1C", "#377EB8", "#4DAF4A")[seq_len(n_mod)],
                        mod_ids)
  } else {
    mod_col <- setNames(
      brewer.pal(max(3, n_mod), "Set2")[seq_len(n_mod)],
      mod_ids
    )
  }
  
  cs_ids <- setdiff(unique(edges$to), mod_ids)
  grid_col <- c(mod_col, setNames(rep("#CCCCCC", length(cs_ids)), cs_ids))
  
  sectors <- c(mod_ids, cs_ids)
  gap.after <- c(rep(3, length(sectors) - 1), 10)
  
  circos.clear()
  circos.par(
    start.degree    = start_deg,
    gap.after       = gap.after,
    track.margin    = c(0.02, 0.02),
    canvas.xlim     = c(-1.3, 1.3),
    canvas.ylim     = c(-1.3, 1.3)
  )
  
  chordDiagram(
    edges,
    order               = sectors,
    grid.col            = grid_col,
    transparency        = 0.2,
    directional         = 0,
    annotationTrack     = "grid",
    annotationTrackHeight = 0.04
  )
  
  # Add labels
  circos.trackPlotRegion(
    track.index = 1,
    bg.border   = NA,
    panel.fun   = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      xlim <- get.cell.meta.data("xlim")
      font_weight <- ifelse(grepl("^C\\d+$", sector.index), 2, 1)
      circos.text(mean(xlim), label_rad, sector.index,
                  facing = "clockwise", niceFacing = TRUE,
                  cex = label_cex, font = font_weight,
                  adj = c(0.5, 0))
    }
  )
  
  mtext(title, side = 3, line = 0.5, cex = title_cex, font = 2)
}


# =============================================================================
# SECTION 6: FUNCTIONAL ENRICHMENT ANALYSIS
# =============================================================================

gene_tbl <- readRDS("../data/example_data/all_marker_gene.rds")

make_gene_lists <- function(hub_tbl, network_tag) {
  gene_tbl %>%
    inner_join(hub_tbl, by = c(cell_state = "hub_state")) %>%
    group_by(module) %>%
    summarise(gene_vec = list(unique(Gene)), .groups = "drop") %>%
    mutate(module_tag = paste0(network_tag, "_M", module)) %>%
    select(module_tag, gene_vec)
}

run_go_enrichment <- function(hub_modules_list) {
  gene_lists <- bind_rows(
    make_gene_lists(hub_modules_list$CXCL13_CXCR5, "CXCL13_CXCR5"),
    make_gene_lists(hub_modules_list$CXCL12_CXCR4, "CXCL12_CXCR4"),
    make_gene_lists(hub_modules_list$CCL19_21_CCR7, "CCL19_21_CCR7")
  )
  
  go_results <- pmap_dfr(gene_lists, function(module_tag, gene_vec) {
    ego <- enrichGO(
      gene          = gene_vec,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
    if (is.null(ego) || nrow(ego) == 0) return(tibble())
    as_tibble(ego) %>% mutate(module = module_tag)
  })
  
  return(go_results)
}

plot_go_enrichment <- function(go_results, topn = 10) {
  
  ordered_levels <- c(
    "CXCL13_CXCR5_C1", "CXCL12_CXCR4_C1", "CCL19_21_CCR7_C1",
    "CXCL13_CXCR5_C2", "CXCL12_CXCR4_C2", "CCL19_21_CCR7_C2"
  )
  
  plot_df <- go_results %>%
    mutate(module = str_replace(module, "_M(\\d+)$", "_C\\1")) %>%
    filter(!grepl("adaptive immune response based on somatic recombination",
                  Description, ignore.case = TRUE)) %>%
    group_by(module) %>%
    slice_min(p.adjust, n = topn, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      module         = factor(module, levels = rev(ordered_levels)),
      Description    = str_wrap(Description, width = 40),
      community_type = ifelse(grepl("_C1$", module), "C1", "C2")
    )
  
  p <- ggplot(plot_df,
              aes(x = Description, y = module,
                  size = Count, colour = -log10(p.adjust))) +
    geom_point() +
    scale_colour_viridis_c(name = "-log10 FDR") +
    scale_size(range = c(3, 10), name = "Gene count") +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      x     = "GO BP term",
      y     = "Hub community",
      title = "Functional enrichment of hub communities"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1,
                                  size = 14, margin = margin(t = 5)),
      axis.text.y  = element_text(
        size = 13,
        colour = ifelse(grepl("_C1$", rev(ordered_levels)),
                        "#E41A1C", "#377EB8")
      ),
      axis.title.x  = element_text(size = 15, face = "bold"),
      axis.title.y  = element_text(size = 15, face = "bold"),
      plot.title     = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.margin    = margin(10, 10, 30, 10),
      panel.grid.major.x = element_line(linewidth = 0.3),
      panel.background    = element_rect(fill = "white", colour = NA)
    )
  
  return(p)
}

# =============================================================================
# SECTION 7: DIRECTIONAL INTERACTION HEATMAP
# =============================================================================

create_significant_interactions_heatmap <- function(axis_data,
                                                    axis_name,
                                                    score_threshold_quantile = 0.8,
                                                    min_samples = 2,
                                                    top_n_pairs = 50) {
  
  cat("=== Creating Significant Interactions Heatmap:", axis_name, "===\n")
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
    
    for (i in rownames(s)) {
      for (j in colnames(s)) {
        pair_name <- paste(i, j, sep = "->")
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
  cat("Score threshold (", score_threshold_quantile * 100,
      "th percentile):", round(score_threshold, 3), "\n")
  
  high_score <- mean_mat >= score_threshold
  high_score[is.na(high_score)] <- FALSE
  
  # Statistical significance test
  significant_mat <- matrix(FALSE, n, n,
                            dimnames = list(cell_types, cell_types))
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
    sender      = rep(cell_types, each = n),
    receiver    = rep(cell_types, times = n),
    score       = as.vector(significant_scores),
    n_samples   = as.vector(cnt_mat),
    significant = as.vector(final_filter),
    stringsAsFactors = FALSE
  ) %>%
    filter(significant & score > 0) %>%
    arrange(desc(score)) %>%
    head(top_n_pairs)
  
  cat("Top", nrow(interaction_df), "significant interactions selected\n")
  
  # Create filtered matrix
  involved_cells <- unique(c(interaction_df$sender, interaction_df$receiver))
  filtered_mat <- matrix(0, length(involved_cells), length(involved_cells),
                         dimnames = list(involved_cells, involved_cells))
  for (k in 1:nrow(interaction_df)) {
    filtered_mat[interaction_df$sender[k],
                 interaction_df$receiver[k]] <- interaction_df$score[k]
  }
  filtered_mat[filtered_mat == 0] <- NA
  
  row_has_data <- apply(filtered_mat, 1, function(x) any(!is.na(x) & x > 0))
  col_has_data <- apply(filtered_mat, 2, function(x) any(!is.na(x) & x > 0))
  if (sum(row_has_data) == 0 || sum(col_has_data) == 0) {
    stop("No valid data remaining. Try relaxing thresholds.")
  }
  filtered_mat <- filtered_mat[row_has_data, col_has_data, drop = FALSE]
  
  # Create annotations
  sender_scores <- rowSums(filtered_mat, na.rm = TRUE)
  row_annotation <- data.frame(Total_Sent = sender_scores,
                               row.names = names(sender_scores))
  receiver_scores <- colSums(filtered_mat, na.rm = TRUE)
  col_annotation <- data.frame(Total_Received = receiver_scores,
                               row.names = names(receiver_scores))
  annotation_colors <- list(
    Total_Sent     = colorRampPalette(c("white", "coral", "red"))(100),
    Total_Received = colorRampPalette(c("white", "lightblue", "darkblue"))(100)
  )
  
  can_cluster_rows <- nrow(filtered_mat) > 1 &&
    sum(!is.na(filtered_mat)) > nrow(filtered_mat)
  can_cluster_cols <- ncol(filtered_mat) > 1 &&
    sum(!is.na(filtered_mat)) > ncol(filtered_mat)
  
  # Axis labels: S = sender, R = receiver
  lab_row <- paste0("S: ", rownames(filtered_mat))
  lab_col <- paste0("R: ", colnames(filtered_mat))
  dimnames(filtered_mat) <- list(lab_row, lab_col)
  
  p <- tryCatch({
    pheatmap(
      filtered_mat,
      color = colorRampPalette(
        c("#f7fbff", "#6baed6", "#2171b5", "#08306b"))(100),
      na_col = "grey95",
      cluster_rows    = can_cluster_rows,
      cluster_cols    = can_cluster_cols,
      annotation_row  = if (can_cluster_rows) row_annotation else NULL,
      annotation_col  = if (can_cluster_cols) col_annotation else NULL,
      annotation_colors = annotation_colors,
      legend       = FALSE,
      main         = paste0("Significant ", axis_name, " Interactions\n",
                            "(p<0.05, >", score_threshold_quantile * 100,
                            "th percentile)\nS: sender / R: receiver"),
      fontsize_row = 8, fontsize_col = 8,
      angle_col    = 45, border_color = "white",
      cellwidth    = 15, cellheight = 15,
      silent       = TRUE
    )
  }, error = function(e) {
    cat("Clustering failed, creating unclustered heatmap...\n")
    pheatmap(
      filtered_mat,
      color = colorRampPalette(
        c("#f7fbff", "#6baed6", "#2171b5", "#08306b"))(100),
      na_col = "grey95",
      cluster_rows = FALSE, cluster_cols = FALSE,
      legend       = FALSE,
      main         = paste0("Significant ", axis_name, " Interactions\n",
                            "(p<0.05, >", score_threshold_quantile * 100,
                            "th percentile)\nS: sender / R: receiver"),
      fontsize_row = 8, fontsize_col = 8,
      angle_col    = 45, border_color = "white",
      cellwidth    = 15, cellheight = 15,
      silent       = TRUE
    )
  })
  
  cat("\nTop 10 interactions:\n")
  print(head(interaction_df[, c("sender", "receiver", "score", "n_samples")], 10))
  
  return(list(
    heatmap        = p,
    interaction_df = interaction_df,
    filtered_matrix = filtered_mat,
    involved_cells  = rownames(filtered_mat)
  ))
}

# =============================================================================
# SECTION 8: HUB LINEAGE DIVERSITY ANALYSIS
# =============================================================================

analyze_hub_lineage_diversity <- function(axis_data, hub_list, axis_name) {
  
  cat("\n=== Analyzing Lineage Diversity for", axis_name, "===\n")
  
  lineage_mapping <- c(
    "B"                          = "B_cell",
    "CD4.T"                      = "T_cell",
    "CD8.T"                      = "T_cell",
    "Dendritic"                  = "Dendritic_cell",
    "Mast"                       = "Mast_cell",
    "PMNs"                       = "Neutrophil",
    "PCs"                        = "Plasma_cell",
    "NK"                         = "NK_cell",
    "Monocytes.and.Macrophages"  = "Macrophage",
    "Endothelial"                = "Endothelial",
    "Epithelial"                 = "Epithelial",
    "Fibroblasts"                = "Fibroblast"
  )
  
  get_lineage <- function(cell_state) {
    for (prefix in names(lineage_mapping)) {
      if (grepl(paste0("^", prefix), cell_state)) {
        return(lineage_mapping[prefix])
      }
    }
    return("Other")
  }
  
  hub_diversity_results <- list()
  
  for (hub in hub_list) {
    cat("\n--- Analyzing", hub, "---\n")
    
    all_partners <- c()
    interaction_counts <- c()
    
    for (sample_id in names(axis_data)) {
      score_mat <- axis_data[[sample_id]]$score
      if (!hub %in% rownames(score_mat) && !hub %in% colnames(score_mat)) next
      
      if (hub %in% rownames(score_mat)) {
        sender_scores <- score_mat[hub, ]
        sender_scores <- sender_scores[sender_scores > 0]
        if (length(sender_scores) > 0) {
          all_partners <- c(all_partners, names(sender_scores))
          interaction_counts <- c(interaction_counts, sender_scores)
        }
      }
      if (hub %in% colnames(score_mat)) {
        receiver_scores <- score_mat[, hub]
        receiver_scores <- receiver_scores[receiver_scores > 0]
        if (length(receiver_scores) > 0) {
          all_partners <- c(all_partners, names(receiver_scores))
          interaction_counts <- c(interaction_counts, receiver_scores)
        }
      }
    }
    
    if (length(all_partners) == 0) {
      cat("  No interaction partners found\n")
      next
    }
    
    partner_lineages <- sapply(all_partners, get_lineage)
    hub_lineage <- get_lineage(hub)
    
    n_cross <- sum(partner_lineages != hub_lineage)
    pct_cross <- round(100 * n_cross / length(partner_lineages), 1)
    
    cat("  Hub lineage:", hub_lineage, "\n")
    cat("  Total interactions:", length(all_partners), "\n")
    cat("  Cross-lineage:", n_cross, "(", pct_cross, "%)\n")
    
    hub_diversity_results[[hub]] <- list(
      hub                  = hub,
      hub_lineage          = hub_lineage,
      total_interactions   = length(all_partners),
      unique_partners      = length(unique(all_partners)),
      unique_lineages      = length(unique(partner_lineages)),
      n_intra_lineage      = sum(partner_lineages == hub_lineage),
      n_cross_lineage      = n_cross,
      pct_cross_lineage    = pct_cross,
      lineage_distribution = as.data.frame(sort(table(partner_lineages),
                                                decreasing = TRUE)),
      all_partners         = all_partners,
      partner_lineages     = partner_lineages
    )
  }
  
  summary_df <- map_dfr(hub_diversity_results, function(x) {
    data.frame(
      hub                = x$hub,
      hub_lineage        = x$hub_lineage,
      total_interactions = x$total_interactions,
      unique_partners    = x$unique_partners,
      unique_lineages    = x$unique_lineages,
      pct_cross_lineage  = x$pct_cross_lineage
    )
  })
  
  cat("\n=== Summary Table ===\n")
  print(summary_df)
  
  return(list(
    detailed_results = hub_diversity_results,
    summary_table    = summary_df
  ))
}


# --- Lineage diversity visualization ---

plot_lineage_diversity <- function(combined_diversity) {
  
  axis_levels <- c("CXCL13-CXCR5", "CXCL12-CXCR4", "CCL19/21-CCR7")
  
  lineage_colors <- c(
    "Neutrophil"     = "#4DBBD5",
    "B_cell"         = "#E64B35",
    "Plasma_cell"    = "#00A087",
    "T_cell"         = "#F4A460",
    "Dendritic_cell" = "#3C5488",
    "Endothelial"    = "#F39B7F",
    "Epithelial"     = "#8491B4",
    "Macrophage"     = "#91D1C2",
    "NK_cell"        = "#DC0000"
  )
  
  make_plot <- function(data_subset, community_name) {
    hub_order <- data_subset %>%
      group_by(hub, community) %>%
      summarise(mean_cross = mean(pct_cross_lineage), .groups = "drop") %>%
      arrange(desc(mean_cross)) %>%
      pull(hub)
    
    data_subset <- data_subset %>%
      mutate(hub_label = gsub("Monocytes.and.Macrophages",
                              "Macrophages", hub))
    hub_label_order <- gsub("Monocytes.and.Macrophages",
                            "Macrophages", hub_order)
    
    hub_lineage_df <- data_subset %>%
      select(hub_label, hub_lineage) %>%
      distinct() %>%
      mutate(hub_label = factor(hub_label, levels = hub_label_order))
    
    x_colors <- lineage_colors[
      hub_lineage_df$hub_lineage[match(hub_label_order,
                                       hub_lineage_df$hub_label)]
    ]
    
    plot_df <- data_subset %>%
      mutate(
        pct_intra = 100 - pct_cross_lineage,
        hub_label = factor(hub_label, levels = hub_label_order),
        axis      = factor(axis, levels = axis_levels)
      ) %>%
      pivot_longer(
        cols      = c(pct_cross_lineage, pct_intra),
        names_to  = "interaction_type",
        values_to = "pct"
      ) %>%
      mutate(
        interaction_type = factor(
          interaction_type,
          levels = c("pct_intra", "pct_cross_lineage"),
          labels = c("Intra-lineage", "Cross-lineage")
        ),
        x_id = as.numeric(hub_label)
      )
    
    label_df <- data_subset %>%
      mutate(
        hub_label = factor(hub_label, levels = hub_label_order),
        axis      = factor(axis, levels = axis_levels),
        label     = paste0(pct_cross_lineage, "%\n(",
                           unique_lineages, " lineages)"),
        x_id      = as.numeric(hub_label)
      )
    
    p <- ggplot(plot_df, aes(x = x_id, y = pct, fill = interaction_type)) +
      geom_bar(stat = "identity", width = 0.65,
               color = "white", linewidth = 0.4) +
      geom_text(
        data = label_df,
        aes(x = x_id, y = 50, label = label),
        inherit.aes = FALSE,
        size = 3.0, fontface = "bold", color = "white", lineheight = 0.9
      ) +
      geom_hline(yintercept = 50, linetype = "dashed",
                 color = "gray40", linewidth = 0.5, alpha = 0.7) +
      facet_wrap(~ axis, ncol = 1, strip.position = "right") +
      scale_fill_manual(
        values = c("Cross-lineage" = "#2171B5", "Intra-lineage" = "#D9D9D9"),
        name   = "Interaction Type"
      ) +
      scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
      scale_x_continuous(
        breaks = seq_along(hub_label_order),
        labels = hub_label_order,
        expand = c(0.05, 0.05)
      ) +
      labs(
        title    = paste("Community", community_name,
                         "- Cross-lineage Interaction Diversity"),
        subtitle = paste0("Bars show proportion of intra- vs. ",
                          "cross-lineage interactions;\nannotations show ",
                          "cross-lineage % and number of distinct ",
                          "partner lineages"),
        x = NULL,
        y = "% of Total Interactions"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title       = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(size = 11, hjust = 0.5,
                                        color = "gray40", lineheight = 1.2),
        axis.text.x      = element_text(size = 11, color = x_colors,
                                        face = "bold"),
        axis.text.y      = element_text(size = 11),
        axis.title.y     = element_text(size = 12, face = "bold"),
        strip.text       = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "gray95", color = NA),
        legend.position  = "bottom",
        legend.title     = element_text(size = 11, face = "bold"),
        legend.text      = element_text(size = 11),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin = margin(15, 15, 10, 15)
      )
    
    return(p)
  }
  
  p_c1 <- make_plot(combined_diversity %>% filter(community == "C1"), "C1")
  p_c2 <- make_plot(combined_diversity %>% filter(community == "C2"), "C2")
  
  return(list(C1 = p_c1, C2 = p_c2))
}


# =============================================================================
# SECTION 9: RUN ALL ANALYSES
# =============================================================================

# --- 9.1 Hubness analysis ---
hubness_results <- list(
  CXCL13_CXCR5  = plot_hubness_analysis(LR_analysis$CXCL13_CXCR5,
                                        "CXCL13-CXCR5"),
  CXCL12_CXCR4  = plot_hubness_analysis(LR_analysis$CXCL12_CXCR4,
                                        "CXCL12-CXCR4"),
  CCL19_21_CCR7 = plot_hubness_analysis(LR_analysis$CCL19_21_CCR7,
                                        "CCL19/21-CCR7")
)

ggsave(hubness_results$CXCL13_CXCR5$plot,
       file = "../output/03_CXCL13-CXCR5_Network_Prominence.pdf",
       width = 10, height = 8, dpi = 300)
ggsave(hubness_results$CXCL12_CXCR4$plot,
       file = "../output/03_CXCL12_CXCR4_Network_Prominence.pdf",
       width = 10, height = 8, dpi = 300)
ggsave(hubness_results$CCL19_21_CCR7$plot,
       file = "../output/03_CCL19_21_CCR7_Network_Prominence.pdf",
       width = 10, height = 8, dpi = 300)

# --- 9.2 Hub community identification ---
hub_modules <- list(
  CXCL13_CXCR5  = get_hub_modules(LR_analysis$CXCL13_CXCR5, top_n = 15),
  CXCL12_CXCR4  = get_hub_modules(LR_analysis$CXCL12_CXCR4, top_n = 15),
  CCL19_21_CCR7 = get_hub_modules(LR_analysis$CCL19_21_CCR7, top_n = 15)
)

# --- 9.3 Chord diagrams ---
pdf(file = "../output/03_CXCL13-CXCR5_Hub_community.pdf",
    width = 10, height = 8)
suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CXCL13_CXCR5,
                 title = "CXCL13-CXCR5: Hub Community to Cell State")
))
dev.off()

pdf(file = "../output/03_CXCL12_CXCR4_Hub_community.pdf",
    width = 10, height = 8)
suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CXCL12_CXCR4,
                 title = "CXCL12-CXCR4: Hub Community to Cell State")
))
dev.off()

pdf(file = "../output/03_CCL19_21_CCR7_Hub_community.pdf",
    width = 10, height = 8)
suppressWarnings(suppressMessages(
  draw_hub_chord(hub_modules$CCL19_21_CCR7,
                 title = "CCL19/21-CCR7: Hub Community to Cell State")
))
dev.off()

# --- 9.4 GO enrichment ---
go_results <- run_go_enrichment(hub_modules)
go_plot <- plot_go_enrichment(go_results, topn = 10)

ggsave(go_plot,
       file   = "../output/03_GO_enrichment_hub_communities.pdf",
       width  = 14,
       height = 10,
       dpi    = 300)

# --- 9.5 Directional interaction heatmaps ---
interaction_results <- list(
  CXCL13_CXCR5 = create_significant_interactions_heatmap(
    LR_analysis$CXCL13_CXCR5, "CXCL13-CXCR5"
  ),
  CXCL12_CXCR4 = create_significant_interactions_heatmap(
    LR_analysis$CXCL12_CXCR4, "CXCL12-CXCR4"
  ),
  CCL19_21_CCR7 = create_significant_interactions_heatmap(
    LR_analysis$CCL19_21_CCR7, "CCL19/21-CCR7"
  )
)

# --- 9.6 Lineage diversity ---
c1_hubs <- c("B_S01", "CD4.T_S02", "Dendritic_S01", "PMNs_S03", "PCs_S01")
c2_hubs <- c("Endothelial_S02", "Epithelial_S04",
             "Monocytes.and.Macrophages_S01", "NK_S01", "NK_S03")
all_hubs <- c(c1_hubs, c2_hubs)

diversity_CXCL13 <- analyze_hub_lineage_diversity(
  LR_analysis$CXCL13_CXCR5, all_hubs, "CXCL13-CXCR5")
diversity_CXCL12 <- analyze_hub_lineage_diversity(
  LR_analysis$CXCL12_CXCR4, all_hubs, "CXCL12-CXCR4")
diversity_CCL19 <- analyze_hub_lineage_diversity(
  LR_analysis$CCL19_21_CCR7, all_hubs, "CCL19/21-CCR7")

combined_diversity <- bind_rows(
  diversity_CXCL13$summary_table %>% mutate(axis = "CXCL13-CXCR5"),
  diversity_CXCL12$summary_table %>% mutate(axis = "CXCL12-CXCR4"),
  diversity_CCL19$summary_table  %>% mutate(axis = "CCL19/21-CCR7")
) %>%
  mutate(community = case_when(
    hub %in% c1_hubs ~ "C1",
    hub %in% c2_hubs ~ "C2",
    TRUE ~ "Other"
  ))

cat("\n=== Overall Summary Across All Axes ===\n")
cat("Median cross-lineage %:", median(combined_diversity$pct_cross_lineage), "\n")
cat("Median unique lineages:", median(combined_diversity$unique_lineages), "\n")

plots <- plot_lineage_diversity(combined_diversity)

ggsave(plots$C1,
       file = "../output/03_Hub_Lineage_Diversity_C1.pdf",
       width = 10, height = 8, dpi = 300)
ggsave(plots$C2,
       file = "../output/03_Hub_Lineage_Diversity_C2.pdf",
       width = 10, height = 8, dpi = 300)

# --- 9.7 C2 members TLS vs NonTLS ---
c2_members <- c("Endothelial_S02", "Epithelial_S04",
                "Monocytes.and.Macrophages_S01", "NK_S01", "NK_S03")

c2_abundance_list <- lapply(names(combined_obj), function(sid) {
  meta <- combined_obj[[sid]]$cellstate_norm
  available_states <- intersect(c2_members, colnames(meta))
  if (length(available_states) == 0 || !"Label" %in% colnames(meta)) return(NULL)
  
  meta %>%
    select(ID, Label, all_of(available_states)) %>%
    pivot_longer(cols = all_of(available_states),
                 names_to  = "cell_state",
                 values_to = "abundance") %>%
    mutate(sample = sid)
})

c2_abundance_df <- bind_rows(c2_abundance_list) %>%
  filter(!is.na(Label), !is.na(abundance))

# Remove outliers beyond 99th percentile per cell state (for visualization)
c2_plot_df <- c2_abundance_df %>%
  group_by(cell_state) %>%
  mutate(upper_99 = quantile(abundance, 0.99, na.rm = TRUE)) %>%
  filter(abundance <= upper_99) %>%
  ungroup()

# Calculate fold change and Wilcoxon p-value per cell state
c2_stats <- c2_abundance_df %>%
  group_by(cell_state) %>%
  summarise(
    fc = mean(abundance[Label == "TLS"], na.rm = TRUE) /
      (mean(abundance[Label == "NonTLS"], na.rm = TRUE) + 1e-10),
    p_val = tryCatch(
      wilcox.test(abundance[Label == "TLS"],
                  abundance[Label == "NonTLS"],
                  exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    fc_label = paste0("FC: ", round(fc, 2)),
    sig_label = case_when(
      p_val < 0.0001 ~ "****",
      p_val < 0.001  ~ "***",
      p_val < 0.01   ~ "**",
      p_val < 0.05   ~ "*",
      TRUE           ~ "ns"
    ),
    cell_state_label = gsub("Monocytes\\.and\\.Macrophages",
                            "Macrophages", cell_state)
  )

c2_plot_df <- c2_plot_df %>%
  left_join(c2_stats %>% select(cell_state, cell_state_label,
                                fc_label, sig_label, p_val),
            by = "cell_state") %>%
  mutate(
    cell_state_label = factor(
      cell_state_label,
      levels = unique(gsub("Monocytes\\.and\\.Macrophages",
                           "Macrophages", c2_members))
    ),
    Label = factor(Label, levels = c("NonTLS", "TLS"))
  )

# Significance label position
c2_sig_pos <- c2_plot_df %>%
  group_by(cell_state_label) %>%
  summarise(
    sig_y     = max(abundance, na.rm = TRUE) * 1.08,
    fc_y      = max(abundance, na.rm = TRUE) * 1.00,
    sig_label = unique(sig_label),
    fc_label  = unique(fc_label),
    .groups   = "drop"
  )

p_c2 <- ggplot(c2_plot_df,
               aes(x = Label, y = abundance, fill = Label)) +
  geom_violin(trim = TRUE, alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.25, outlier.shape = NA, alpha = 0.8) +
  geom_text(data = c2_sig_pos,
            aes(x = 1.5, y = sig_y, label = sig_label),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  geom_text(data = c2_sig_pos,
            aes(x = 1.5, y = fc_y, label = fc_label),
            inherit.aes = FALSE, size = 3.8,
            color = "#003399", fontface = "bold") +
  facet_wrap(~ cell_state_label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("NonTLS" = "#AEC6E8", "TLS" = "#F4B8B0"),
                    name = "Label") +
  labs(
    title = "C2 Hub Community Members: TLS vs NonTLS",
    y     = "Abundance",
    x     = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title       = element_text(size = 15, face = "bold", hjust = 0),
    strip.text       = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", color = "black"),
    axis.text.x      = element_text(size = 11),
    axis.text.y      = element_text(size = 10),
    axis.title.y     = element_text(size = 12, face = "bold"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(p_c2,
       file   = "../output/03_C2_TLS_vs_NonTLS.pdf",
       width  = 10,
       height = 8,
       dpi    = 300)

# =============================================================================
# END OF SCRIPT
# =============================================================================