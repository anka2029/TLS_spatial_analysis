############################################################
# Script 0: Prepare aggregated TLS / Non-TLS expression
#           profiles as EcoTyper input
#
# PURPOSE:
# This script generates pseudo-bulk gene expression profiles
# by averaging Visium spot-level expression within TLS and
# Non-TLS regions for each sample.
#
# The resulting matrices are used as INPUT for EcoTyper
# (bulk mode) to identify Carcinoma EcoTypes (CEs)
# representing shared transcriptional programs across samples.
#
# IMPORTANT:
# - This script is OPTIONAL.
# - It is intended for users who want to analyze
#   cross-sample, aggregated TLS vs Non-TLS characteristics.
# - If you already ran EcoTyper using Visium (spot-level) mode,
#   you can SKIP this script and start from Script 1.
#
############################################################


########################
# Required Input
########################

#   located in:
#     data/combined_obj.RData
#
# A named list "combined_obj", where each element corresponds
# to one Visium sample and contains:
#
# combined_obj[[sample_name]]$hs
#   - sparse gene expression matrix (genes x barcodes)
#
# combined_obj[[sample_name]]$cellstate_raw
#   - data.frame with at least:
#       ID    : Visium barcode
#       Label : "TLS" or "NonTLS"

# combined_obj[[sample_name]]$cellstate_norm (cell state abundance normalized)
#   - data.frame with at least:
#       ID    : Visium barcode
#       Label : "TLS" or "NonTLS"
#
# NOTE:
# TLS labels are assumed to be defined externally
# (e.g., based on pathology annotation of H&E images).
#
########################


########################
# Output
########################

# 1. gene_avg_by_sample.txt
#    - Gene expression matrix (genes x samples)
#    - Columns: sample_TLS / sample_NonTLS
#
# 2. sample_histology.txt
#    - Metadata file mapping each column to:
#      TLS or NonTLS
#
# These two files are used as direct input for EcoTyper bulk mode.
#
########################


########################
# Load packages
########################

library(Matrix)
library(dplyr)
library(purrr)
library(tidyr)

########################
# Load data
########################

# ---- Example data (default) ----
load("../data/example_data/combined_obj.RData")



########################
# Step 1: Define common gene set
########################

all_genes <- Reduce(
  intersect,
  lapply(combined_obj, function(x) rownames(x$hs))
)


########################
# Step 2: Initialize output matrix
########################

sample_names <- names(combined_obj)

col_names <- as.vector(rbind(
  paste0(sample_names, "_TLS"),
  paste0(sample_names, "_NonTLS")
))

out_mat <- matrix(
  NA_real_,
  nrow = length(all_genes),
  ncol = length(col_names),
  dimnames = list(all_genes, col_names)
)


########################
# Step 3: Compute TLS / NonTLS averages per sample
########################

for (s in sample_names) {
  
  message("Processing sample: ", s)
  
  expr <- combined_obj[[s]]$hs[all_genes, , drop = FALSE]
  meta <- combined_obj[[s]]$cellstate_norm
  
  # Extract barcodes
  tls_bcs <- meta %>%
    filter(Label == "TLS") %>%
    pull(ID) %>%
    intersect(colnames(expr))
  
  nontls_bcs <- meta %>%
    filter(Label != "TLS") %>%
    pull(ID) %>%
    intersect(colnames(expr))
  
  # Compute averages (sparse-aware)
  tls_avg <- if (length(tls_bcs) > 0) {
    Matrix::rowMeans(expr[, tls_bcs, drop = FALSE])
  } else {
    rep(NA_real_, length(all_genes))
  }
  
  nontls_avg <- if (length(nontls_bcs) > 0) {
    Matrix::rowMeans(expr[, nontls_bcs, drop = FALSE])
  } else {
    rep(NA_real_, length(all_genes))
  }
  
  out_mat[, paste0(s, "_TLS")]    <- tls_avg
  out_mat[, paste0(s, "_NonTLS")] <- nontls_avg
}


########################
# Step 4: Remove empty columns
########################

keep_cols <- colSums(!is.na(out_mat)) > 0
out_mat   <- out_mat[, keep_cols, drop = FALSE]


########################
# Step 5: Write EcoTyper input files
########################

expr_df <- data.frame(
  Gene = rownames(out_mat),
  out_mat,
  check.names = FALSE
)

meta_df <- data.frame(
  ID        = colnames(out_mat),
  Histology = ifelse(grepl("_TLS$", colnames(out_mat)), "TLS", "NonTLS"),
  stringsAsFactors = FALSE
)

write.table(
  expr_df,
  file = "../output/00_gene_avg_by_sample.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  meta_df,
  file = "../output/00_sample_histology.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("EcoTyper input files generated successfully.")
