#!/usr/bin/env Rscript
# fgsea-based GSEA for Perturb-seq Knockdowns
#
# Computes GSEA enrichment scores, NES, and p-values for all knockdowns
# against all gene sets using the fast fgsea algorithm.
#
# Usage: Rscript run_fgsea_perturbseq.R [--dataset RPE1|K562_GWPS] [--nperm 10000]

library(fgsea)
library(data.table)
library(rhdf5)
library(Matrix)

#%%
# Configuration
PROJECT_ROOT <- "/home/ysuhail/work/Tannin-AMD"
RESULTS_DIR <- file.path(PROJECT_ROOT, "results/perturbseq-analysis")
CACHE_DIR <- file.path(RESULTS_DIR, "cache")

# Data paths
PERTURBSEQ_DATA <- list(
  RPE1 = file.path(PROJECT_ROOT, "data/external/perturbseq/rpe1_normalized_bulk_01.h5ad"),
  K562_GWPS = file.path(PROJECT_ROOT, "data/external/perturbseq/K562_gwps_normalized_bulk_01.h5ad")
)

#%%
# Load gene sets
load_gene_sets <- function() {
  gene_sets <- list()
  
  # Age gene sets
  age_up <- fread(file.path(RESULTS_DIR, "gene-sets/Age-UP.csv"))$gene
  age_down <- fread(file.path(RESULTS_DIR, "gene-sets/Age-DOWN.csv"))$gene
  are <- fread(file.path(RESULTS_DIR, "gene-sets/ARE.csv"))$gene
  
  gene_sets[["Age-UP"]] <- age_up
  gene_sets[["Age-DOWN"]] <- age_down
  gene_sets[["ARE"]] <- are
  
  # AMD gene sets
  amd_macula_up <- fread(file.path(RESULTS_DIR, "gene-sets/AMD-Macula-UP.csv"))$gene
  amd_macula_down <- fread(file.path(RESULTS_DIR, "gene-sets/AMD-Macula-DOWN.csv"))$gene
  amd_nonmacula_up <- fread(file.path(RESULTS_DIR, "gene-sets/AMD-nonMacula-UP.csv"))$gene
  amd_nonmacula_down <- fread(file.path(RESULTS_DIR, "gene-sets/AMD-nonMacula-DOWN.csv"))$gene
  
  gene_sets[["AMD-Macula-UP"]] <- amd_macula_up
  gene_sets[["AMD-Macula-DOWN"]] <- amd_macula_down
  gene_sets[["AMD-nonMacula-UP"]] <- amd_nonmacula_up
  gene_sets[["AMD-nonMacula-DOWN"]] <- amd_nonmacula_down
  
  # Axis gene sets
  axis_dir <- file.path(PROJECT_ROOT, "results/three-axes/gene-sets")
  gene_sets[["Senescence-PRO"]] <- fread(file.path(axis_dir, "senescence_pro_disease.csv"))$gene
  gene_sets[["Senescence-ANTI"]] <- fread(file.path(axis_dir, "senescence_anti_disease.csv"))$gene
  gene_sets[["Redox-PRO"]] <- fread(file.path(axis_dir, "oxidative_pro_disease.csv"))$gene
  gene_sets[["Redox-ANTI"]] <- fread(file.path(axis_dir, "oxidative_anti_disease.csv"))$gene
  gene_sets[["Inflammation-PRO"]] <- fread(file.path(axis_dir, "inflammatory_pro_disease.csv"))$gene
  gene_sets[["Inflammation-ANTI"]] <- fread(file.path(axis_dir, "inflammatory_anti_disease.csv"))$gene
  
  cat("Loaded", length(gene_sets), "gene sets\n")
  for (name in names(gene_sets)) {
    cat("  ", name, ":", length(gene_sets[[name]]), "genes\n")
  }
  
  return(gene_sets)
}

#%%
# Load Perturb-seq data from cached pickle (via Python export)
# We'll load the pre-processed CSVs instead for simplicity
load_perturbseq_expression <- function(dataset_name) {
  cat("Loading expression data for", dataset_name, "...\n")
  
  # We need to export from Python first - check if exported
  expr_file <- file.path(CACHE_DIR, paste0(dataset_name, "_expression.csv.gz"))
  meta_file <- file.path(CACHE_DIR, paste0(dataset_name, "_metadata.csv"))
  
  if (!file.exists(expr_file)) {
    stop(paste("Expression file not found:", expr_file, 
               "\nRun: python export_for_fgsea.py --dataset", dataset_name))
  }
  
  # Load expression matrix (perturbations x genes)
  expr_dt <- fread(expr_file)
  gene_names <- colnames(expr_dt)[-1]  # First col is perturbation ID
  perturbation_ids <- expr_dt[[1]]
  expr_matrix <- as.matrix(expr_dt[, -1, with=FALSE])
  rownames(expr_matrix) <- perturbation_ids
  colnames(expr_matrix) <- gene_names
  
  # Load metadata
  meta <- fread(meta_file)
  
  cat("  Loaded:", nrow(expr_matrix), "perturbations x", ncol(expr_matrix), "genes\n")
  
  return(list(
    expression = expr_matrix,
    gene_names = gene_names,
    metadata = meta
  ))
}

#%%
# Run fgsea for a single knockdown
run_fgsea_single <- function(ranked_list, gene_sets, nperm = 10000) {
  result <- fgsea(
    pathways = gene_sets,
    stats = ranked_list,
    minSize = 5,
    maxSize = 3000,
    nperm = nperm,
    nproc = 1
  )
  return(result)
}

#%%
# Run fgsea for all knockdowns
run_fgsea_all <- function(expr_matrix, gene_names, metadata, gene_sets, nperm = 10000) {
  n_perts <- nrow(expr_matrix)
  
  cat("Running fgsea on", n_perts, "perturbations...\n")
  
  # Initialize results
  all_results <- list()
  
  # Progress tracking
  start_time <- Sys.time()
  
  for (i in seq_len(n_perts)) {
    if (i %% 100 == 0) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      rate <- i / as.numeric(elapsed)
      remaining <- (n_perts - i) / rate
      cat(sprintf("\r  Progress: %d/%d (%.1f min elapsed, ~%.1f min remaining)", 
                  i, n_perts, elapsed, remaining))
    }
    
    # Skip non-targeting controls
    target_gene <- metadata$target_gene[i]
    if (grepl("non-targeting|nontargeting|control", tolower(target_gene))) {
      next
    }
    
    # Create ranked list (gene names -> expression values)
    expr_vector <- expr_matrix[i, ]
    ranked_list <- setNames(as.numeric(expr_vector), gene_names)
    
    # Remove NA/Inf values
    ranked_list <- ranked_list[is.finite(ranked_list)]
    
    if (length(ranked_list) < 100) {
      next  # Skip if too few valid genes
    }
    
    ranked_list <- sort(ranked_list, decreasing = TRUE)
    
    # Run fgsea
    result <- tryCatch({
      suppressWarnings(run_fgsea_single(ranked_list, gene_sets, nperm))
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(result)) next
    
    # Add perturbation info
    result$perturbation_id <- metadata$perturbation_id[i]
    result$target_gene <- target_gene
    
    all_results[[i]] <- result
  }
  
  cat("\n  Done!\n")
  
  # Combine results
  combined <- rbindlist(all_results, fill = TRUE)
  
  return(combined)
}

#%%
# Main function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Parse arguments
  dataset <- "RPE1"
  nperm <- 10000
  
  for (i in seq_along(args)) {
    if (args[i] == "--dataset" && i < length(args)) {
      dataset <- args[i + 1]
    }
    if (args[i] == "--nperm" && i < length(args)) {
      nperm <- as.integer(args[i + 1])
    }
  }
  
  cat("=== fgsea Perturb-seq Analysis ===\n")
  cat("Dataset:", dataset, "\n")
  cat("Permutations:", nperm, "\n\n")
  
  # Load gene sets
  gene_sets <- load_gene_sets()
  
  # Load expression data
  data <- load_perturbseq_expression(dataset)
  
  # Run fgsea
  results <- run_fgsea_all(
    data$expression, 
    data$gene_names, 
    data$metadata, 
    gene_sets, 
    nperm
  )
  
  # Save results
  output_file <- file.path(RESULTS_DIR, "fgsea", paste0("fgsea_", dataset, "_results.csv"))
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  fwrite(results, output_file)
  
  cat("\nResults saved to:", output_file, "\n")
  cat("Total rows:", nrow(results), "\n")
  
  # Summary
  cat("\n=== Summary ===\n")
  summary_dt <- results[, .(
    n_significant = sum(padj < 0.05, na.rm = TRUE),
    mean_NES = mean(NES, na.rm = TRUE),
    median_pval = median(pval, na.rm = TRUE)
  ), by = pathway]
  print(summary_dt)
}

#%%
# Run if called directly
if (!interactive()) {
  main()
}
