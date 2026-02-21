#%%
# GSE135092 Bimodality Analysis
# replicate methodology from GSE29801: Test high-variance genes for bimodal expression in AMD samples
# using mclust (BIC difference > 10)

library(dplyr)
library(readr)
library(mclust)
library(tidyr)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")
bimod_dir <- file.path(results_dir, "bimodality")
dir.create(bimod_dir, showWarnings = FALSE)

# Load Data
cat("Loading data...\n")
# Expression and Metadata
cache_file <- file.path(base_dir, "results/cohort-GSE135092/expression_matrix_cache.rds")
if (!file.exists(cache_file)) stop("Cache file not found!")
data_cache <- readRDS(cache_file)

# CORRECT SLOT NAMES
expr_mat <- data_cache$expr_mat # Log2 CPM
meta <- data_cache$meta         # Metadata

# Variance Results (Macula & non-Macula)
res_mac <- read_csv(file.path(results_dir, "RPE_Macula_all_results_symbols.csv"), show_col_types=FALSE)
res_nonmac <- read_csv(file.path(results_dir, "RPE_nonMacula_all_results_symbols.csv"), show_col_types=FALSE)

cat("Data loaded.\n")

#%%
# Function: Test Bimodality
test_bimodality_gene <- function(gene_expr_amd) {
  if (length(gene_expr_amd) < 10) return(list(is_bimodal=FALSE, bic_diff=NA))
  
  tryCatch({
    fit1 <- Mclust(gene_expr_amd, G=1, verbose=FALSE)
    fit2 <- Mclust(gene_expr_amd, G=2, verbose=FALSE)
    
    if (is.null(fit1) || is.null(fit2)) return(list(is_bimodal=FALSE, bic_diff=NA))
    
    bic_diff <- fit2$bic - fit1$bic
    is_bimodal <- bic_diff > 10 # Strong evidence
    
    return(list(is_bimodal=is_bimodal, bic_diff=bic_diff))
  }, error = function(e) {
    return(list(is_bimodal=FALSE, bic_diff=NA))
  })
}

# Wrapper for list of genes
run_bimodality_analysis <- function(res_df, expr_data, metadata, tissue_type, tissue_string_match, top_n=100) {
  cat(sprintf("\nAnalyzing %s (Tag: %s)...\n", tissue_type, tissue_string_match))
  
  # Filter metadata using 'tissue' column directly
  # Because 'region' derived column has a bug (Macula grep matches non-Macula)
  
  samples_keep <- metadata %>% 
    filter(grepl(tissue_string_match, tissue, fixed=TRUE)) %>%
    pull(gsm)
  
  cat(sprintf("  Total samples: %d\n", length(samples_keep)))
  
  # Separate by Disease Status (amd_status column)
  meta_sub <- metadata %>% filter(gsm %in% samples_keep)
  
  amd_samples <- meta_sub %>% 
    filter(grepl("AMD", amd_status, ignore.case=TRUE)) %>% 
    pull(gsm)
  
  cat(sprintf("  AMD Samples: %d\n", length(amd_samples)))
  
  if (length(amd_samples) < 5) {
    cat("  Not enough AMD samples.\n")
    return(NULL)
  }
  
  # Select Top Variance Genes (Increased Variance Only?)
  # Filter for valid symbols
  top_genes <- res_df %>%
    filter(levene_pval < 0.01) %>%
    filter(!is.na(gene_symbol)) %>%
    arrange(levene_pval) %>%
    head(top_n)
  
  cat(sprintf("  Testing top %d variance genes...\n", nrow(top_genes)))
  
  results_list <- list()
  
  for (i in 1:nrow(top_genes)) {
    sym <- top_genes$gene_symbol[i]
    ens <- top_genes$gene[i]
    
    # Check if in expr matrix
    if (ens %in% rownames(expr_data)) {
      gene_id <- ens
    } else if (sym %in% rownames(expr_data)) {
      gene_id <- sym
    } else {
      next
    }
    
    # Get expression
    # Use intersect just in case sample names mismatch slightly
    valid_samples <- intersect(amd_samples, colnames(expr_data))
    if (length(valid_samples) < 5) next
    
    expr_vals <- expr_data[gene_id, valid_samples]
    
    # Test
    res <- test_bimodality_gene(expr_vals)
    
    if (!is.na(res$bic_diff)) {
      results_list[[length(results_list)+1]] <- data.frame(
        gene_symbol = sym,
        gene_id = ens,
        levene_pval = top_genes$levene_pval[i],
        var_ratio = top_genes$var_ratio[i],
        is_bimodal = res$is_bimodal,
        bic_diff = res$bic_diff
      )
    }
  }
  
  if (length(results_list) > 0) {
    out_df <- do.call(rbind, results_list)
    return(out_df)
  } else {
    return(NULL)
  }
}

#%%
# Run Analysis
# Use exact strings or unique substrings for tissue column
# "RPE, Macula"
# "RPE, non-Macula"

bimod_mac <- run_bimodality_analysis(res_mac, expr_mat, meta, "Macula", "RPE, Macula", top_n=200)
bimod_nonmac <- run_bimodality_analysis(res_nonmac, expr_mat, meta, "Non-Macula", "RPE, non-Macula", top_n=200)

#%%
# Save Results
if (!is.null(bimod_mac)) {
  write_csv(bimod_mac, file.path(bimod_dir, "macular_bimodality.csv"))
  n_bi <- sum(bimod_mac$is_bimodal)
  cat(sprintf("  Macula: %d bimodal genes identified (out of %d tested)\n", n_bi, nrow(bimod_mac)))
}

if (!is.null(bimod_nonmac)) {
  write_csv(bimod_nonmac, file.path(bimod_dir, "nonmacular_bimodality.csv"))
  n_bi <- sum(bimod_nonmac$is_bimodal)
  cat(sprintf("  Non-Macula: %d bimodal genes identified (out of %d tested)\n", n_bi, nrow(bimod_nonmac)))
}

cat("\nBimodality analysis complete.\n")
