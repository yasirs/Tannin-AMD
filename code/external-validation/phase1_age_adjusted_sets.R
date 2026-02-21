# Phase 1: Complete Age-Adjusted Variance Gene Sets
# Define LOW, HIGH, and BIMODAL variance genes using age-adjusted data

#%%
library(dplyr)
library(mclust)
library(pheatmap)
library(RColorBrewer)

#%%
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- file.path(BASE_DIR, "results/external-validation/gse29801_prg4_comparison")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "gene_sets"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)

FONT_FAMILY <- "Palatino Linotype"

#%%
# Load data
message("=== Loading Data ===")
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

# Load age-adjusted variance results
macular_var_adj <- read.csv(file.path(DATA_DIR, "macular_variance_age_adjusted.csv"))
extramacular_var_adj <- read.csv(file.path(DATA_DIR, "extramacular_variance_age_adjusted.csv"))

# Load original variance (for var_fc)
macular_var <- read.csv(file.path(DATA_DIR, "macular_variance_analysis_annotated.csv"))
extramacular_var <- read.csv(file.path(DATA_DIR, "extramacular_variance_analysis_annotated.csv"))

# Merge to get var_fc and gene symbols
macular_var_adj <- macular_var_adj %>%
  left_join(macular_var %>% select(gene, gene_symbol, var_fc, var_control, var_amd),
            by = "gene")

extramacular_var_adj <- extramacular_var_adj %>%
  left_join(extramacular_var %>% select(gene, gene_symbol, var_fc, var_control, var_amd),
            by = "gene")

#%%
# ========================================
# DEFINE AGE-ADJUSTED GENE SETS
# ========================================

message("\n=== Defining Age-Adjusted Gene Sets ===")

define_variance_sets <- function(var_adj_results, tissue_name) {
  
  message(sprintf("\nTissue: %s", tissue_name))
  
  # LOW variance (decreased in AMD)
  low_var <- var_adj_results %>%
    filter(f_fdr_adj < 0.05 & var_fc < 0) %>%
    arrange(var_fc)
  
  # HIGH variance (increased in AMD)
  high_var <- var_adj_results %>%
    filter(f_fdr_adj < 0.05 & var_fc > 0) %>%
    arrange(desc(var_fc))
  
  message(sprintf("  Low variance genes: %d", nrow(low_var)))
  message(sprintf("  High variance genes: %d", nrow(high_var)))
  
  return(list(
    low = low_var,
    high = high_var
  ))
}

macular_sets <- define_variance_sets(macular_var_adj, "Macular")
extramacular_sets <- define_variance_sets(extramacular_var_adj, "Extramacular")

#%%
# ========================================
# BIMODALITY TESTING ON AGE-ADJUSTED DATA
# ========================================

message("\n=== Bimodality Testing (Age-Adjusted) ===")

test_bimodality_age_adjusted <- function(expr_data, metadata, var_results, tissue_name, top_n = 100) {
  
  message(sprintf("\nTesting %s tissue...", tissue_name))
  
  # Get high-variance genes (age-adjusted)
  high_var_genes <- var_results %>%
    filter(f_fdr_adj < 0.05 & var_fc > 0) %>%
    arrange(f_pvalue_adj) %>%
    head(top_n)
  
  if (nrow(high_var_genes) == 0) {
    message("  No high-variance genes to test")
    return(NULL)
  }
  
  message(sprintf("  Testing %d genes", nrow(high_var_genes)))
  
  # Get AMD sample indices
  amd_idx <- which(metadata$disease_status == "AMD")
  
  # Test each gene
  bimodality_results <- data.frame()
  
  for (i in 1:nrow(high_var_genes)) {
    gene_id <- as.character(high_var_genes$gene[i])
    gene_symbol <- high_var_genes$gene_symbol[i]
    
    if (gene_id %in% colnames(expr_data)) {
      gene_expr <- expr_data[, gene_id]
      
      # Get AMD expression
      amd_expr <- gene_expr[amd_idx]
      amd_expr <- amd_expr[!is.na(amd_expr)]
      
      # Remove age effect
      age_clean <- metadata$age[amd_idx][!is.na(amd_expr)]
      amd_expr_clean <- amd_expr[!is.na(age_clean)]
      age_clean <- age_clean[!is.na(age_clean)]
      
      if (length(amd_expr_clean) < 10) next
      
      # Regress out age
      age_model <- lm(amd_expr_clean ~ age_clean)
      residuals <- residuals(age_model)
      
      # Test bimodality on residuals
      tryCatch({
        fit1 <- Mclust(residuals, G = 1, verbose = FALSE)
        fit2 <- Mclust(residuals, G = 2, verbose = FALSE)
        
        if (!is.null(fit1) && !is.null(fit2)) {
          bic_diff <- fit2$bic - fit1$bic
          is_bimodal <- bic_diff > 10
          
          bimodality_results <- rbind(bimodality_results, data.frame(
            gene = gene_id,
            gene_symbol = ifelse(is.na(gene_symbol), "", gene_symbol),
            var_fc = high_var_genes$var_fc[i],
            f_fdr_adj = high_var_genes$f_fdr_adj[i],
            is_bimodal = is_bimodal,
            bic_improvement = bic_diff
          ))
        }
      }, error = function(e) {})
    }
  }
  
  n_bimodal <- sum(bimodality_results$is_bimodal, na.rm = TRUE)
  message(sprintf("  Bimodal genes (BIC > 10): %d / %d", n_bimodal, nrow(bimodality_results)))
  
  return(bimodality_results)
}

# Run bimodality testing
macular_bimod_adj <- test_bimodality_age_adjusted(
  macular_data$expr, 
  macular_data$metadata,
  macular_var_adj,
  "Macular"
)

extramacular_bimod_adj <- test_bimodality_age_adjusted(
  extramacular_data$expr,
  extramacular_data$metadata, 
  extramacular_var_adj,
  "Extramacular"
)

# Save bimodality results
write.csv(macular_bimod_adj,
          file.path(OUTPUT_DIR, "gene_sets", "macular_bimodality_age_adjusted.csv"),
          row.names = FALSE)
write.csv(extramacular_bimod_adj,
          file.path(OUTPUT_DIR, "gene_sets", "extramacular_bimodality_age_adjusted.csv"),
          row.names = FALSE)

#%%
# Save gene sets
message("\n=== Saving Gene Sets ===")

save_gene_set <- function(gene_df, filename) {
  write.csv(gene_df, 
            file.path(OUTPUT_DIR, "gene_sets", filename),
            row.names = FALSE)
  message(sprintf("  Saved: %s (%d genes)", filename, nrow(gene_df)))
}

# Macular
save_gene_set(macular_sets$low, "macular_low_variance_age_adj.csv")
save_gene_set(macular_sets$high, "macular_high_variance_age_adj.csv")

# Extramacular
save_gene_set(extramacular_sets$low, "extramacular_low_variance_age_adj.csv")
save_gene_set(extramacular_sets$high, "extramacular_high_variance_age_adj.csv")

message("\n=== Phase 1 Complete ===")
message(sprintf("Gene sets saved to: %s/gene_sets/", OUTPUT_DIR))
