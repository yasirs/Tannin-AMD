# Phase 2-4: GSE29801-PRG4 GSEA Cross-Validation
# Comprehensive enrichment analysis of GSE29801 gene sets in PRG4 data

#%%
library(dplyr)
library(fgsea)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyr)

#%%
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
OUTPUT_DIR <- file.path(BASE_DIR, "results/external-validation/gse29801_prg4_comparison")
GENE_SET_DIR <- file.path(OUTPUT_DIR, "gene_sets")
GSEA_DIR <- file.path(OUTPUT_DIR, "gsea_results")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")

dir.create(GSEA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, showWarnings = FALSE)

FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

#%%
# ========================================
# PHASE 2: PREPARE ALL GENE SETS
# ========================================

message("=== Phase 2: Preparing Gene Sets ===\n")

# Load GSE29801 results
macular_de <- read.csv(file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de/macular_variance_analysis_annotated.csv"))
extramacular_de <- read.csv(file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de/extramacular_variance_analysis_annotated.csv"))
macular_age <- read.csv(file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de/macular_age_age_main.csv"))
extramacular_age <- read.csv(file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de/extramacular_age_age_main.csv"))

# Load age-adjusted variance gene sets
macular_low_adj <- read.csv(file.path(GENE_SET_DIR, "macular_low_variance_age_adj.csv"))
macular_high_adj <- read.csv(file.path(GENE_SET_DIR, "macular_high_variance_age_adj.csv"))
macular_bimod_adj <- read.csv(file.path(GENE_SET_DIR, "macular_bimodality_age_adjusted.csv"))
extramacular_low_adj <- read.csv(file.path(GENE_SET_DIR, "extramacular_low_variance_age_adj.csv"))
extramacular_high_adj <- read.csv(file.path(GENE_SET_DIR, "extramacular_high_variance_age_adj.csv"))
extramacular_bimod_adj <- read.csv(file.path(GENE_SET_DIR, "extramacular_bimodality_age_adjusted.csv"))

# Define all gene sets
define_all_gene_sets <- function(tissue_name, de_results, age_results, low_var, high_var, bimod) {
  
  message(sprintf("Tissue: %s", tissue_name))
  
  # 1. DE genes (p < 0.01) - use mean_control/mean_amd columns to get logFC
  de_genes <- de_results %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    mutate(logFC = log2((mean_amd + 0.01) / (mean_control + 0.01))) %>%
    mutate(f_pvalue = ifelse(is.na(f_pvalue), 1, f_pvalue)) %>%
    filter(f_pvalue < 0.01) %>%
    pull(gene_symbol) %>%
    unique()
  
  # 2. Low variance (age-adjusted, FDR < 0.05)
  low_var_genes <- low_var %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    pull(gene_symbol) %>%
    unique()
  
  # 3. High variance (age-adjusted, FDR < 0.05)
  high_var_genes <- high_var %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    pull(gene_symbol) %>%
    unique()
  
  # 4. Bimodal (age-adjusted)
  bimod_genes <- bimod %>%
    filter(is_bimodal == TRUE & !is.na(gene_symbol) & gene_symbol != "") %>%
    pull(gene_symbol) %>%
    unique()
  
  # 5. Age DE (FDR < 0.05) - need to merge with DE results to get gene symbols
  age_with_symbols <- age_results %>%
    left_join(de_results %>% select(gene, gene_symbol), by = "gene")
  
  age_genes <- age_with_symbols %>%
    filter(adj.P.Val < 0.05 & !is.na(gene_symbol) & gene_symbol != "") %>%
    pull(gene_symbol) %>%
    unique()
  
  message(sprintf("  DE (p<0.01): %d genes", length(de_genes)))
  message(sprintf("  Low Var (age-adj): %d genes", length(low_var_genes)))
  message(sprintf("  High Var (age-adj): %d genes", length(high_var_genes)))
  message(sprintf("  Bimodal (age-adj): %d genes", length(bimod_genes)))
  message(sprintf("  Age DE: %d genes", length(age_genes)))
  
  return(list(
    de = de_genes,
    low_var = low_var_genes,
    high_var = high_var_genes,
    bimodal = bimod_genes,
    age_de = age_genes
  ))
}

macular_all_sets <- define_all_gene_sets(
  "Macular", macular_de, macular_age,
  macular_low_adj, macular_high_adj, macular_bimod_adj
)

extramacular_all_sets <- define_all_gene_sets(
  "Extramacular", extramacular_de, extramacular_age,
  extramacular_low_adj, extramacular_high_adj, extramacular_bimod_adj
)

# Save gene set summary
gene_set_summary <- data.frame(
  Tissue = rep(c("Macular", "Extramacular"), each = 5),
  GeneSet = rep(c("DE (p<0.01)", "Low Variance", "High Variance", "Bimodal", "Age DE"), 2),
  NumGenes = c(
    length(macular_all_sets$de), length(macular_all_sets$low_var),
    length(macular_all_sets$high_var), length(macular_all_sets$bimodal),
    length(macular_all_sets$age_de),
    length(extramacular_all_sets$de), length(extramacular_all_sets$low_var),
    length(extramacular_all_sets$high_var), length(extramacular_all_sets$bimodal),
    length(extramacular_all_sets$age_de)
  )
)

write.csv(gene_set_summary,
          file.path(TABLE_DIR, "gene_set_sizes.csv"),
          row.names = FALSE)

message(sprintf("\n✓ Gene set summary saved"))

#%%
# ========================================
# PHASE 3: LOAD PRG4 BULK DATA
# ========================================

message("\n=== Phase 3: Loading PRG4 Bulk Data ===\n")

# Load DE results
prg4_data <- read.csv(file.path(BASE_DIR, "data/RPE_DE_results.csv"))

message(sprintf("Loaded PRG4 data: %d genes", nrow(prg4_data)))

# Create ranked gene lists for GSEA
# 1. H2O2 vs Control
h2o2_ranked <- prg4_data %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "" & !is.na(H2O2_vs_CTRL.log2FoldChange)) %>%
  arrange(desc(H2O2_vs_CTRL.log2FoldChange)) %>%
  select(hgnc_symbol, H2O2_vs_CTRL.log2FoldChange)

h2o2_ranked <- h2o2_ranked[!duplicated(h2o2_ranked$hgnc_symbol), ]
h2o2_ranks <- setNames(h2o2_ranked$H2O2_vs_CTRL.log2FoldChange, h2o2_ranked$hgnc_symbol)

# 2. PRG4+H2O2 vs H2O2 (rescue)
prg4_rescue_ranked <- prg4_data %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "" & !is.na(H2O2PRG4_vs_H2O2.log2FoldChange)) %>%
  arrange(desc(H2O2PRG4_vs_H2O2.log2FoldChange)) %>%
  select(hgnc_symbol, H2O2PRG4_vs_H2O2.log2FoldChange)

prg4_rescue_ranked <- prg4_rescue_ranked[!duplicated(prg4_rescue_ranked$hgnc_symbol), ]
prg4_rescue_ranks <- setNames(prg4_rescue_ranked$H2O2PRG4_vs_H2O2.log2FoldChange, prg4_rescue_ranked$hgnc_symbol)

message(sprintf("H2O2 vs Control: %d genes", length(h2o2_ranks)))
message(sprintf("PRG4 Rescue: %d genes", length(prg4_rescue_ranks)))

#%%
# ========================================
# PHASE 4: GSEA ENRICHMENT ANALYSIS
# ========================================

message("\n=== Phase 4: GSEA Enrichment Analysis ===\n")

# Function to run GSEA
run_gsea_gse29801 <- function(gene_set_list, ranked_genes, comparison_name, tissue_name) {
  
  message(sprintf("\n%s - %s:", tissue_name, comparison_name))
  
  results_list <- list()
  
  for (set_name in names(gene_set_list)) {
    gene_set <- gene_set_list[[set_name]]
    
    if (length(gene_set) < 5) {
      message(sprintf("  %s: Skipped (too few genes: %d)", set_name, length(gene_set)))
      next
    }
    
    # Calculate overlap
    overlap <- sum(gene_set %in% names(ranked_genes))
    overlap_pct <- 100 * overlap / length(gene_set)
    
    if (overlap < 5 || overlap_pct < 50) {
      message(sprintf("  %s: Skipped (low overlap: %d/%d = %.1f%%)", 
                      set_name, overlap, length(gene_set), overlap_pct))
      next
    }
    
    # Run GSEA
    pathway_list <- list()
    pathway_list[[set_name]] <- gene_set
    
    gsea_result <- fgsea(
      pathways = pathway_list,
      stats = ranked_genes,
      minSize = 5,
      nproc = 1
    )
    
    gsea_result$gene_set <- set_name
    gsea_result$tissue <- tissue_name
    gsea_result$comparison <- comparison_name
    gsea_result$overlap_genes <- overlap
    gsea_result$total_genes <- length(gene_set)
    gsea_result$overlap_pct <- overlap_pct
    
    results_list[[set_name]] <- gsea_result
    
    message(sprintf("  %s: NES=%.2f, p=%.3f, overlap=%d/%d (%.1f%%)",
                    set_name, gsea_result$NES, gsea_result$pval,
                    overlap, length(gene_set), overlap_pct))
  }
  
  if (length(results_list) > 0) {
    return(bind_rows(results_list))
  } else {
    return(NULL)
  }
}

# Run GSEA for all combinations
all_gsea_results <- list()

# Macular - H2O2
all_gsea_results[["macular_h2o2"]] <- run_gsea_gse29801(
  macular_all_sets, h2o2_ranks, "H2O2_vs_Control", "Macular"
)

# Macular - PRG4 Rescue
all_gsea_results[["macular_prg4"]] <- run_gsea_gse29801(
  macular_all_sets, prg4_rescue_ranks, "PRG4_Rescue", "Macular"
)

# Extramacular - H2O2
all_gsea_results[["extramacular_h2o2"]] <- run_gsea_gse29801(
  extramacular_all_sets, h2o2_ranks, "H2O2_vs_Control", "Extramacular"
)

# Extramacular - PRG4 Rescue
all_gsea_results[["extramacular_prg4"]] <- run_gsea_gse29801(
  extramacular_all_sets, prg4_rescue_ranks, "PRG4_Rescue", "Extramacular"
)

# Combine all results
gsea_combined <- bind_rows(all_gsea_results)

# Convert leadingEdge list to string
gsea_combined <- gsea_combined %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(head(x, 10), collapse=", ")))

# Save complete GSEA results
write.csv(gsea_combined,
          file.path(TABLE_DIR, "gsea_all_results.csv"),
          row.names = FALSE)

message(sprintf("\n✓ GSEA complete: %d enrichment tests", nrow(gsea_combined)))

#%%
# Create summary table
gsea_summary <- gsea_combined %>%
  select(tissue, gene_set, comparison, NES, pval, padj, overlap_genes, total_genes, overlap_pct) %>%
  arrange(tissue, comparison, pval)

write.csv(gsea_summary,
          file.path(TABLE_DIR, "gsea_summary.csv"),
          row.names = FALSE)

print(gsea_summary)

message("\n=== Phases 2-4 Complete ===")
