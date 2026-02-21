#%%
# GSE135092 Directional Variance Analysis - CORRECTED (Using ORA)
# Replacing incorrec GSEA on filtered sets with proper ORA

library(dplyr)
library(readr)
library(msigdbr)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")

cat("GSE135092 Directional Variance Analysis - ORA Correction\n")
cat(rep("=", 70), "\n\n", sep = "")

#%%
# Load results
results_mac <- read_csv(file.path(results_dir, "RPE_Macula_all_results_symbols.csv"), 
                        show_col_types = FALSE)
results_nonmac <- read_csv(file.path(results_dir, "RPE_nonMacula_all_results_symbols.csv"),
                           show_col_types = FALSE)

#%%
# Download MSigDB gene sets
cat("1. Downloading MSigDB gene sets...\n")
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)

msig_gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gobp_sets <- split(msig_gobp$gene_symbol, msig_gobp$gs_name)

msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
c2_sets <- split(msig_c2$gene_symbol, msig_c2$gs_name)

cat(sprintf("  Hallmark: %d sets\n", length(hallmark_sets)))
cat(sprintf("  GO BP: %d sets\n", length(gobp_sets)))
cat(sprintf("  C2 Full: %d sets\n\n", length(c2_sets)))

#%%
# Function to perform ORA (hypergeometric test)
perform_ora <- function(gene_set, pathway_sets, background_genes, set_name, db_name) {
  
  cat(sprintf("  Running ORA: %s - %s\n", set_name, db_name))
  
  # Remove genes not in background
  gene_set <- intersect(gene_set, background_genes)
  
  if (length(gene_set) == 0) {
    cat("    No genes in set!\n")
    return(data.frame())
  }
  
  cat(sprintf("    Gene set size: %d\n", length(gene_set)))
  cat(sprintf("    Background size: %d\n", length(background_genes)))
  
  # Initialize results
  results <- data.frame(
    pathway = character(),
    pathway_size = integer(),
    overlap = integer(),
    pval = numeric(),
    odds_ratio = numeric(),
    gene_list = character(),
    stringsAsFactors = FALSE
  )
  
  # Test each pathway
  for (pathway_name in names(pathway_sets)) {
    pathway_genes <- pathway_sets[[pathway_name]]
    pathway_genes <- intersect(pathway_genes, background_genes)
    
    if (length(pathway_genes) < 15 || length(pathway_genes) > 500) next
    
    # Overlap
    overlap_genes <- intersect(gene_set, pathway_genes)
    overlap_size <- length(overlap_genes)
    
    if (overlap_size == 0) next
    
    # Hypergeometric test
    # White balls (pathway genes), black balls (non-pathway), 
    # Draws (gene set size), white balls drawn (overlap)
    p_val <- phyper(
      q = overlap_size - 1,
      m = length(pathway_genes),
      n = length(background_genes) - length(pathway_genes),
      k = length(gene_set),
      lower.tail = FALSE
    )
    
    # Odds ratio
    # (overlap / (gene_set - overlap)) / ((pathway - overlap) / (background - gene_set - pathway + overlap))
    a <- overlap_size
    b <- length(gene_set) - overlap_size
    c <- length(pathway_genes) - overlap_size
    d <- length(background_genes) - length(gene_set) - length(pathway_genes) + overlap_size
    
    or <- (a * d) / (b * c + 0.001)
    
    results <- rbind(results, data.frame(
      pathway = pathway_name,
      pathway_size = length(pathway_genes),
      overlap = overlap_size,
      pval = p_val,
      odds_ratio = or,
      gene_list = paste(head(overlap_genes, 10), collapse = ", "),
      stringsAsFactors = FALSE
    ))
  }
  
  # FDR correction
  if (nrow(results) > 0) {
    results$fdr <- p.adjust(results$pval, method = "BH")
    results <- results %>% arrange(pval)
    
    cat(sprintf("    Pathways tested: %d\n", nrow(results)))
    cat(sprintf("    Significant (FDR<0.25): %d\n", sum(results$fdr < 0.25)))
  }
  
  return(results)
}

#%%
# Extract gene sets
cat("2. Extracting directional variance gene sets...\n")

# Background: all genes with symbols
background_mac <- results_mac %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

background_nonmac <- results_nonmac %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

# Increased variance genes (p<0.01)
inc_var_mac <- results_mac %>%
  filter(!is.na(var_ratio), var_ratio > 1, 
         !is.na(levene_pval), levene_pval < 0.01,
         !is.na(gene_symbol), gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

inc_var_nonmac <- results_nonmac %>%
  filter(!is.na(var_ratio), var_ratio > 1,
         !is.na(levene_pval), levene_pval < 0.01,
         !is.na(gene_symbol), gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

# Decreased variance genes (p<0.01)
dec_var_mac <- results_mac %>%
  filter(!is.na(var_ratio), var_ratio < 1,
         !is.na(levene_pval), levene_pval < 0.01,
         !is.na(gene_symbol), gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

dec_var_nonmac <- results_nonmac %>%
  filter(!is.na(var_ratio), var_ratio < 1,
         !is.na(levene_pval), levene_pval < 0.01,
         !is.na(gene_symbol), gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

cat(sprintf("  Background genes:\n"))
cat(sprintf("    Macula: %d\n", length(background_mac)))
cat(sprintf("    non-Macula: %d\n", length(background_nonmac)))

cat(sprintf("\n  Increased variance (p<0.01):\n"))
cat(sprintf("    Macula: %d\n", length(inc_var_mac)))
cat(sprintf("    non-Macula: %d\n", length(inc_var_nonmac)))

cat(sprintf("\n  Decreased variance (p<0.01):\n"))
cat(sprintf("    Macula: %d\n", length(dec_var_mac)))
cat(sprintf("    non-Macula: %d\n\n", length(dec_var_nonmac)))

#%%
# Run ORA for Increased Variance
cat("3. ORA for Increased Variance Genes\n")
cat(rep("-", 70), "\n", sep = "")

ora_inc_hall_mac <- perform_ora(inc_var_mac, hallmark_sets, background_mac,
                                "Macula Increased", "Hallmark")
ora_inc_hall_nonmac <- perform_ora(inc_var_nonmac, hallmark_sets, background_nonmac,
                                   "non-Macula Increased", "Hallmark")

ora_inc_gobp_mac <- perform_ora(inc_var_mac, gobp_sets, background_mac,
                                "Macula Increased", "GO BP")
ora_inc_gobp_nonmac <- perform_ora(inc_var_nonmac, gobp_sets, background_nonmac,
                                   "non-Macula Increased", "GO BP")

ora_inc_c2_mac <- perform_ora(inc_var_mac, c2_sets, background_mac,
                              "Macula Increased", "C2")
ora_inc_c2_nonmac <- perform_ora(inc_var_nonmac, c2_sets, background_nonmac,
                                 "non-Macula Increased", "C2")

#%%
# Run ORA for Decreased Variance
cat("\n4. ORA for Decreased Variance Genes\n")
cat(rep("-", 70), "\n", sep = "")

ora_dec_hall_mac <- perform_ora(dec_var_mac, hallmark_sets, background_mac,
                                "Macula Decreased", "Hallmark")
ora_dec_hall_nonmac <- perform_ora(dec_var_nonmac, hallmark_sets, background_nonmac,
                                   "non-Macula Decreased", "Hallmark")

ora_dec_gobp_mac <- perform_ora(dec_var_mac, gobp_sets, background_mac,
                                "Macula Decreased", "GO BP")
ora_dec_gobp_nonmac <- perform_ora(dec_var_nonmac, gobp_sets, background_nonmac,
                                   "non-Macula Decreased", "GO BP")

ora_dec_c2_mac <- perform_ora(dec_var_mac, c2_sets, background_mac,
                              "Macula Decreased", "C2")
ora_dec_c2_nonmac <- perform_ora(dec_var_nonmac, c2_sets, background_nonmac,
                                 "non-Macula Decreased", "C2")

#%%
# Save results
cat("\n5. Saving ORA results...\n")

write_csv(ora_inc_hall_mac, file.path(results_dir, "gsea/ORA_increased_Hallmark_macula.csv"))
write_csv(ora_inc_hall_nonmac, file.path(results_dir, "gsea/ORA_increased_Hallmark_nonmacula.csv"))
write_csv(ora_inc_gobp_mac, file.path(results_dir, "gsea/ORA_increased_GOBP_macula.csv"))
write_csv(ora_inc_gobp_nonmac, file.path(results_dir, "gsea/ORA_increased_GOBP_nonmacula.csv"))
write_csv(ora_inc_c2_mac, file.path(results_dir, "gsea/ORA_increased_C2_macula.csv"))
write_csv(ora_inc_c2_nonmac, file.path(results_dir, "gsea/ORA_increased_C2_nonmacula.csv"))

write_csv(ora_dec_hall_mac, file.path(results_dir, "gsea/ORA_decreased_Hallmark_macula.csv"))
write_csv(ora_dec_hall_nonmac, file.path(results_dir, "gsea/ORA_decreased_Hallmark_nonmacula.csv"))
write_csv(ora_dec_gobp_mac, file.path(results_dir, "gsea/ORA_decreased_GOBP_macula.csv"))
write_csv(ora_dec_gobp_nonmac, file.path(results_dir, "gsea/ORA_decreased_GOBP_nonmacula.csv"))
write_csv(ora_dec_c2_mac, file.path(results_dir, "gsea/ORA_decreased_C2_macula.csv"))
write_csv(ora_dec_c2_nonmac, file.path(results_dir, "gsea/ORA_decreased_C2_nonmacula.csv"))

#%%
# Create summary table
cat("\n6. Creating summary...\n")

summary_ora <- data.frame(
  Tissue = rep(c("RPE Macula", "RPE non-Macula"), each = 2),
  Direction = rep(c("Increased", "Decreased"), 2),
  Gene_Count = c(
    length(inc_var_mac), length(dec_var_mac),
    length(inc_var_nonmac), length(dec_var_nonmac)
  ),
  Hallmark_FDR025 = c(
    sum(ora_inc_hall_mac$fdr < 0.25, na.rm = TRUE),
    sum(ora_dec_hall_mac$fdr < 0.25, na.rm = TRUE),
    sum(ora_inc_hall_nonmac$fdr < 0.25, na.rm = TRUE),
    sum(ora_dec_hall_nonmac$fdr < 0.25, na.rm = TRUE)
  ),
  GOBP_FDR025 = c(
    sum(ora_inc_gobp_mac$fdr < 0.25, na.rm = TRUE),
    sum(ora_dec_gobp_mac$fdr < 0.25, na.rm = TRUE),
    sum(ora_inc_gobp_nonmac$fdr < 0.25, na.rm = TRUE),
    sum(ora_dec_gobp_nonmac$fdr < 0.25, na.rm = TRUE)
  ),
  C2_FDR025 = c(
    sum(ora_inc_c2_mac$fdr < 0.25, na.rm = TRUE),
    sum(ora_dec_c2_mac$fdr < 0.25, na.rm = TRUE),
    sum(ora_inc_c2_nonmac$fdr < 0.25, na.rm = TRUE),
    sum(ora_dec_c2_nonmac$fdr < 0.25, na.rm = TRUE)
  )
)

write_csv(summary_ora, file.path(results_dir, "ORA_directional_summary.csv"))

cat("\n", rep("=", 70), "\n", sep = "")
cat("ORA Analysis Complete!\n")
cat(rep("=", 70), "\n\n", sep = "")

print(summary_ora)

cat("\nFiles saved to:", file.path(results_dir, "gsea/"), "\n")
