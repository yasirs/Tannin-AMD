#%%
# Fix GSE135092 Gene Symbol Mapping
# Convert Ensembl IDs to Gene Symbols using biomaRt, then re-run analyses

library(biomaRt)
library(dplyr)
library(readr)
library(fgsea)
library(msigdbr)
library(ggplot2)
library(VennDiagram)
library(tidyr)
library(pheatmap)
library(grid)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")

cat("GSE135092 Gene Symbol Mapping Fix\n")
cat(rep("=", 70), "\n\n", sep = "")

#%%
# 1. Load results with Ensembl IDs
cat("1. Loading results with Ensembl IDs...\n")
results_mac <- read_csv(file.path(results_dir, "RPE_Macula_all_results.csv"), 
                        show_col_types = FALSE)
results_nonmac <- read_csv(file.path(results_dir, "RPE_nonMacula_all_results.csv"),
                           show_col_types = FALSE)

cat(sprintf("  RPE Macula: %d genes\n", nrow(results_mac)))
cat(sprintf("  RPE non-Macula: %d genes\n", nrow(results_nonmac)))
cat(sprintf("  Example gene ID: %s\n\n", results_mac$gene[1]))

#%%
# 2. Get unique Ensembl IDs
all_ensembl_ids <- unique(c(results_mac$gene, results_nonmac$gene))
cat(sprintf("2. Found %d unique Ensembl IDs\n\n", length(all_ensembl_ids)))

#%%
# 3. Map Ensembl IDs to Gene Symbols using biomaRt
cat("3. Mapping Ensembl IDs to gene symbols (this may take a few minutes)...\n")

# Connect to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Get mapping
gene_map <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
  filters = 'ensembl_gene_id',
  values = all_ensembl_ids,
  mart = ensembl
)

cat(sprintf("  Retrieved %d mappings\n", nrow(gene_map)))
cat(sprintf("  %.1f%% of genes have symbols\n", 
            100 * sum(gene_map$hgnc_symbol != "") / nrow(gene_map)))

# Clean up: use Ensembl ID if no symbol
gene_map <- gene_map %>%
  mutate(gene_symbol = ifelse(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol))

cat("\n")

#%%
# 4. Add gene symbols to results
cat("4. Adding gene symbols to results...\n")

results_mac_sym <- results_mac %>%
  left_join(gene_map, by = c("gene" = "ensembl_gene_id")) %>%
  mutate(gene_symbol = coalesce(gene_symbol, gene)) # Use Ensembl ID if mapping failed

results_nonmac_sym <- results_nonmac %>%
  left_join(gene_map, by = c("gene" = "ensembl_gene_id")) %>%
  mutate(gene_symbol = coalesce(gene_symbol, gene))

# Save updated results
write_csv(results_mac_sym, file.path(results_dir, "RPE_Macula_all_results_symbols.csv"))
write_csv(results_nonmac_sym, file.path(results_dir, "RPE_nonMacula_all_results_symbols.csv"))

cat("  Saved updated results with gene symbols\n\n")

#%%
# 5. Check for AMD genes
cat("5. Checking for known AMD genes...\n")

amd_genes <- c("CFH", "ARMS2", "C2", "C3", "CFB", "CFI", "HTRA1", 
               "TIMP3", "COL8A1", "BEST1", "RPE65", "ABCA4", 
               "EFEMP1", "APOE", "CETP", "LIPC", "VEGFA", "HMCN1")

amd_in_mac <- results_mac_sym %>%
  filter(gene_symbol %in% amd_genes) %>%
  select(gene_symbol, logFC, P.Value, adj.P.Val, levene_pval, levene_fdr, 
         var_ratio, variance_direction)

amd_in_nonmac <- results_nonmac_sym %>%
  filter(gene_symbol %in% amd_genes) %>%
  select(gene_symbol, logFC, P.Value, adj.P.Val, levene_pval, levene_fdr,
         var_ratio, variance_direction)

cat(sprintf("  Found %d AMD genes in Macula\n", nrow(amd_in_mac)))
cat(sprintf("  Found %d AMD genes in non-Macula\n", nrow(amd_in_nonmac)))

# Combine and save
amd_combined <- bind_rows(
  amd_in_mac %>% mutate(tissue = "RPE Macula"),
  amd_in_nonmac %>% mutate(tissue = "RPE non-Macula")
)

write_csv(amd_combined, file.path(results_dir, "AMD_genes_variance_FIXED.csv"))
cat("\n")

#%%
# 6. Re-run GSEA with gene symbols
cat("6. Re-running GSEA with gene symbols...\n")

# Download MSigDB
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)

msig_gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gobp_sets <- split(msig_gobp$gene_symbol, msig_gobp$gs_name)

msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP")
c2_sets <- split(msig_c2$gene_symbol, msig_c2$gs_name)

# Function to run GSEA
run_gsea_fixed <- function(results_df, gene_sets, tissue_name, db_name) {
  cat(sprintf("  Running GSEA: %s - %s\n", tissue_name, db_name))
  
  # Create ranked list using gene SYMBOLS
  ranked_list <- results_df %>%
    filter(!is.na(levene_pval), !is.na(var_ratio), gene_symbol != "") %>%
    mutate(rank_stat = -log10(levene_pval) * sign(var_ratio - 1)) %>%
    arrange(desc(rank_stat))
  
  ranks <- setNames(ranked_list$rank_stat, ranked_list$gene_symbol)
  
  # Run GSEA
  gsea_res <- fgsea(pathways = gene_sets,
                    stats = ranks,
                    minSize = 15,
                    maxSize = 500)
  
  gsea_res <- gsea_res %>%
    filter(padj < 0.25) %>%
    arrange(padj, desc(abs(NES)))
  
  cat(sprintf("    Found %d pathways (FDR<0.25)\n", nrow(gsea_res)))
  
  return(gsea_res)
}

# Run GSEA
gsea_h_mac <- run_gsea_fixed(results_mac_sym, hallmark_sets, "RPE Macula", "Hallmark")
gsea_h_nonmac <- run_gsea_fixed(results_nonmac_sym, hallmark_sets, "RPE non-Macula", "Hallmark")

gsea_gobp_mac <- run_gsea_fixed(results_mac_sym, gobp_sets, "RPE Macula", "GO BP")
gsea_gobp_nonmac <- run_gsea_fixed(results_nonmac_sym, gobp_sets, "RPE non-Macula", "GO BP")

gsea_c2_mac <- run_gsea_fixed(results_mac_sym, c2_sets, "RPE Macula", "C2")
gsea_c2_nonmac <- run_gsea_fixed(results_nonmac_sym, c2_sets, "RPE non-Macula", "C2")

# Save
write_csv(gsea_h_mac, file.path(results_dir, "gsea/variance_GSEA_Hallmark_macula_FIXED.csv"))
write_csv(gsea_h_nonmac, file.path(results_dir, "gsea/variance_GSEA_Hallmark_nonmacula_FIXED.csv"))
write_csv(gsea_gobp_mac, file.path(results_dir, "gsea/variance_GSEA_GOBP_macula_FIXED.csv"))
write_csv(gsea_gobp_nonmac, file.path(results_dir, "gsea/variance_GSEA_GOBP_nonmacula_FIXED.csv"))
write_csv(gsea_c2_mac, file.path(results_dir, "gsea/variance_GSEA_C2_macula_FIXED.csv"))
write_csv(gsea_c2_nonmac, file.path(results_dir, "gsea/variance_GSEA_C2_nonmacula_FIXED.csv"))

cat("\n")

#%%
# 7. Meta-analysis with gene symbols
cat("7. Running meta-analysis...\n")

fisher_combined_p <- function(p_values) {
  chi_sq <- -2 * sum(log(p_values))
  df <- 2 * length(p_values)
  p_combined <- pchisq(chi_sq, df, lower.tail = FALSE)
  return(p_combined)
}

meta_variance <- results_mac_sym %>%
  select(gene_symbol, levene_pval_mac = levene_pval, var_ratio_mac = var_ratio) %>%
  inner_join(
    results_nonmac_sym %>% select(gene_symbol, levene_pval_nonmac = levene_pval, var_ratio_nonmac = var_ratio),
    by = "gene_symbol"
  ) %>%
  filter(!is.na(levene_pval_mac), !is.na(levene_pval_nonmac)) %>%
  rowwise() %>%
  mutate(
    meta_pval = fisher_combined_p(c(levene_pval_mac, levene_pval_nonmac)),
    avg_var_ratio = (var_ratio_mac + var_ratio_nonmac) / 2
  ) %>%
  ungroup() %>%
  mutate(meta_fdr = p.adjust(meta_pval, method = "BH")) %>%
  arrange(meta_pval)

meta_de <- results_mac_sym %>%
  select(gene_symbol, pval_mac = P.Value, logfc_mac = logFC) %>%
  inner_join(
    results_nonmac_sym %>% select(gene_symbol, pval_nonmac = P.Value, logfc_nonmac = logFC),
    by = "gene_symbol"
  ) %>%
  filter(!is.na(pval_mac), !is.na(pval_nonmac)) %>%
  rowwise() %>%
  mutate(
    meta_pval = fisher_combined_p(c(pval_mac, pval_nonmac)),
    avg_logfc = (logfc_mac + logfc_nonmac) / 2
  ) %>%
  ungroup() %>%
  mutate(meta_fdr = p.adjust(meta_pval, method = "BH")) %>%
  arrange(meta_pval)

write_csv(meta_variance, file.path(results_dir, "meta_analysis_variance_FIXED.csv"))
write_csv(meta_de, file.path(results_dir, "meta_analysis_DE_FIXED.csv"))

cat(sprintf("  Meta variance: %d genes at FDR<0.05\n", sum(meta_variance$meta_fdr < 0.05)))
cat(sprintf("  Meta DE: %d genes at FDR<0.05\n", sum(meta_de$meta_fdr < 0.05)))

cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE!\n")
cat("Summary:\n")
cat(sprintf("  - AMD genes found: %d\n", nrow(amd_combined)))
cat(sprintf("  - Hallmark pathways: Mac=%d, nonMac=%d\n", nrow(gsea_h_mac), nrow(gsea_h_nonmac)))
cat(sprintf("  - GO BP pathways: Mac=%d, nonMac=%d\n", nrow(gsea_gobp_mac), nrow(gsea_gobp_nonmac)))
cat(sprintf("  - C2 pathways: Mac=%d, nonMac=%d\n", nrow(gsea_c2_mac), nrow(gsea_c2_nonmac)))
cat(sprintf("  - Meta-analysis variance genes (FDR<0.05): %d\n", sum(meta_variance$meta_fdr < 0.05)))
cat(sprintf("  - Meta-analysis DE genes (FDR<0.05): %d\n", sum(meta_de$meta_fdr < 0.05)))
cat(rep("=", 70), "\n", sep = "")
