#%%
# GSE135092 Variance Gene Follow-Up Analysis
# Comprehensive downstream analysis following GSE29801 framework

library(dplyr)
library(readr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(mclust)
library(VennDiagram)
library(grid)
library(gridExtra)
library(patchwork)
library(tidyr)

# Manual Fisher's combined probability test
fisher_combined_p <- function(p_values) {
  chi_sq <- -2 * sum(log(p_values))
  df <- 2 * length(p_values)
  p_combined <- pchisq(chi_sq, df, lower.tail = FALSE)
  return(p_combined)
}

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")
out_dir <- results_dir  # Add to existing directory

# Create output subdirectories
dir.create(file.path(out_dir, "gsea"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "bimodality"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

cat("GSE135092 Variance Follow-Up Analysis\n")
cat(rep("=", 70), "\n\n", sep = "")

#%%
# Custom theme
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "serif"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

#%%
# 1. Load existing results
cat("1. Loading existing results...\n")
results_mac <- read_csv(file.path(results_dir, "RPE_Macula_all_results.csv"), 
                        show_col_types = FALSE)
results_nonmac <- read_csv(file.path(results_dir, "RPE_nonMacula_all_results.csv"),
                           show_col_types = FALSE)

cat(sprintf("  RPE Macula: %d genes\n", nrow(results_mac)))
cat(sprintf("  RPE non-Macula: %d genes\n", nrow(results_nonmac)))

# Extract variance genes (p < 0.01)
var_genes_mac <- results_mac %>%
  filter(levene_pval < 0.01, !is.na(levene_pval)) %>%
  arrange(levene_pval)

var_genes_nonmac <- results_nonmac %>%
  filter(levene_pval < 0.01, !is.na(levene_pval)) %>%
  arrange(levene_pval)

cat(sprintf("\n  Variance genes (p<0.01): Macula=%d, non-Macula=%d\n", 
            nrow(var_genes_mac), nrow(var_genes_nonmac)))

# Extract DE genes (p < 0.01) for comparison
de_genes_mac <- results_mac %>%
  filter(P.Value < 0.01, !is.na(P.Value)) %>%
  arrange(P.Value)

de_genes_nonmac <- results_nonmac %>%
  filter(P.Value < 0.01, !is.na(P.Value)) %>%
  arrange(P.Value)

cat(sprintf("  DE genes (p<0.01): Macula=%d, non-Macula=%d\n\n", 
            nrow(de_genes_mac), nrow(de_genes_nonmac)))

#%%
# 2. Download MSigDB gene sets
cat("2. Downloading MSigDB gene sets...\n")

# Hallmark
msig_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_sets <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)
cat(sprintf("  Hallmark: %d gene sets\n", length(hallmark_sets)))

# GO Biological Process (subset to avoid too many terms)
msig_gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
gobp_sets <- split(msig_gobp$gene_symbol, msig_gobp$gs_name)
cat(sprintf("  GO BP: %d gene sets\n", length(gobp_sets)))

# C2 Curated (KEGG, Reactome, BioCarta, etc.)
msig_c2 <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP")
c2_sets <- split(msig_c2$gene_symbol, msig_c2$gs_name)
cat(sprintf("  C2 Curated: %d gene sets\n\n", length(c2_sets)))

#%%
# 3. Function to run GSEA
run_gsea_variance <- function(results_df, gene_sets, tissue_name, db_name) {
  cat(sprintf("Running GSEA: %s - %s\n", tissue_name, db_name))
  
  # Create ranked list based on variance p-value
  # Use -log10(p) * sign(var_ratio - 1) for ranking
  ranked_list <- results_df %>%
    filter(!is.na(levene_pval), !is.na(var_ratio)) %>%
    mutate(
      rank_stat = -log10(levene_pval) * sign(var_ratio - 1)
    ) %>%
    arrange(desc(rank_stat))
  
  ranks <- setNames(ranked_list$rank_stat, ranked_list$gene)
  
  # Run GSEA
  gsea_res <- fgsea(pathways = gene_sets,
                    stats = ranks,
                    minSize = 15,
                    maxSize = 500,
                    nperm = 10000)
  
  # Filter and sort
  gsea_res <- gsea_res %>%
    filter(padj < 0.25) %>%  # Relaxed threshold for initial exploration
    arrange(padj, desc(abs(NES)))
  
  cat(sprintf("  Found %d pathways (FDR<0.25)\n", nrow(gsea_res)))
  
  return(gsea_res)
}

#%%
# 4. Run GSEA on all combinations
cat("\n3. Running GSEA on variance genes...\n")

# Hallmark
gsea_hallmark_mac <- run_gsea_variance(results_mac, hallmark_sets, 
                                       "RPE Macula", "Hallmark")
gsea_hallmark_nonmac <- run_gsea_variance(results_nonmac, hallmark_sets,
                                          "RPE non-Macula", "Hallmark")

write_csv(gsea_hallmark_mac, file.path(out_dir, "gsea/variance_GSEA_Hallmark_macula.csv"))
write_csv(gsea_hallmark_nonmac, file.path(out_dir, "gsea/variance_GSEA_Hallmark_nonmacula.csv"))

# GO BP (limit for runtime)
gsea_gobp_mac <- run_gsea_variance(results_mac, gobp_sets,
                                   "RPE Macula", "GO BP")
gsea_gobp_nonmac <- run_gsea_variance(results_nonmac, gobp_sets,
                                      "RPE non-Macula", "GO BP")

write_csv(gsea_gobp_mac, file.path(out_dir, "gsea/variance_GSEA_GOBP_macula.csv"))
write_csv(gsea_gobp_nonmac, file.path(out_dir, "gsea/variance_GSEA_GOBP_nonmacula.csv"))

# C2 Curated
gsea_c2_mac <- run_gsea_variance(results_mac, c2_sets,
                                 "RPE Macula", "C2 Curated")
gsea_c2_nonmac <- run_gsea_variance(results_nonmac, c2_sets,
                                    "RPE non-Macula", "C2 Curated")

write_csv(gsea_c2_mac, file.path(out_dir, "gsea/variance_GSEA_C2_macula.csv"))
write_csv(gsea_c2_nonmac, file.path(out_dir, "gsea/variance_GSEA_C2_nonmacula.csv"))

cat("\n")

#%%
# 5. Known AMD genes variance analysis
cat("4. Analyzing known AMD genes variance...\n")

amd_genes <- c("CFH", "ARMS2", "C2", "C3", "CFB", "CFI", "HTRA1", 
               "TIMP3", "COL8A1", "BEST1", "RPE65", "ABCA4", 
               "EFEMP1", "APOE", "CETP", "LIPC", "VEGFA", "HMCN1")

amd_var_analysis <- function(results_df, tissue_name) {
  results_df %>%
    filter(gene %in% amd_genes) %>%
    select(gene, ctrl_var, amd_var, var_ratio, levene_pval, levene_fdr) %>%
    mutate(
      tissue = tissue_name,
      var_change = ifelse(var_ratio > 1, "Increased", "Decreased"),
      var_fc = log2(var_ratio)
    ) %>%
    arrange(levene_pval)
}

amd_var_mac <- amd_var_analysis(results_mac, "RPE Macula")
amd_var_nonmac <- amd_var_analysis(results_nonmac, "RPE non-Macula")

amd_var_combined <- bind_rows(amd_var_mac, amd_var_nonmac)
write_csv(amd_var_combined, file.path(out_dir, "AMD_genes_variance.csv"))

cat(sprintf("  Found %d AMD genes with data\n", 
            length(unique(amd_var_combined$gene))))
cat(sprintf("  Significant (p<0.01): Macula=%d, non-Macula=%d\n\n",
            sum(amd_var_mac$levene_pval < 0.01, na.rm = TRUE),
            sum(amd_var_nonmac$levene_pval < 0.01, na.rm = TRUE)))

#%%
# 6. Variance vs DE overlap analysis
cat("5. Variance vs DE overlap analysis...\n")

overlap_analysis <- function(var_df, de_df, tissue_name) {
  var_set <- var_df$gene
  de_set <- de_df$gene
  overlap <- intersect(var_set, de_set)
  
  # Fisher's exact test
  total_genes <- 58302  # Total genes in dataset
  contingency <- matrix(c(
    length(overlap),  # Both variance and DE
    length(var_set) - length(overlap),  # Variance only
    length(de_set) - length(overlap),  # DE only
    total_genes - length(var_set) - length(de_set) + length(overlap)  # Neither
  ), nrow = 2)
  
  fisher_test <- fisher.test(contingency)
  
  list(
    tissue = tissue_name,
    n_variance = length(var_set),
    n_de = length(de_set),
    n_overlap = length(overlap),
    overlap_genes = overlap,
    fisher_pval = fisher_test$p.value,
    fisher_OR = fisher_test$estimate
  )
}

overlap_mac <- overlap_analysis(var_genes_mac, de_genes_mac, "RPE Macula")
overlap_nonmac <- overlap_analysis(var_genes_nonmac, de_genes_nonmac, "RPE non-Macula")

# Save overlap summary
overlap_summary <- data.frame(
  Tissue = c(overlap_mac$tissue, overlap_nonmac$tissue),
  Variance_genes = c(overlap_mac$n_variance, overlap_nonmac$n_variance),
  DE_genes = c(overlap_mac$n_de, overlap_nonmac$n_de),
  Overlap = c(overlap_mac$n_overlap, overlap_nonmac$n_overlap),
  Fisher_pval = c(overlap_mac$fisher_pval, overlap_nonmac$fisher_pval),
  Fisher_OR = c(overlap_mac$fisher_OR, overlap_nonmac$fisher_OR)
)

write_csv(overlap_summary, file.path(out_dir, "variance_DE_overlap.csv"))

cat(sprintf("  Macula: %d variance, %d DE, %d overlap (OR=%.2f, p=%.3f)\n",
            overlap_mac$n_variance, overlap_mac$n_de, overlap_mac$n_overlap,
            overlap_mac$fisher_OR, overlap_mac$fisher_pval))
cat(sprintf("  non-Macula: %d variance, %d DE, %d overlap (OR=%.2f, p=%.3f)\n\n",
            overlap_nonmac$n_variance, overlap_nonmac$n_de, overlap_nonmac$n_overlap,
            overlap_nonmac$fisher_OR, overlap_nonmac$fisher_pval))

#%%
# 7. Bimodality testing
cat("6. Bimodality testing (top 50 variance genes per tissue)...\n")

# Load expression data for bimodality testing
expr_cache <- readRDS(file.path(base_dir, "results/cohort-GSE135092/expression_matrix_cache.rds"))
expr_mat <- expr_cache$expr_mat
meta <- read_csv(file.path(base_dir, "results/cohort-GSE135092/GSE135092_sample_metadata.csv"),
                 show_col_types = FALSE)

test_bimodality <- function(var_genes_df, expr_mat, meta, tissue_name, n_top = 50) {
  cat(sprintf("  Testing %s...\n", tissue_name))
  
  # Get top N variance genes
  top_genes <- head(var_genes_df, n_top)$gene
  
  # Filter metadata for this tissue and AMD samples only
  tissue_pattern <- tissue_name
  meta_tissue <- meta %>%
    filter(tissue == tissue_pattern, amd_status == "AMD", !is.na(age))
  
  results <- data.frame(
    gene = character(),
    bic_1comp = numeric(),
    bic_2comp = numeric(),
    bic_improvement = numeric(),
    is_bimodal = logical(),
    stringsAsFactors = FALSE
  )
  
  for (gene in top_genes) {
    if (!gene %in% rownames(expr_mat)) next
    
    expr_vals <- as.numeric(expr_mat[gene, meta_tissue$gsm])
    
    tryCatch({
      # 1-component model
      fit1 <- Mclust(expr_vals, G = 1, verbose = FALSE)
      # 2-component model
      fit2 <- Mclust(expr_vals, G = 2, verbose = FALSE)
      
      bic1 <- fit1$bic
      bic2 <- fit2$bic
      bic_improve <- bic2 - bic1
      
      results <- rbind(results, data.frame(
        gene = gene,
        bic_1comp = bic1,
        bic_2comp = bic2,
        bic_improvement = bic_improve,
        is_bimodal = bic_improve > 10,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      # Skip genes that fail
    })
  }
  
  return(results)
}

bimodal_mac <- test_bimodality(var_genes_mac, expr_mat, meta, "RPE, Macula")
bimodal_nonmac <- test_bimodality(var_genes_nonmac, expr_mat, meta, "RPE, non-Macula")

write_csv(bimodal_mac, file.path(out_dir, "bimodality/bimodality_results_macula.csv"))
write_csv(bimodal_nonmac, file.path(out_dir, "bimodality/bimodality_results_nonmacula.csv"))

cat(sprintf("  Macula: %d/%d genes bimodal (BIC improvement > 10)\n",
            sum(bimodal_mac$is_bimodal, na.rm = TRUE), nrow(bimodal_mac)))
cat(sprintf("  non-Macula: %d/%d genes bimodal\n\n",
            sum(bimodal_nonmac$is_bimodal, na.rm = TRUE), nrow(bimodal_nonmac)))

#%%
# 8. Variance direction analysis
cat("7. Variance direction analysis...\n")

var_direction_summary <- data.frame(
  Tissue = c("RPE Macula", "RPE non-Macula"),
  Increased_variance = c(
    sum(var_genes_mac$var_ratio > 1, na.rm = TRUE),
    sum(var_genes_nonmac$var_ratio > 1, na.rm = TRUE)
  ),
  Decreased_variance = c(
    sum(var_genes_mac$var_ratio < 1, na.rm = TRUE),
    sum(var_genes_nonmac$var_ratio < 1, na.rm = TRUE)
  )
)

write_csv(var_direction_summary, file.path(out_dir, "variance_direction_summary.csv"))
print(var_direction_summary)

cat("\n")

#%%
# 9. Meta-analysis
cat("8. Meta-analysis (combining both RPE tissues)...\n")

# For variance: combine p-values using Fisher's method
meta_variance <- results_mac %>%
  select(gene, levene_pval_mac = levene_pval, var_ratio_mac = var_ratio) %>%
  inner_join(
    results_nonmac %>% select(gene, levene_pval_nonmac = levene_pval, var_ratio_nonmac = var_ratio),
    by = "gene"
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

write_csv(meta_variance, file.path(out_dir, "meta_analysis_variance.csv"))

meta_var_sig <- meta_variance %>% filter(meta_fdr < 0.05)
cat(sprintf("  Meta-analysis variance: %d genes at FDR<0.05\n", nrow(meta_var_sig)))

# For DE: similar approach
meta_de <- results_mac %>%
  select(gene, pval_mac = P.Value, logfc_mac = logFC) %>%
  inner_join(
    results_nonmac %>% select(gene, pval_nonmac = P.Value, logfc_nonmac = logFC),
    by = "gene"
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

write_csv(meta_de, file.path(out_dir, "meta_analysis_DE.csv"))

meta_de_sig <- meta_de %>% filter(meta_fdr < 0.05)
cat(sprintf("  Meta-analysis DE: %d genes at FDR<0.05\n\n", nrow(meta_de_sig)))

cat("\n", rep("=", 70), "\n", sep = "")
cat("Analysis complete! Results saved to:\n")
cat(sprintf("  %s\n", out_dir))
cat("\nNext: Creating visualizations...\n")
cat(rep("=", 70), "\n\n", sep = "")
