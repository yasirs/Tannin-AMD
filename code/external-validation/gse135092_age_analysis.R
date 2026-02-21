#%%
# GSE135092 Age DE and C2 GSEA Analysis
# Focus: Extract Age effects (UP with Age) and run GSEA (Hallmark, GO BP, C2)

library(limma)
library(dplyr)
library(readr)
library(fgsea)
library(msigdbr)
library(ggplot2)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092")
out_dir <- file.path(results_dir, "rpe_covariate_de")

# Load cached expression and metadata
cat("Loading data...\n")
expr_data <- readRDS(file.path(results_dir, "expression_matrix_cache.rds"))
expr_mat <- expr_data$expr_mat
meta_all <- read_csv(file.path(results_dir, "GSE135092_sample_metadata.csv"), show_col_types = FALSE)

#%%
# Function to run Limma ONLY (Fast)
run_limma_age <- function(expr, meta, tissue_name) {
  cat(sprintf("\nRunning Limma for %s...\n", tissue_name))
  
  meta_filt <- meta %>%
    filter(!is.na(age), amd_status %in% c("Control", "AMD"))
  
  expr_filt <- expr[, meta_filt$gsm]
  
  # Design: Age + AMD_Status
  design <- model.matrix(~ age + amd_status, data = meta_filt)
  
  fit <- lmFit(expr_filt, design)
  fit <- eBayes(fit)
  
  # Extract Age coefficient
  res_age <- topTable(fit, coef = "age", number = Inf, sort.by = "none")
  res_age$gene <- rownames(res_age)
  
  # Map Symbols (using existing mapping file to save time)
  mapping <- read_csv(file.path(out_dir, "RPE_Macula_all_results_symbols.csv"), 
                      col_select = c("gene", "gene_symbol"), show_col_types = FALSE)
  
  res_age <- res_age %>%
    left_join(mapping, by = "gene") %>%
    filter(!is.na(gene_symbol), gene_symbol != "")
    
  return(res_age)
}

#%%
# 1. Extract Age DE Results
cat("1. Extracting Age DE Results...\n")

meta_mac <- meta_all %>% filter(tissue == "RPE, Macula")
meta_nonmac <- meta_all %>% filter(tissue == "RPE, non-Macula")

age_mac <- run_limma_age(expr_mat[, meta_mac$gsm], meta_mac, "RPE Macula")
age_nonmac <- run_limma_age(expr_mat[, meta_nonmac$gsm], meta_nonmac, "RPE non-Macula")

write_csv(age_mac, file.path(out_dir, "age_de_macula.csv"))
write_csv(age_nonmac, file.path(out_dir, "age_de_nonmacula.csv"))

cat(sprintf("  Macula: %d genes\n", nrow(age_mac)))
cat(sprintf("  non-Macula: %d genes\n\n", nrow(age_nonmac)))

#%%
# 2. Prepare GSEA Ranks (UP with Age)
# Rank by t-statistic for Age coefficient
ranks_mac <- setNames(age_mac$t, age_mac$gene_symbol)
ranks_mac <- sort(ranks_mac, decreasing = TRUE)

ranks_nonmac <- setNames(age_nonmac$t, age_nonmac$gene_symbol)
ranks_nonmac <- sort(ranks_nonmac, decreasing = TRUE)

#%%
# 3. Load Gene Sets (Hallmark, GO BP, C2)
cat("2. Loading MSigDB Sets (H, C5, C2)...\n")

msig_h <- msigdbr(species = "Homo sapiens", collection = "H")
sets_h <- split(msig_h$gene_symbol, msig_h$gs_name)

msig_gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
sets_gobp <- split(msig_gobp$gene_symbol, msig_gobp$gs_name)

msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
sets_c2 <- split(msig_c2$gene_symbol, msig_c2$gs_name)

cat(sprintf("  Hallmark: %d\n", length(sets_h)))
cat(sprintf("  GO BP: %d\n", length(sets_gobp)))
cat(sprintf("  C2 Curated: %d\n\n", length(sets_c2)))

#%%
# 4. Run GSEA
run_gsea_age <- function(ranks, sets, name) {
  cat(sprintf("  Running GSEA: %s...\n", name))
  res <- fgsea(pathways = sets, stats = ranks, minSize = 15, maxSize = 500, nproc = 1)
  res <- res %>% arrange(pval) %>% filter(padj < 0.25)
  return(res)
}

cat("3. Running GSEA for Age UP (positive enrichment)...\n")

# Macula
gsea_mac_h <- run_gsea_age(ranks_mac, sets_h, "Macula - Hallmark")
gsea_mac_gobp <- run_gsea_age(ranks_mac, sets_gobp, "Macula - GO BP")
gsea_mac_c2 <- run_gsea_age(ranks_mac, sets_c2, "Macula - C2")

# non-Macula
gsea_nonmac_h <- run_gsea_age(ranks_nonmac, sets_h, "non-Macula - Hallmark")
gsea_nonmac_gobp <- run_gsea_age(ranks_nonmac, sets_gobp, "non-Macula - GO BP")
gsea_nonmac_c2 <- run_gsea_age(ranks_nonmac, sets_c2, "non-Macula - C2")

# Save Results
cat("\nSaving GSEA Results...\n")
out_gsea <- file.path(out_dir, "gsea")
write_csv(gsea_mac_h, file.path(out_gsea, "Age_GSEA_Hallmark_macula.csv"))
write_csv(gsea_mac_gobp, file.path(out_gsea, "Age_GSEA_GOBP_macula.csv"))
write_csv(gsea_mac_c2, file.path(out_gsea, "Age_GSEA_C2_macula.csv"))

write_csv(gsea_nonmac_h, file.path(out_gsea, "Age_GSEA_Hallmark_nonmacula.csv"))
write_csv(gsea_nonmac_gobp, file.path(out_gsea, "Age_GSEA_GOBP_nonmacula.csv"))
write_csv(gsea_nonmac_c2, file.path(out_gsea, "Age_GSEA_C2_nonmacula.csv"))

#%%
# 5. Visualization (Barplots)
cat("4. Creating Visualizations...\n")

plot_gsea <- function(df, title) {
  if(nrow(df) == 0) return(NULL)
  
  # Filter for positive NES (UP with Age)
  df_up <- df %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(10)
  
  if(nrow(df_up) == 0) return(NULL)
  
  p <- ggplot(df_up, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "red", high = "blue", name = "FDR") +
    labs(title = title, x = "", y = "NES") +
    theme_minimal()
  return(p)
}

figs_dir <- file.path(out_dir, "figures")

# Macula Plots
p1 <- plot_gsea(gsea_mac_h, "Age UP - Macula (Hallmark)")
p2 <- plot_gsea(gsea_mac_c2, "Age UP - Macula (C2)")

# non-Macula Plots
p3 <- plot_gsea(gsea_nonmac_h, "Age UP - non-Macula (Hallmark)")
p4 <- plot_gsea(gsea_nonmac_c2, "Age UP - non-Macula (C2)")

# Save Combined
library(gridExtra)
plist <- list(p1, p2, p3, p4)
plist <- plist[!sapply(plist, is.null)]

if(length(plist) > 0) {
  p_combined <- grid.arrange(grobs = plist, ncol = 2)
  ggsave(file.path(figs_dir, "Age_GSEA_summary.pdf"), p_combined, width = 24, height = 12)
  ggsave(file.path(figs_dir, "Age_GSEA_summary.tiff"), p_combined, 
         width = 24, height = 12, dpi = 300, compression = "lzw")
  cat("Saved Age_GSEA_summary [.pdf, .tiff]\n")
} else {
  cat("No significant Age UP pathways found to plot.\n")
}

cat("\nAnalysis Complete!\n")
