#%%
# GSE135092 vs PRG4 Bulk Data Comparison
# Cross-validation of cohort signatures in H2O2/PRG4 model
# Method: GSEA of Cohort Sets vs PRG4 Ranks

library(dplyr)
library(readr)
library(fgsea)
library(ggplot2)
library(tidyr)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")
output_dir <- file.path(results_dir, "prg4_comparison")
dir.create(output_dir, showWarnings = FALSE)

cat("GSE135092 vs PRG4 Comparison\n")
cat(rep("=", 70), "\n\n", sep = "")

#%%
# 1. Load PRG4 Data & Create Ranks
cat("1. Loading PRG4 bulk data...\n")
prg4_data <- read.csv(file.path(base_dir, "data/RPE_DE_results.csv"))

# H2O2 vs Control Rank
h2o2_rank <- prg4_data %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "" & !is.na(H2O2_vs_CTRL.log2FoldChange)) %>%
  arrange(desc(H2O2_vs_CTRL.log2FoldChange)) %>%
  select(hgnc_symbol, H2O2_vs_CTRL.log2FoldChange) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) 
  
ranks_h2o2 <- setNames(h2o2_rank$H2O2_vs_CTRL.log2FoldChange, h2o2_rank$hgnc_symbol)

# PRG4 Rescue Rank (PRG4+H2O2 vs H2O2)
rescue_rank <- prg4_data %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "" & !is.na(H2O2PRG4_vs_H2O2.log2FoldChange)) %>%
  arrange(desc(H2O2PRG4_vs_H2O2.log2FoldChange)) %>%
  select(hgnc_symbol, H2O2PRG4_vs_H2O2.log2FoldChange) %>%
  distinct(hgnc_symbol, .keep_all = TRUE)

ranks_rescue <- setNames(rescue_rank$H2O2PRG4_vs_H2O2.log2FoldChange, rescue_rank$hgnc_symbol)

cat(sprintf("  H2O2 ranks: %d genes\n", length(ranks_h2o2)))
cat(sprintf("  Rescue ranks: %d genes\n\n", length(ranks_rescue)))

#%%
# 2. Define Cohort Gene Sets
cat("2. Defining GSE135092 Gene Sets...\n")

# Load results
res_mac <- read_csv(file.path(results_dir, "RPE_Macula_all_results_symbols.csv"), show_col_types=F)
res_nonmac <- read_csv(file.path(results_dir, "RPE_nonMacula_all_results_symbols.csv"), show_col_types=F)

# Helper to get genes
get_genes <- function(df, type, direction=NULL) {
  if (type == "DE_UP") {
    return(df %>% filter(P.Value < 0.01, logFC > 0) %>% pull(gene_symbol) %>% unique())
  } else if (type == "DE_DN") {
    return(df %>% filter(P.Value < 0.01, logFC < 0) %>% pull(gene_symbol) %>% unique())
  } else if (type == "VAR_INC") {
    return(df %>% filter(levene_pval < 0.01, var_ratio > 1) %>% pull(gene_symbol) %>% unique())
  } else if (type == "VAR_DEC") {
    return(df %>% filter(levene_pval < 0.01, var_ratio < 1) %>% pull(gene_symbol) %>% unique())
  }
}

gene_sets <- list(
  "Macula_DE_UP" = get_genes(res_mac, "DE_UP"),
  "Macula_DE_DN" = get_genes(res_mac, "DE_DN"),
  "Macula_VAR_INC" = get_genes(res_mac, "VAR_INC"),
  "Macula_VAR_DEC" = get_genes(res_mac, "VAR_DEC"),
  
  "nonMacula_DE_UP" = get_genes(res_nonmac, "DE_UP"),
  "nonMacula_DE_DN" = get_genes(res_nonmac, "DE_DN"),
  "nonMacula_VAR_INC" = get_genes(res_nonmac, "VAR_INC"),
  "nonMacula_VAR_DEC" = get_genes(res_nonmac, "VAR_DEC")
)

# Filter empty or small sets
gene_sets <- gene_sets[sapply(gene_sets, length) >= 5]

cat("  Gene Sets Defined:\n")
print(sapply(gene_sets, length))
cat("\n")

#%%
# 3. Run GSEA
cat("3. Running GSEA (Cohort Sets vs Bulk Ranks)...\n")

run_comparison <- function(ranks, rank_name) {
  gsea_res <- fgsea(pathways = gene_sets, stats = ranks, minSize=5, nproc=1)
  gsea_res$rank_name <- rank_name
  return(gsea_res)
}

res_h2o2 <- run_comparison(ranks_h2o2, "H2O2_vs_Ctrl")
res_rescue <- run_comparison(ranks_rescue, "PRG4_Rescue")

combined_res <- bind_rows(res_h2o2, res_rescue) %>%
  arrange(padj) %>%
  select(pathway, rank_name, NES, pval, padj, size)

print(combined_res)

write_csv(combined_res, file.path(output_dir, "gse135092_prg4_gsea_results.csv"))
cat("\n  Saved gse135092_prg4_gsea_results.csv\n")

#%%
# 4. Visualization
cat("4. Creating Visualization...\n")

# filter for significant or interesting results
plot_data <- combined_res %>%
  mutate(sig = ifelse(padj < 0.05, "*", "")) %>%
  mutate(pathway = factor(pathway, levels = sort(unique(pathway), decreasing = TRUE)))

p <- ggplot(plot_data, aes(x = rank_name, y = pathway, fill = NES)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig), vjust = 0.8, size = 6) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(title = "GSE135092 Signatures in H2O2/PRG4 Model",
       subtitle = "* FDR < 0.05",
       x = "Bulk RNA-seq Contrast",
       y = "GSE135092 Gene Set",
       fill = "NES") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "gse135092_prg4_enrichment_heatmap.pdf"), p, width = 8, height = 6)

ggsave(file.path(output_dir, "gse135092_prg4_enrichment_heatmap.tiff"), p, 
       width = 8, height = 6, dpi = 300, compression = "lzw")
cat("  Saved gse135092_prg4_enrichment_heatmap [.pdf, .tiff]\n")

cat("\nComparison Complete!\n")
