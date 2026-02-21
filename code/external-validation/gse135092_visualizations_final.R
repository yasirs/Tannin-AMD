#%%
# GSE135092 Visualization Script - Comprehensive
# 1. GSEA Barplots (DE genes)
# 2. ORA Barplots (Directional Variance)
# 3. Venn Diagrams (Variance vs DE)
# 4. AMD Genes Heatmap
# 5. Comprehensive Summary Figure

library(ggplot2)
library(dplyr)
library(readr)
library(pheatmap)
library(grid)
library(gridExtra)
library(VennDiagram)
library(RColorBrewer)
library(cowplot)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")
plot_dir <- file.path(results_dir, "figures")
dir.create(plot_dir, showWarnings = FALSE)

# Load Results
cat("Loading results...\n")

# DE GSEA Results (Hallmark)
de_gsea_mac <- read_csv(file.path(results_dir, "gsea/DE_GSEA_Hallmark_macula.csv"), show_col_types = FALSE)
de_gsea_nonmac <- read_csv(file.path(results_dir, "gsea/DE_GSEA_Hallmark_nonmacula.csv"), show_col_types = FALSE)

# Variance ORA Results (Decreased - Hallmark & GO BP)
ora_dec_hall_mac <- read_csv(file.path(results_dir, "gsea/ORA_decreased_Hallmark_macula.csv"), show_col_types = FALSE)
ora_dec_hall_nonmac <- read_csv(file.path(results_dir, "gsea/ORA_decreased_Hallmark_nonmacula.csv"), show_col_types = FALSE)
ora_dec_gobp_mac <- read_csv(file.path(results_dir, "gsea/ORA_decreased_GOBP_macula.csv"), show_col_types = FALSE)
ora_dec_gobp_nonmac <- read_csv(file.path(results_dir, "gsea/ORA_decreased_GOBP_nonmacula.csv"), show_col_types = FALSE)

# AMD Variance Genes
amd_genes <- read_csv(file.path(results_dir, "AMD_genes_variance_FIXED.csv"), show_col_types = FALSE)

# Age GSEA Results (Hallmark)
age_gsea_mac <- read_csv(file.path(results_dir, "gsea/Age_GSEA_Hallmark_macula.csv"), show_col_types = FALSE)
age_gsea_nonmac <- read_csv(file.path(results_dir, "gsea/Age_GSEA_Hallmark_nonmacula.csv"), show_col_types = FALSE)

# Variance & DE Gene Lists (for Venn)
summary_ora <- read_csv(file.path(results_dir, "ORA_directional_summary.csv"), show_col_types = FALSE)

cat("Results loaded.\n\n")

#%%
# Helper function: GSEA Barplot
plot_gsea_bar <- function(df, title, top_n = 10, fill_var = "NES") {
  df_top <- df %>%
    head(top_n) %>%
    mutate(pathway = gsub("HALLMARK_", "", pathway)) %>%
    mutate(pathway = gsub("GO_", "", pathway)) %>%
    mutate(pathway = gsub("_", " ", pathway)) %>%
    mutate(pathway = reorder(pathway, !!sym(fill_var)))
  
  p <- ggplot(df_top, aes(x = pathway, y = !!sym(fill_var), fill = padj)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "red", high = "blue", name = "FDR") +
    labs(title = title, x = "", y = fill_var) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"))
  return(p)
}

# Helper function: ORA Barplot
plot_ora_bar <- function(df, title, top_n = 10) {
  df_top <- df %>%
    head(top_n) %>%
    mutate(pathway = gsub("HALLMARK_", "", pathway)) %>%
    mutate(pathway = gsub("GO_", "", pathway)) %>%
    mutate(pathway = gsub("_", " ", pathway)) %>%
    mutate(log_p = -log10(pval)) %>%
    mutate(pathway = reorder(pathway, log_p))
  
  p <- ggplot(df_top, aes(x = pathway, y = log_p, fill = fdr)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "red", high = "blue", name = "FDR") +
    labs(title = title, x = "", y = "-log10(p-value)") +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"))
  return(p)
}

#%%
# 1. Generate Enrichment Plots
cat("1. Generating Enrichment Plots...\n")

p_de_mac <- plot_gsea_bar(de_gsea_mac, "DE (Mean) - Macula\n(Hallmark)", fill_var = "NES")
p_de_nonmac <- plot_gsea_bar(de_gsea_nonmac, "DE (Mean) - non-Macula\n(Hallmark)", fill_var = "NES")

p_ora_mac <- plot_ora_bar(ora_dec_hall_mac, "Decr. Variance - Macula\n(Hallmark)")
p_ora_nonmac <- plot_ora_bar(ora_dec_gobp_nonmac, "Decr. Variance - non-Macula\n(GO BP - Top 10)")

p_age_mac <- plot_gsea_bar(age_gsea_mac, "Age Up - Macula\n(Hallmark)", fill_var = "NES")
p_age_nonmac <- plot_gsea_bar(age_gsea_nonmac, "Age Up - non-Macula\n(Hallmark)", fill_var = "NES")

# Save combined enrichment plot
# Row 1: Disease DE (Mean)
# Row 2: Disease Variance (Decreased)
# Row 3: Age Effect (Up with Age)
p_grid <- plot_grid(
  p_de_mac, p_de_nonmac, 
  p_ora_mac, p_ora_nonmac, 
  p_age_mac, p_age_nonmac,
  ncol = 2, labels = "AUTO"
)

# Save as PDF (Vector)
ggsave(file.path(plot_dir, "enrichment_summary_v2.pdf"), p_grid, width = 21, height = 18)

# Save as TIFF (Raster, High Res, LZW)
ggsave(file.path(plot_dir, "enrichment_summary_v2.tiff"), p_grid, 
       width = 21, height = 18, dpi = 300, compression = "lzw")
cat("  Saved enrichment_summary [.pdf, .tiff]\n")

#%%
# 2. AMD Genes Variance Heatmap
cat("2. Generating AMD Genes Heatmap...\n")

# Prepare data matrix
amd_mat <- amd_genes %>%
  select(gene_symbol, tissue, var_ratio) %>%
  mutate(log2_var_ratio = log2(var_ratio)) %>%
  tidyr::pivot_wider(names_from = tissue, values_from = log2_var_ratio) 

mat_data <- as.matrix(amd_mat[, -1])
rownames(mat_data) <- amd_mat$gene_symbol

# Save heatmap (TIFF)
tiff(file.path(plot_dir, "AMD_genes_variance_heatmap.tiff"), width = 6, height = 8, units = "in", res = 300, compression = "lzw")
pheatmap(mat_data, 
         main = "AMD Risk Genes\nLog2 Variance Ratio (AMD/Control)",
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         display_numbers = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-2, 2, length.out = 101))
dev.off()
cat("  Saved AMD_genes_variance_heatmap.tiff\n")

#%%
# 3. Venn Diagrams (Variance vs DE)
cat("3. Generating Venn Diagrams...\n")

# Need to reload original gene lists for exact counts
res_mac <- read_csv(file.path(results_dir, "RPE_Macula_all_results_symbols.csv"), show_col_types=F)
res_nonmac <- read_csv(file.path(results_dir, "RPE_nonMacula_all_results_symbols.csv"), show_col_types=F)

# Filter lists (p < 0.01)
get_lists <- function(df) {
  var_genes <- df %>% filter(levene_pval < 0.01) %>% pull(gene_symbol)
  de_genes <- df %>% filter(P.Value < 0.01) %>% pull(gene_symbol)
  return(list(var = unique(var_genes), de = unique(de_genes)))
}

lists_mac <- get_lists(res_mac)
lists_nonmac <- get_lists(res_nonmac)

# Macula Venn (TIFF)
tiff(file.path(plot_dir, "venn_macula.tiff"), width = 6, height = 6, units = "in", res = 300, compression = "lzw")
grid.draw(venn.diagram(
  x = list(Variance = lists_mac$var, DE = lists_mac$de),
  filename = NULL,
  fill = c("blue", "red"), alpha = 0.5,
  main = "RPE Macula: Variance vs DE"
))
dev.off()

# non-Macula Venn (TIFF)
tiff(file.path(plot_dir, "venn_nonmacula.tiff"), width = 6, height = 6, units = "in", res = 300, compression = "lzw")
grid.draw(venn.diagram(
  x = list(Variance = lists_nonmac$var, DE = lists_nonmac$de),
  filename = NULL,
  fill = c("blue", "red"), alpha = 0.5,
  main = "RPE non-Macula: Variance vs DE"
))
dev.off()
cat("  Saved venn diagrams [.tiff]\n")

#%%
# Summary
cat("\nVisualization script completed successfully!\n")
