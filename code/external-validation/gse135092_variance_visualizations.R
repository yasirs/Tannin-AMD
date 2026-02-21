#%%
# GSE135092 Variance Follow-Up Visualizations
# Creates multi-panel combined figures

library(dplyr)
library(readr)
library(ggplot2)
library(VennDiagram)
library(grid)
library(gridExtra)
library(patchwork)
library(tidyr)
library(pheatmap)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")

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

cat("Creating visualizations...\n\n")

#%%
# 1. GSEA barplot (top pathways from each database)
cat("1. GSEA pathway enrichment barplot...\n")

# Load GSEA results
gsea_h_mac <- read_csv(file.path(results_dir, "gsea/variance_GSEA_Hallmark_macula.csv"), show_col_types = FALSE)
gsea_h_nonmac <- read_csv(file.path(results_dir, "gsea/variance_GSEA_Hallmark_nonmacula.csv"), show_col_types = FALSE)

# Combine and select top pathways
top_n <- 10
gsea_combined <- bind_rows(
  gsea_h_mac %>% slice_head(n = top_n) %>% mutate(Tissue = "RPE Macula", Database = "Hallmark"),
  gsea_h_nonmac %>% slice_head(n = top_n) %>% mutate(Tissue = "RPE non-Macula", Database = "Hallmark")
)

if (nrow(gsea_combined) > 0) {
  p_gsea <- ggplot(gsea_combined, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ Tissue, scales = "free_y") +
    scale_fill_gradient(low = "#D6604D", high = "#4393C3", name = "FDR") +
    labs(
      title = "Top Variance-Associated Pathways (GSEA)",
      x = NULL,
      y = "Normalized Enrichment Score (NES)"
    ) +
    theme_pub() +
    theme(
      axis.text.y = element_text(size = 8),
      strip.background = element_rect(fill = "white", color = "black")
    )
  
  ggsave(file.path(results_dir, "figures/variance_GSEA_barplot.pdf"), p_gsea, 
         width = 10, height = 8)
  ggsave(file.path(results_dir, "figures/variance_GSEA_barplot.tiff"), p_gsea,
         width = 10, height = 8, dpi = 300, compression = "lzw")
  
  cat("  Saved GSEA barplot\n")
} else {
  cat("  No significant pathways to plot\n")
}

#%%
# 2. Variance vs DE Venn diagrams
cat("2. Variance vs DE overlap Venn diagrams...\n")

overlap_summary <- read_csv(file.path(results_dir, "variance_DE_overlap.csv"), show_col_types = FALSE)

# Function to create Venn diagram
create_venn <- function(n_var, n_de, n_overlap, tissue_name) {
  grid.newpage()
  venn_plot <- draw.pairwise.venn(
    area1 = n_var,
    area2 = n_de,
    cross.area = n_overlap,
    category = c("Variance\n(p<0.01)", "DE\n(p<0.01)"),
    fill = c("#FC8D62", "#8DA0CB"),
    alpha = 0.5,
    cat.pos = c(0, 0),
    cat.dist = 0.03,
    cex = 1.5,
    cat.cex = 1.2,
    fontfamily = rep("serif", 3),
    cat.fontfamily = rep("serif", 2)
  )
  return(venn_plot)
}

# Macula
pdf(file.path(results_dir, "figures/variance_DE_venn_macula.pdf"), width = 6, height = 6)
create_venn(
  overlap_summary$Variance_genes[1],
  overlap_summary$DE_genes[1],
  overlap_summary$Overlap[1],
  "RPE Macula"
)
dev.off()

# non-Macula
pdf(file.path(results_dir, "figures/variance_DE_venn_nonmacula.pdf"), width = 6, height = 6)
create_venn(
  overlap_summary$Variance_genes[2],
  overlap_summary$DE_genes[2],
  overlap_summary$Overlap[2],
  "RPE non-Macula"
)
dev.off()

cat("  Saved Venn diagrams\n")

#%%
# 3. AMD genes variance heatmap
cat("3. AMD genes variance heatmap...\n")

amd_var <- read_csv(file.path(results_dir, "AMD_genes_variance.csv"), show_col_types = FALSE)

if (nrow(amd_var) > 0) {
  # Reshape for heatmap
  amd_mat <- amd_var %>%
    select(gene, tissue, var_fc) %>%
    pivot_wider(names_from = tissue, values_from = var_fc) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # Create annotation for significance
  amd_annot <- amd_var %>%
    mutate(sig = ifelse(levene_pval < 0.01, "*", "")) %>%
    select(gene, tissue, sig) %>%
    pivot_wider(names_from = tissue, values_from = sig) %>%
    column_to_rownames("gene")
  
  pdf(file.path(results_dir, "figures/AMD_genes_variance_heatmap.pdf"), width = 5, height = 8)
  pheatmap(amd_mat,
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           color = colorRampPalette(c("#4393C3", "white", "#D6604D"))(100),
           breaks = seq(-2, 2, length.out = 101),
           main = "AMD Risk Genes: Variance Change in AMD",
           fontsize_row = 10,
           fontsize_col = 10,
           display_numbers = amd_annot,
           number_color = "black",
           fontsize_number = 12)
  dev.off()
  
  cat("  Saved AMD genes heatmap\n")
} else {
  cat("  No AMD genes to plot\n")
}

#%%
# 4. Bimodality summary
cat("4. Bimodality results visualization...\n")

bimodal_mac <- read_csv(file.path(results_dir, "bimodality/bimodality_results_macula.csv"), 
                        show_col_types = FALSE)
bimodal_nonmac <- read_csv(file.path(results_dir, "bimodality/bimodality_results_nonmacula.csv"),
                           show_col_types = FALSE)

bimodal_combined <- bind_rows(
  bimodal_mac %>% mutate(Tissue = "RPE Macula"),
  bimodal_nonmac %>% mutate(Tissue = "RPE non-Macula")
)

if (nrow(bimodal_combined) > 0) {
  p_bimodal <- ggplot(bimodal_combined, aes(x = bic_improvement, fill = Tissue)) +
    geom_histogram(bins = 20, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
    facet_wrap(~ Tissue, ncol = 1) +
    scale_fill_manual(values = c("RPE Macula" = "#FC8D62", "RPE non-Macula" = "#8DA0CB")) +
    labs(
      title = "Bimodality Testing Results",
      subtitle = "Top 50 variance genes per tissue",
      x = "BIC Improvement (2-comp vs 1-comp)",
      y = "Number of Genes"
    ) +
    theme_pub() +
    theme(legend.position = "none")
  
  ggsave(file.path(results_dir, "figures/bimodality_summary.pdf"), p_bimodal,
         width = 8, height = 6)
  ggsave(file.path(results_dir, "figures/bimodality_summary.tiff"), p_bimodal,
         width = 8, height = 6, dpi = 300, compression = "lzw")
  
  cat("  Saved bimodality summary\n")
}

#%%
# 5. Variance direction summary
cat("5. Variance direction barplot...\n")

var_dir <- read_csv(file.path(results_dir, "variance_direction_summary.csv"), show_col_types = FALSE)

var_dir_long <- var_dir %>%
  pivot_longer(cols = c(Increased_variance, Decreased_variance),
               names_to = "Direction",
               values_to = "Count") %>%
  mutate(Direction = gsub("_variance", "", Direction))

p_var_dir <- ggplot(var_dir_long, aes(x = Tissue, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Increased" = "#D6604D", "Decreased" = "#4393C3")) +
  labs(
    title = "Variance Direction in AMD",
    subtitle = "Genes with p < 0.01",
    x = NULL,
    y = "Number of Genes"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(results_dir, "figures/variance_direction.pdf"), p_var_dir,
       width = 6, height = 5)
ggsave(file.path(results_dir, "figures/variance_direction.tiff"), p_var_dir,
       width = 6, height = 5, dpi = 300, compression = "lzw")

cat("  Saved variance direction plot\n")

#%%
cat("\n", rep("=", 70), "\n", sep = "")
cat("Visualizations complete!\n")
cat("All figures saved to:\n")
cat(sprintf("  %s/figures/\n", results_dir))
cat(rep("=", 70), "\n\n", sep = "")
