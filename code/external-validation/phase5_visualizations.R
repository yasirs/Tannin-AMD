# Phase 5: Comprehensive Visualization
# Creates faceted barplot and pathway enrichment panels

#%%
library(dplyr)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(gridExtra)
library(tidyr)

#%%
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
OUTPUT_DIR <- file.path(BASE_DIR, "results/external-validation/gse29801_prg4_comparison")
GENE_SET_DIR <- file.path(OUTPUT_DIR, "gene_sets")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")

FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

#%%
# ========================================
# FIGURE 1: FACETED BARPLOT
# ========================================

message("=== Creating Figure 1: Faceted Barplot ===\n")

# Load GSEA results
gsea_results <- read.csv(file.path(TABLE_DIR, "gsea_summary.csv"))

# Prepare data for plotting
plot_data <- gsea_results %>%
  mutate(
    gene_set_label = factor(gene_set, 
                           levels = c("de", "low_var", "high_var", "bimodal", "age_de"),
                           labels = c("DE (p<0.01)", "Low Variance", "High Variance", "Bimodal", "Age DE")),
    comparison_label = factor(comparison,
                             levels = c("H2O2_vs_Control", "PRG4_Rescue"),
                             labels = c("H2O2 vs Control", "PRG4 Rescue")),
    tissue_label = factor(tissue, levels = c("Macular", "Extramacular")),
    signif = case_when(
      pval < 0.05 ~ "p < 0.05",
      pval < 0.20 ~ "p < 0.20",
      TRUE ~ "NS"
    ),
    signif = factor(signif, levels = c("p < 0.05", "p < 0.20", "NS"))
  )

# Create faceted barplot
p1 <- ggplot(plot_data, aes(x = gene_set_label, y = NES, fill = signif)) +
  geom_col(position = "dodge", color = "black", size = 0.3) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  facet_grid(tissue_label ~ comparison_label, scales = "free_x") +
  scale_fill_manual(
    values = c("p < 0.05" = "#d73027", "p < 0.20" = "#fc8d59", "NS" = "gray70"),
    name = "Significance"
  ) +
  labs(
    title = "GSE29801 Gene Set Enrichment in PRG4 Data",
    subtitle = "Normalized Enrichment Score (NES) by Gene Set Type",
    x = "Gene Set",
    y = "NES"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top",
    panel.grid.major.y = element_line(color = "gray90", size = 0.3)
  ) +
  coord_cartesian(ylim = c(-1.5, 1.5))

ggsave(file.path(FIGURE_DIR, "figure1_gene_set_enrichment_barplot.pdf"),
       p1, width = 12, height = 8)
ggsave(file.path(FIGURE_DIR, "figure1_gene_set_enrichment_barplot.tiff"),
       p1, width = 12, height = 8, compression = "lzw", dpi = 300)

message("✓ Figure 1 saved\n")

#%%
# ========================================
# FIGURE 2: PATHWAY ENRICHMENT FOR GENE SETS
# ========================================

message("=== Creating Figure 2: Pathway Enrichment Panels ===\n")

# Load all gene sets (simpler approach - read from original results)
DE_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")

gene_sets <- list(
  macular_de = read.csv(file.path(DE_DIR, "macular_variance_analysis_annotated.csv")) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "" & f_pvalue < 0.01) %>% 
    pull(gene_symbol) %>% unique(),
  
  macular_low = read.csv(file.path(GENE_SET_DIR, "macular_low_variance_age_adj.csv")) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>% pull(gene_symbol) %>% unique(),
  
  macular_high = read.csv(file.path(GENE_SET_DIR, "macular_high_variance_age_adj.csv")) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>% pull(gene_symbol) %>% unique(),
  
  extramacular_de = read.csv(file.path(DE_DIR, "extramacular_variance_analysis_annotated.csv")) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "" & f_pvalue < 0.01) %>% 
    pull(gene_symbol) %>% unique(),
  
  extramacular_low = read.csv(file.path(GENE_SET_DIR, "extramacular_low_variance_age_adj.csv")) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>% pull(gene_symbol) %>% unique(),
  
  extramacular_high = read.csv(file.path(GENE_SET_DIR, "extramacular_high_variance_age_adj.csv")) %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>% pull(gene_symbol) %>% unique()
)

message(sprintf("Loaded %d gene sets", length(gene_sets)))
for (name in names(gene_sets)) {
  message(sprintf("  %s: %d genes", name, length(gene_sets[[name]])))
}

# Load MSigDB gene sets
message("Loading MSigDB collections...")
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
go_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
curated <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP")

# Filter GO BP to reasonable size
go_bp_filtered <- go_bp %>%
  group_by(gs_name) %>%
  summarise(size = n()) %>%
  filter(size >= 15 & size <= 500) %>%
  pull(gs_name)

go_bp <- go_bp %>% filter(gs_name %in% go_bp_filtered)

# Create pathway lists
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
go_bp_list <- split(go_bp$gene_symbol, go_bp$gs_name)
curated_list <- split(curated$gene_symbol, curated$gs_name)

# Function to run ORA for one gene set
run_ora_for_geneset <- function(gene_set, gene_set_name, pathway_collection, collection_name, background_size = 20000) {
  
  if (length(gene_set) < 5) {
    return(NULL)
  }
  
  results_list <- list()
  
  for (pathway_name in names(pathway_collection)) {
    pathway_genes <- pathway_collection[[pathway_name]]
    
    # 2x2 contingency table
    in_both <- sum(gene_set %in% pathway_genes)
    in_geneset_not_pathway <- length(gene_set) - in_both
    in_pathway_not_geneset <- length(pathway_genes) - in_both
    in_neither <- background_size - in_geneset_not_pathway - in_pathway_not_geneset - in_both
    
    if (in_both == 0) next
    
    # Fisher's exact test
    fisher_result <- fisher.test(matrix(c(in_both, in_geneset_not_pathway,
                                          in_pathway_not_geneset, in_neither),
                                        nrow = 2))
    
    results_list[[pathway_name]] <- data.frame(
      pathway = pathway_name,
      overlap = in_both,
      pathway_size = length(pathway_genes),
      gene_set_size = length(gene_set),
      pvalue = fisher_result$p.value,
      odds_ratio = as.numeric(fisher_result$estimate),
      collection = collection_name,
      gene_set = gene_set_name
    )
  }
  
  if (length(results_list) == 0) return(NULL)
  
  results_df <- bind_rows(results_list)
  results_df$fdr <- p.adjust(results_df$pvalue, method = "BH")
  
  return(results_df)
}

# Run ORA for key gene sets (focus on those with enough genes)
message("\nRunning ORA for gene sets...")

ora_results_list <- list()

# Extramacular high variance (largest set)
ora_results_list[["extramacular_high_hallmark"]] <- run_ora_for_geneset(
  gene_sets$extramacular_high, "Extramacular High-Var", hallmark_list, "Hallmark"
)
ora_results_list[["extramacular_high_gobp"]] <- run_ora_for_geneset(
  gene_sets$extramacular_high, "Extramacular High-Var", go_bp_list, "GO BP"
)

# Extramacular DE
ora_results_list[["extramacular_de_hallmark"]] <- run_ora_for_geneset(
  gene_sets$extramacular_de, "Extramacular DE", hallmark_list, "Hallmark"
)
ora_results_list[["extramacular_de_gobp"]] <- run_ora_for_geneset(
  gene_sets$extramacular_de, "Extramacular DE", go_bp_list, "GO BP"
)

# Extramacular low variance
ora_results_list[["extramacular_low_hallmark"]] <- run_ora_for_geneset(
  gene_sets$extramacular_low, "Extramacular Low-Var", hallmark_list, "Hallmark"
)

# Macular DE
ora_results_list[["macular_de_hallmark"]] <- run_ora_for_geneset(
  gene_sets$macular_de, "Macular DE", hallmark_list, "Hallmark"
)

# Combine all ORA results
ora_combined <- bind_rows(ora_results_list[!sapply(ora_results_list, is.null)])

# Get top pathways for each gene set/collection
top_pathways <- ora_combined %>%
  group_by(gene_set, collection) %>%
  arrange(pvalue) %>%
  slice_head(n = 5) %>%
  ungroup()

# Create multi-panel dotplot
p2 <- ggplot(top_pathways, aes(x = gene_set, y = reorder(pathway, -pvalue), 
                               size = overlap, color = -log10(pvalue))) +
  geom_point() +
  facet_wrap(~ collection, scales = "free_y", ncol = 2) +
  scale_color_gradient(low = "gray70", high = "#d73027", name = "-log10(p)") +
  scale_size_continuous(range = c(2, 8), name = "Overlap") +
  labs(
    title = "Top Enriched Pathways per GSE29801 Gene Set",
    subtitle = "Over-Representation Analysis (Fisher's Exact Test)",
    x = "Gene Set",
    y = "Pathway"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right"
  )

ggsave(file.path(FIGURE_DIR, "figure2_pathway_enrichment_panels.pdf"),
       p2, width = 14, height = 12)
ggsave(file.path(FIGURE_DIR, "figure2_pathway_enrichment_panels.tiff"),
       p2, width = 14, height = 12, compression = "lzw", dpi = 300)

# Save ORA results
write.csv(ora_combined,
          file.path(TABLE_DIR, "pathway_enrichment_ora_results.csv"),
          row.names = FALSE)

message("✓ Figure 2 saved\n")

message("\n=== Phase 5 Complete ===")
message("Generated:")
message("  - figure1_gene_set_enrichment_barplot.pdf/.tiff")
message("  - figure2_pathway_enrichment_panels.pdf/.tiff")
message("  - pathway_enrichment_ora_results.csv")
