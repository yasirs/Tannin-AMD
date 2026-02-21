# ==============================================================================
# Part 5: AMD Patient Data Overlay
# ==============================================================================
# Integrate AMD patient differential expression (GSE135092) to validate 
# H2O2/PRG4 findings in human disease context

print("Loading AMD patient data...")
amd_de <- read_csv("results/cohort-GSE135092/GSE135092_DE_results.csv", show_col_types = FALSE)

# Clean and prepare AMD data
amd_de <- amd_de %>%
  filter(!is.na(gene_symbol), gene_symbol != "N/A") %>%
  select(gene_symbol, AMD_logFC = logfc, AMD_pval = pvalue)

# Merge with H2O2/PRG4 data
de_res_amd <- de_res %>%
  left_join(amd_de, by = c("hgnc_symbol" = "gene_symbol"))

# ------------------------------------------------------------------------------
# 5.1: Axis Check Figures with AMD Column
# ------------------------------------------------------------------------------
plot_axis_bubble_amd <- function(axis_name, keys) {
    pattern <- paste(keys, collapse="|")
    subset <- merged_nes %>% 
        filter(grepl(pattern, pathway, ignore.case=TRUE)) %>%
        filter(Padj_H2O2 < 0.1 | Padj_Rescue < 0.1)
    
    if (nrow(subset) == 0) return(NULL)
    
    # Run GSEA on AMD data for these pathways
    amd_ranks <- sign(amd_de$AMD_logFC) * -log10(ifelse(amd_de$AMD_pval==0, 1e-300, amd_de$AMD_pval))
    names(amd_ranks) <- amd_de$gene_symbol
    amd_ranks <- amd_ranks[!is.na(names(amd_ranks))]
    amd_ranks <- sort(amd_ranks, decreasing = TRUE)
    
    # Get pathway subset for AMD GSEA
    pathway_subset <- all_pathways[names(all_pathways) %in% subset$pathway]
    fgsea_amd <- fgseaMultilevel(pathways = pathway_subset, stats = amd_ranks, minSize=15, maxSize=500, eps=0)
    
    subset <- subset %>% 
        left_join(fgsea_amd %>% select(pathway, NES_AMD = NES, Padj_AMD = padj), by = "pathway")
    
    # Top 15
    subset <- subset %>% arrange(desc(abs(NES_H2O2))) %>% head(15)
    
    # Prepare for plot
    dat_long <- subset %>% 
        select(pathway, NES_H2O2, NES_Rescue, NES_AMD) %>%
        pivot_longer(-pathway, names_to = "Condition", values_to = "NES") %>%
        mutate(Condition = case_when(
            Condition == "NES_H2O2" ~ "H2O2 Induction",
            Condition == "NES_Rescue" ~ "PRG4 Rescue",
            Condition == "NES_AMD" ~ "AMD Patients"
        ))
        
    pvals <- subset %>% 
        select(pathway, Padj_H2O2, Padj_Rescue, Padj_AMD) %>%
        pivot_longer(-pathway, names_to = "Condition", values_to = "Padj") %>%
        mutate(Condition = case_when(
            Condition == "Padj_H2O2" ~ "H2O2 Induction",
            Condition == "Padj_Rescue" ~ "PRG4 Rescue",
            Condition == "Padj_AMD" ~ "AMD Patients"
        ))
    
    dat_long$Padj <- pvals$Padj
    dat_long$LogP <- -log10(dat_long$Padj)
    dat_long$LogP[is.na(dat_long$LogP)] <- 0
    dat_long$LogP[is.infinite(dat_long$LogP)] <- 5
    
    p <- ggplot(dat_long, aes(x = Condition, y = pathway, size = LogP, fill = NES)) +
        geom_point(shape = 21, stroke = 0.5) +
        geom_point(data = subset(dat_long, Padj < 0.05), shape = 21, color = "black", stroke = 1) +
        geom_point(data = subset(dat_long, Padj >= 0.05), shape = 21, color = "transparent", stroke = 0) +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
        theme_pub() +
        labs(title = paste(axis_name, "Axis: Model vs Patients"), x = "", y = "") +
        theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(angle = 45, hjust = 1))
        
    ggsave(file.path(out_dir, paste0("Fig_Axis_Check_AMD_", axis_name, ".pdf")), p, width = 8, height = 6)
}

for (ax in names(keywords)) {
    plot_axis_bubble_amd(ax, keywords[[ax]])
}

# ------------------------------------------------------------------------------
# 5.2: Top Reversed Pathways with AMD NES Column
# ------------------------------------------------------------------------------
print("Adding AMD NES to Top Reversed Pathways...")

# Run GSEA on AMD for top reversed pathways
amd_ranks <- sign(amd_de$AMD_logFC) * -log10(ifelse(amd_de$AMD_pval==0, 1e-300, amd_de$AMD_pval))
names(amd_ranks) <- amd_de$gene_symbol
amd_ranks <- amd_ranks[!is.na(names(amd_ranks))]
amd_ranks <- sort(amd_ranks, decreasing = TRUE)

pathway_subset <- all_pathways[names(all_pathways) %in% top_reversed$pathway]
fgsea_amd_top <- fgseaMultilevel(pathways = pathway_subset, stats = amd_ranks, minSize=15, maxSize=500, eps=0)

top_reversed_amd <- top_reversed %>%
    left_join(fgsea_amd_top %>% select(pathway, NES_AMD = NES, Padj_AMD = padj), by = "pathway")

heatmap_dat_amd <- top_reversed_amd %>% select(pathway, NES_H2O2, NES_Rescue, NES_AMD) %>%
  pivot_longer(-pathway, names_to = "Contrast", values_to = "NES")

heatmap_padj_amd <- top_reversed_amd %>% select(pathway, Padj_H2O2, Padj_Rescue, Padj_AMD) %>%
    pivot_longer(-pathway, names_to = "Contrast", values_to = "Padj") %>%
    mutate(Contrast = case_when(
        Contrast == "Padj_H2O2" ~ "NES_H2O2",
        Contrast == "Padj_Rescue" ~ "NES_Rescue",
        Contrast == "Padj_AMD" ~ "NES_AMD"
    ))

heatmap_dat_amd$Padj <- heatmap_padj_amd$Padj[match(paste(heatmap_dat_amd$pathway, heatmap_dat_amd$Contrast), paste(heatmap_padj_amd$pathway, heatmap_padj_amd$Contrast))]

p_heat_amd <- ggplot(heatmap_dat_amd, aes(x = Contrast, y = pathway, fill = NES)) +
  geom_tile(aes(color = ifelse(Padj < 0.05, "black", NA)), size = 0.3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  scale_color_identity() +
  theme_pub() +
  labs(title = "Top Reversed Pathways: Model & Patients", x = "", y = "") +
  theme(axis.text.y = element_text(size = 5))

ggsave(file.path(out_dir, "Fig_Top_Reversed_AMD.pdf"), p_heat_amd, width = 9, height = 8)

# ------------------------------------------------------------------------------
# 5.3: Gene Scatter with AMD Overlay
# ------------------------------------------------------------------------------
print("Creating gene scatter with AMD overlay...")

responsive_genes_amd <- responsive_genes %>%
    left_join(amd_de, by = c("hgnc_symbol" = "gene_symbol")) %>%
    filter(!is.na(AMD_logFC))

# Color by AMD dysregulation magnitude
p_scatter_amd <- ggplot(responsive_genes_amd, aes(x = H2O2PRG4_vs_H2O2.log2FoldChange, y = AMD_logFC, color = abs(AMD_logFC))) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_density_2d(color = "black", alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  scale_color_gradient(low = "lightgray", high = "darkred", name = "AMD |LogFC|") +
  theme_pub() +
  labs(
    title = "PRG4 Rescue vs AMD Dysregulation",
    subtitle = "Genes responsive to H2O2",
    x = "PRG4 Rescue Effect (Log2FC)",
    y = "AMD Patient Dysregulation (LogFC)"
  )
ggsave(file.path(out_dir, "Fig_PRG4_vs_AMD_Scatter.pdf"), p_scatter_amd, width = 7, height = 6)

cor_amd <- cor.test(responsive_genes_amd$H2O2PRG4_vs_H2O2.log2FoldChange, responsive_genes_amd$AMD_logFC)
cat("PRG4 Rescue vs AMD Correlation: r =", cor_amd$estimate, "p =", cor_amd$p.value, "\n")

# ------------------------------------------------------------------------------
# 5.4: Triple Overlap Venn Diagram
# ------------------------------------------------------------------------------
print("Creating triple overlap Venn diagram...")

# Define gene sets
h2o2_induced <- de_res %>% filter(H2O2_vs_CTRL.log2FoldChange > 0.5, H2O2_vs_CTRL.pvalue < 0.05) %>% pull(hgnc_symbol)
prg4_reversed <- de_res %>% filter(H2O2PRG4_vs_H2O2.log2FoldChange < -0.5, H2O2PRG4_vs_H2O2.pvalue < 0.05) %>% pull(hgnc_symbol)
amd_dysreg <- amd_de %>% filter(abs(AMD_logFC) > 0.5, AMD_pval < 0.05) %>% pull(gene_symbol)

# Calculate overlaps
library(VennDiagram)

venn_data <- list(
    "H2O2 Induced" = h2o2_induced,
    "PRG4 Reversed" = prg4_reversed,
    "AMD Dysregulated" = amd_dysreg
)

venn.plot <- venn.diagram(
    x = venn_data,
    filename = NULL,
    category.names = c("H2O2\nInduced", "PRG4\nReversed", "AMD\nDysregulated"),
    fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
    alpha = 0.5,
    cex = 1.5,
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontfamily = "sans"
)

pdf(file.path(out_dir, "Fig_Triple_Overlap_Venn.pdf"), width = 8, height = 8)
grid.draw(venn.plot)
dev.off()

# ------------------------------------------------------------------------------
# 5.5: Clinical Relevance Heatmap
# ------------------------------------------------------------------------------
# Already done in 5.2 above (Top Reversed with AMD column)

# ------------------------------------------------------------------------------
# 5.6: PRG4 Rescue vs AMD Correlation (Detailed)
# ------------------------------------------------------------------------------
# This is essentially 5.3 but let's add a version with significance overlay

responsive_genes_amd <- responsive_genes_amd %>%
    mutate(
        Significant = case_when(
            H2O2PRG4_vs_H2O2.pvalue < 0.05 & AMD_pval < 0.05 ~ "Both Significant",
            H2O2PRG4_vs_H2O2.pvalue < 0.05 ~ "PRG4 Only",
            AMD_pval < 0.05 ~ "AMD Only",
            TRUE ~ "Neither"
        )
    )

p_scatter_amd_sig <- ggplot(responsive_genes_amd, aes(x = H2O2PRG4_vs_H2O2.log2FoldChange, y = AMD_logFC, color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", color = "black", se = TRUE, linetype = "dashed") +
  scale_color_manual(values = c("Both Significant" = "darkred", "PRG4 Only" = "blue", 
                                 "AMD Only" = "green", "Neither" = "lightgray")) +
  theme_pub() +
  labs(
    title = "PRG4 Rescue vs AMD Dysregulation (Significance Overlay)",
    x = "PRG4 Rescue Effect (Log2FC)",
    y = "AMD Patient Dysregulation (LogFC)"
  )
ggsave(file.path(out_dir, "Fig_PRG4_vs_AMD_Significance.pdf"), p_scatter_amd_sig, width = 8, height = 6)

print("AMD overlay analysis complete!")
