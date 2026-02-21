
source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
# library(fgsea) # Plan mentioned fgsea. Check if installed. If not, Fisher test (hypergeometric) is already done in 03/04.
# Use enrichment barplot logic based on hypergeometric results for simplicity if fgsea unavailable?
# Or try simple Wilcoxon rank sum (GSEA-like) if fgsea missing.
# Let's try to load fgsea or clusterProfiler. If missing, manual GSEA or enrichment only.
# "perform enrichment analysis 07_axis_enrichment.R"

# Constants
FDR_THRESH <- 0.05
LFC_THRESH <- 0.5
de_file <- "data/RPE_DE_results.csv"
gene_sets_dir <- "results/three-axes/gene-sets"
out_dir <- "results/three-axes/enrichment"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

de_res <- read_csv(de_file)

load_genes <- function(name) {
  f <- file.path(gene_sets_dir, paste0(name, "_mapped.csv"))
  if (file.exists(f)) read_csv(f, show_col_types = FALSE)$gene_symbol else c()
}

axes <- list(
  Oxidative = load_genes("oxidative_pro_disease"),
  Inflammatory = load_genes("inflammatory_pro_disease"),
  Senescence = load_genes("senescence_pro_disease")
)

# Define PRG4 Rescue Signature
# H2O2 vs CTRL UP (Disease) AND H2O2PRG4 vs H2O2 DOWN (Rescue)
# But simplified: Just rank genes by PRG4 Rescue Effect (H2O2PRG4 vs H2O2 log2FC)
# Expected result: Pro-disease genes should be at the BOTTOM (Downregulated by PRG4).
# So we test for enrichment in BOTTOM of list (Negative enrichment score).

# Rank genes
vals <- de_res$H2O2PRG4_vs_H2O2.log2FoldChange
names(vals) <- de_res$hgnc_symbol
vals <- sort(vals, decreasing = TRUE)
vals <- vals[!is.na(vals)]

# Perform Enrichment
results <- data.frame()

for (ax in names(axes)) {
    genes <- axes[[ax]]
    # Filter genes present in ranking
    genes <- intersect(genes, names(vals))
    
    # 1. Wilcoxon Rank Sum Test (Simple GSEA equivalent)
    # Test if genes in set have lower ranks (more negative LFC) than background
    # Use "greater" if looking for UP, "less" if looking for DOWN.
    # We expect rescue -> DOWN. So genes should have lower values.
    # Alternative = "less"
    
    bg_vals <- vals[!names(vals) %in% genes]
    test_res <- wilcox.test(vals[genes], bg_vals, alternative = "less")
    
    results <- rbind(results, data.frame(
        Axis = ax,
        PValue = test_res$p.value,
        MedianLFC_Set = median(vals[genes]),
        MedianLFC_Bg = median(bg_vals)
    ))
}

results$FDR <- p.adjust(results$PValue, method = "BH")
results$LogP <- -log10(results$PValue)

print(results)
write.csv(results, file.path(out_dir, "enrichment_results.csv"), row.names = FALSE)

# Plot
p <- ggplot(results, aes(x = Axis, y = LogP, fill = Axis)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed") +
  theme_pub() +
  labs(title = "Enrichment of Axes in PRG4 Rescue Signature", y = "-log10(P-value) [Wilcoxon Less]", subtitle = "Test: Axis genes are downregulated by PRG4")

ggsave(file.path(out_dir, "Fig_Enrichment_barplot.pdf"), p, width = 4, height = 4)
