
# 03_h2o2_axis_induction.R
# Purpose: Determine if H2O2 induces the three AMD axes using GSEA.

source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(fgsea)
library(tibble)

# Constants
out_dir <- "results/three-axes/h2o2-induction"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load DE Results
de_res <- read_csv("data/RPE_DE_results.csv", show_col_types = FALSE)

# Prepare Ranked List for H2O2 vs CTRL
# Rank by signed statistic: sign(logFC) * -log10(pvalue)
# Handle p=0 by using min non-zero or just large number
min_p <- min(de_res$H2O2_vs_CTRL.pvalue[de_res$H2O2_vs_CTRL.pvalue > 0], na.rm=TRUE)
de_res <- de_res %>%
  mutate(
    pval_safe = ifelse(H2O2_vs_CTRL.pvalue == 0, min_p, H2O2_vs_CTRL.pvalue),
    stat = sign(H2O2_vs_CTRL.log2FoldChange) * -log10(pval_safe)
  )

# Remove NAs and sort
ranks <- de_res %>% 
  filter(!is.na(stat), !is.na(hgnc_symbol)) %>%
  arrange(desc(stat)) %>%
  select(hgnc_symbol, stat) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>% # Handle dups
  deframe()

# Load Axis Gene Sets (Pro-Disease)
gene_sets_dir <- "results/three-axes/gene-sets"
load_set <- function(name) {
  f <- file.path(gene_sets_dir, paste0(name, "_mapped.csv"))
  if (file.exists(f)) read_csv(f, show_col_types = FALSE)$gene_symbol else NULL
}

pathways <- list(
  Oxidative = load_set("oxidative_pro_disease"),
  Inflammatory = load_set("inflammatory_pro_disease"),
  Senescence = load_set("senescence_pro_disease")
)

# Run GSEA
set.seed(42)
fgsea_res <- fgseaMultilevel(pathways = pathways, stats = ranks, minSize=15, maxSize=500)

# Save Results
fgsea_res <- fgsea_res %>% select(pathway, pval, padj, NES, size) %>% arrange(padj)
write_csv(fgsea_res, file.path(out_dir, "h2o2_axis_gsea.csv"))
print("GSEA H2O2 Induction:")
print(fgsea_res)

# Plot Enrichment (One per axis)
for (ax in names(pathways)) {
  if (ax %in% fgsea_res$pathway) {
    p <- plotEnrichment(pathways[[ax]], ranks) + 
      labs(title = paste(ax, "Axis Induction (H2O2)"), subtitle = paste("NES =", round(fgsea_res$NES[fgsea_res$pathway == ax], 2))) +
      theme_pub()
    
    ggsave(file.path(out_dir, paste0("Fig_GSEA_", ax, ".pdf")), p, width = 5, height = 4)
  }
}

# Also plot Canonical Markers (Kept from old script as it's useful verification)
# Just a simple bar plot of select genes
markers <- c("NQO1", "HMOX1", "IL6", "CXCL8", "CDKN1A", "CDKN2A")
marker_dat <- de_res %>% 
  filter(hgnc_symbol %in% markers) %>%
  select(hgnc_symbol, H2O2_vs_CTRL.log2FoldChange, H2O2_vs_CTRL.pvalue)

p_mark <- ggplot(marker_dat, aes(x = hgnc_symbol, y = H2O2_vs_CTRL.log2FoldChange)) +
  geom_col(fill = "grey") +
  theme_pub() +
  labs(title = "Canonical Marker Induction", y = "Log2FC (H2O2 vs CTRL)") +
  geom_text(aes(label = ifelse(H2O2_vs_CTRL.pvalue < 0.05, "*", "ns")), vjust = -0.5)

ggsave(file.path(out_dir, "Fig_H2O2_canonical_markers.pdf"), p_mark, width = 6, height = 4)
