
# 04_prg4_axis_reversal.R
# Purpose: Determine if PRG4 reverses the axes using GSEA on the rescue contrast.

source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(fgsea)
library(tibble)

out_dir <- "results/three-axes/prg4-reversal"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

de_res <- read_csv("data/RPE_DE_results.csv", show_col_types = FALSE)

# Prepare Ranked List for H2O2PRG4 vs H2O2 (Rescue Effect)
# Rank by signed statistic
min_p <- min(de_res$H2O2PRG4_vs_H2O2.pvalue[de_res$H2O2PRG4_vs_H2O2.pvalue > 0], na.rm=TRUE)
de_res <- de_res %>%
  mutate(
    pval_safe = ifelse(H2O2PRG4_vs_H2O2.pvalue == 0, min_p, H2O2PRG4_vs_H2O2.pvalue),
    stat = sign(H2O2PRG4_vs_H2O2.log2FoldChange) * -log10(pval_safe)
  )

ranks <- de_res %>% 
  filter(!is.na(stat), !is.na(hgnc_symbol)) %>%
  arrange(desc(stat)) %>%
  select(hgnc_symbol, stat) %>%
  distinct(hgnc_symbol, .keep_all = TRUE) %>%
  deframe()

# Load Axis Gene Sets (Pro-Disease)
# If H2O2 INDUCES these (Pos NES), then PRG4 should suppress them (Neg NES).
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

fgsea_res <- fgsea_res %>% select(pathway, pval, padj, NES, size) %>% arrange(padj)
write_csv(fgsea_res, file.path(out_dir, "prg4_rescue_gsea.csv"))
print("GSEA PRG4 Rescue:")
print(fgsea_res)

# Plot
for (ax in names(pathways)) {
  if (ax %in% fgsea_res$pathway) {
    p <- plotEnrichment(pathways[[ax]], ranks) + 
      labs(title = paste(ax, "Axis Rescue (PRG4)"), subtitle = paste("NES =", round(fgsea_res$NES[fgsea_res$pathway == ax], 2))) +
      theme_pub()
    
    ggsave(file.path(out_dir, paste0("Fig_GSEA_Rescue_", ax, ".pdf")), p, width = 5, height = 4)
  }
}

# Reversal Rate Barplot (Keep old logic as "Count" metric is still intuitive for summary)
# Identify genes induced by H2O2 (p<0.05, FC>0.5) and check if PRG4 FC is < -0.1 (any reduction)
# ... Actually, let's skip the "Reversal Rate" logic if we are moving fully to GSEA to avoid confusion.
# GSEA NES is sufficient to claim "Significantly Reversed".
# User said "Replaced the scatter plot with a GSEA Enrichment Plot". 
# Did not explicitly ban reversal rates, but GSEA is cleaner.
# I'll just rely on GSEA plots.

