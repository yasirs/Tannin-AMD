library(ggplot2)
library(extrafont)

# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
COV_DIR <- file.path(BASE_DIR, "results/coverage-analysis")
OUTPUT_DIR <- COV_DIR

# Custom Theme
theme_tannin <- function() {
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )
}

# 1. Coverage Bar Chart
df_cov <- read.csv(file.path(COV_DIR, "coverage_summary.csv"))

# Reshape for ggplot
df_long <- data.frame(
  Contrast = rep(df_cov$Contrast, 3),
  Dataset = c(rep("RPE1 Essential", nrow(df_cov)), rep("K562 Essential", nrow(df_cov)), rep("K562 GWPS", nrow(df_cov))),
  Pct_Covered = c(df_cov$Pct_KD_RPE1_Essential, df_cov$Pct_KD_K562_Essential, df_cov$Pct_KD_K562_GWPS)
)

# Order factors
df_long$Dataset <- factor(df_long$Dataset, levels = c("RPE1 Essential", "K562 Essential", "K562 GWPS"))

p1 <- ggplot(df_long, aes(x = Contrast, y = Pct_Covered, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("RPE1 Essential" = "#A6BDDB", "K562 Essential" = "#74A9CF", "K562 GWPS" = "#0570B0")) +
  theme_tannin() +
  labs(
    title = "Perturb-seq Coverage of Bulk DEGs",
    x = "",
    y = "% of DEGs Perturbed"
  ) +
  coord_flip()

ggsave(file.path(OUTPUT_DIR, "coverage_comparison.pdf"), p1, width = 7, height = 4)
ggsave(file.path(OUTPUT_DIR, "coverage_comparison.tiff"), p1, width = 7, height = 4, device = "tiff", compression = "lzw", dpi = 300)

# 2. Enrichment of Missing Genes (Example: H2O2 missing from RPE1)
# Check if file exists
enr_file <- file.path(COV_DIR, "enrichment_missing_H2O2_Stress_RPE1_Essential.csv")
if (file.exists(enr_file)) {
  df_enr <- read.csv(enr_file)
  if (nrow(df_enr) > 0) {
    df_enr <- head(df_enr[order(df_enr$pvalue), ], 10)
    df_enr$pathway <- factor(df_enr$pathway, levels = rev(df_enr$pathway))
    
    p2 <- ggplot(df_enr, aes(x = pathway, y = -log10(pvalue))) +
      geom_bar(stat = "identity", fill = "#E34A33") +
      coord_flip() +
      theme_tannin() +
      labs(
        title = "Pathways Missing from RPE1 Essential\n(H2O2 Stress Signature)",
        x = "",
        y = "-Log10 P-value"
      )
    
    ggsave(file.path(OUTPUT_DIR, "enrichment_missing_h2o2.pdf"), p2, width = 8, height = 5)
    ggsave(file.path(OUTPUT_DIR, "enrichment_missing_h2o2.tiff"), p2, width = 8, height = 5, device = "tiff", compression = "lzw", dpi = 300)
  }
}

print("Visualizations complete.")
