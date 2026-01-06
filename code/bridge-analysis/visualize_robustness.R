library(ggplot2)
library(extrafont)
library(openxlsx)
library(data.table)

# Load fonts if necessary
# extrafont::font_import(prompt = FALSE) # This takes too long, assuming standard fonts or will fall back gracefully
extrafont::loadfonts(device = "postscript", quiet = TRUE)
extrafont::loadfonts(device = "pdf", quiet = TRUE)

# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
ROBUST_DIR <- file.path(BASE_DIR, "results/robustness-analysis")
OUTPUT_DIR <- ROBUST_DIR

# Custom Theme
theme_tannin <- function() {
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )
}

# 1. Visualize Threshold Comparison
df_comp <- read.csv(file.path(ROBUST_DIR, "threshold_comparison.csv"))

# Filter out Top N as they are not reliable
df_comp_filtered <- df_comp[!grepl("Top", df_comp$threshold), ]

p1 <- ggplot(df_comp_filtered, aes(x = threshold, y = top_correlation, color = contrast, group = contrast)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_tannin() +
  labs(
    title = "Bridge Correlation Stability",
    x = "Threshold",
    y = "Top Spearman Rho",
    color = "Contrast"
  ) +
  scale_color_manual(values = c("H2O2_Stress" = "#A6611A", "PRG4_Baseline" = "#018571", "PRG4_Rescue" = "#DFC27D"))

ggsave(file.path(OUTPUT_DIR, "threshold_sensitivity.pdf"), p1, width = 6, height = 4)
ggsave(file.path(OUTPUT_DIR, "threshold_sensitivity.tiff"), p1, width = 6, height = 4, device = "tiff", compression = "lzw", dpi = 300)

# 2. Visualize Stability Metrics (Jaccard)
df_stab <- read.csv(file.path(ROBUST_DIR, "stability_metrics.csv"))

# Focus on standard thresholds
df_stab_sub <- df_stab[!grepl("Top", df_stab$threshold_1) & !grepl("Top", df_stab$threshold_2), ]

p2 <- ggplot(df_stab_sub, aes(x = interaction(threshold_1, threshold_2, sep = " vs "), y = jaccard_top50, fill = contrast)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_tannin() +
  coord_flip() +
  labs(
    title = "Top 50 KD List Stability",
    x = "Threshold Comparison",
    y = "Jaccard Similarity",
    fill = "Contrast"
  ) +
  scale_fill_manual(values = c("H2O2_Stress" = "#A6611A", "PRG4_Baseline" = "#018571", "PRG4_Rescue" = "#DFC27D"))

ggsave(file.path(OUTPUT_DIR, "jaccard_stability.pdf"), p2, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "jaccard_stability.tiff"), p2, width = 7, height = 5, device = "tiff", compression = "lzw", dpi = 300)

# 3. Enrichment Summary
df_enr <- read.csv(file.path(ROBUST_DIR, "signature_enrichment_summary.csv"))

# Plot top 10 for H2O2 p < 0.01 (the recommended one)
df_enr_sub <- df_enr[df_enr$signature == "H2O2_p01", ]
df_enr_sub <- head(df_enr_sub[order(df_enr_sub$pvalue), ], 10)
df_enr_sub$pathway <- factor(df_enr_sub$pathway, levels = rev(df_enr_sub$pathway))

p3 <- ggplot(df_enr_sub, aes(x = pathway, y = -log10(pvalue))) +
  geom_bar(stat = "identity", fill = "#8C510A") +
  coord_flip() +
  theme_tannin() +
  labs(
    title = "Top Enrichment: H2O2 Stress (p < 0.01)",
    x = "",
    y = "-Log10 P-value"
  )

ggsave(file.path(OUTPUT_DIR, "enrichment_h2o2_p01.pdf"), p3, width = 8, height = 5)
ggsave(file.path(OUTPUT_DIR, "enrichment_h2o2_p01.tiff"), p3, width = 8, height = 5, device = "tiff", compression = "lzw", dpi = 300)

print("Visualizations complete.")
