library(ggplot2)
library(grid)
library(gridExtra)

# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_FILE <- file.path(BASE_DIR, "results/external-validation/gse129964_comparison.csv")
OUTPUT_DIR <- file.path(BASE_DIR, "results/gsea-analysis")

# Load Data
df <- read.csv(DATA_FILE)

# 1. Prepare Ranked List
# Rank by Serum Starvation effect
df <- df[order(df$logfc_serum_stress, decreasing = TRUE), ]
gene_ranks <- df$logfc_serum_stress
names(gene_ranks) <- df$gene_symbol

# 2. Define Gene Sets
# PRG4 Rescue UP (p < 0.05, FC > 0)
set_up <- df$gene_symbol[df$H2O2PRG4_vs_H2O2.pvalue < 0.05 & df$H2O2PRG4_vs_H2O2.log2FoldChange > 0]
set_up <- set_up[!is.na(set_up)]

# PRG4 Rescue DOWN (p < 0.05, FC < 0)
set_down <- df$gene_symbol[df$H2O2PRG4_vs_H2O2.pvalue < 0.05 & df$H2O2PRG4_vs_H2O2.log2FoldChange < 0]
set_down <- set_down[!is.na(set_down)]

# Function to calculate Running ES
calc_running_es <- function(ranked_stats, gene_set) {
  N <- length(ranked_stats)
  Nh <- length(gene_set)
  Nm <- N - Nh
  
  # Hit vector
  hits <- names(ranked_stats) %in% gene_set
  
  # Weights (standard GSEA uses weight=1 for correlation, or 0 for classic KS)
  # fgsea default is 1.0
  scores <- abs(ranked_stats)
  scores[!hits] <- 0
  
  P_hit <- cumsum(scores) / sum(scores)
  P_miss <- cumsum(!hits) / Nm
  
  RES <- P_hit - P_miss
  return(data.frame(x = 1:N, y = RES, hit = hits))
}

plot_gsea <- function(res_df, title, filename) {
  
  # Find ES (max deviation)
  max_dev_idx <- which.max(abs(res_df$y))
  es <- res_df$y[max_dev_idx]
  
  # Plot Curve
  p_curve <- ggplot(res_df, aes(x = x, y = y)) +
    geom_line(color = "#4DAF4A", linewidth = 1.5) + # Green standard
    geom_hline(yintercept = 0, color = "black") +
    geom_vline(xintercept = max_dev_idx, linetype = "dashed", color = "grey50") +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(b = 0),
      axis.text = element_text(color = "black")
    ) +
    labs(
      title = title,
      subtitle = paste0("Enrichment Score: ", round(es, 3)),
      y = "Enrichment Score"
    )
  
  # Plot Rug (Barcode)
  p_rug <- ggplot(res_df[res_df$hit, ], aes(x = x, y = 1)) +
    geom_errorbar(aes(ymin = 0, ymax = 1), color = "black", width = 0, linewidth = 0.2) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(limits = c(1, nrow(res_df)), expand = c(0, 0)) +
    theme_void() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(t = 0)
    ) +
    labs(x = "Rank in Ordered Dataset (Serum Starvation)")
  
  # Arrange
  g <- arrangeGrob(p_curve, p_rug, heights = c(3, 0.5))
  ggsave(filename, g, width = 6, height = 4)
}

# Run
print("Calculating ES for Rescue UP...")
res_up <- calc_running_es(gene_ranks, set_up)
plot_gsea(res_up, "GSEA: PRG4 Rescue UP vs Serum Starvation", file.path(OUTPUT_DIR, "gsea_standard_rescue_up.pdf"))

print("Calculating ES for Rescue DOWN...")
res_down <- calc_running_es(gene_ranks, set_down)
plot_gsea(res_down, "GSEA: PRG4 Rescue DOWN vs Serum Starvation", file.path(OUTPUT_DIR, "gsea_standard_rescue_down.pdf"))

print("Plots saved.")
