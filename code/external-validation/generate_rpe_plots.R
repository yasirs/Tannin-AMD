# Quick script to generate comparison plots from existing results
library(ggplot2)
library(dplyr)
library(readr)

base_dir <- "/home/ysuhail/work/Tannin-AMD"
out_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")

# Load results
results_rpe_mac <- read_csv(file.path(out_dir, "RPE_Macula_all_results.csv"))
results_rpe_nonmac <- read_csv(file.path(out_dir, "RPE_nonMacula_all_results.csv"))

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

# Count function
count_sig_genes <- function(results_df, tissue_label) {
  data.frame(
    Tissue = rep(tissue_label, 6),
    Method = c(
      "Mean DE (limma)", "Mean DE (limma)",
      "Variance (Levene)", "Variance (Levene)",
      "Residual Variance", "Residual Variance"
    ),
    Threshold = rep(c("p < 0.01", "FDR < 0.05"), 3),
    Count = c(
      sum(results_df$P.Value < 0.01, na.rm = TRUE),
      sum(results_df$adj.P.Val < 0.05, na.rm = TRUE),
      sum(results_df$levene_pval < 0.01, na.rm = TRUE),
      sum(results_df$levene_fdr < 0.05, na.rm = TRUE),
      sum(results_df$resid_var_pval < 0.01, na.rm = TRUE),
      sum(results_df$resid_var_fdr < 0.05, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

counts_mac <- count_sig_genes(results_rpe_mac, "RPE Macula")
counts_nonmac <- count_sig_genes(results_rpe_nonmac, "RPE non-Macula")
all_counts <- rbind(counts_mac, counts_nonmac)

write_csv(all_counts, file.path(out_dir, "gene_counts_summary.csv"))

# Create plot
p_counts <- ggplot(all_counts, aes(x = Method, y = Count, fill = Threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Tissue, scales = "free_y") +
  scale_fill_manual(values = c("p < 0.01" = "#FC8D62", "FDR < 0.05" = "#8DA0CB")) +
  labs(
    title = "Significant Genes Detected by Method and Threshold",
    subtitle = "RPE samples with age covariate",
    x = NULL,
    y = "Number of Significant Genes"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "Fig_gene_counts_comparison.pdf"), p_counts, width = 8, height = 5)
ggsave(file.path(out_dir, "Fig_gene_counts_comparison.tiff"), p_counts, 
       width = 8, height = 5, dpi = 300, compression = "lzw")

# Print summary
cat("\nSUMMARY OF SIGNIFICANT GENES\n")
cat(rep("=", 70), "\n", sep = "")
print(all_counts)
cat("\nPlots saved!\n")
