library(ggplot2)
library(extrafont)

# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
RES_DIR <- file.path(BASE_DIR, "results/concordance-analysis")
MAP_FILE <- file.path(BASE_DIR, "results/baseline-expression/all_baseline_expression.csv")
OUTPUT_DIR <- RES_DIR

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

# Load Data
df_res <- read.csv(file.path(RES_DIR, "gene_concordance.csv"))
df_map <- read.csv(MAP_FILE)

# Map Symbols
# df_map has ensembl_gene_id, gene_name
df_res$gene_symbol <- df_map$gene_name[match(df_res$target_ensembl_id, df_map$ensembl_gene_id)]

# 1. Distribution Plot
p1 <- ggplot(df_res, aes(x = concordance_r)) +
  geom_histogram(binwidth = 0.05, fill = "#2B8CBE", color = "white") +
  geom_vline(xintercept = median(df_res$concordance_r), linetype = "dashed", color = "red") +
  theme_tannin() +
  labs(
    title = "Distribution of RPE1-K562 Concordance",
    subtitle = paste0("Median Pearson r = ", round(median(df_res$concordance_r), 3)),
    x = "Pearson Correlation (r)",
    y = "Count of Knockdowns"
  )

ggsave(file.path(OUTPUT_DIR, "concordance_distribution.pdf"), p1, width = 6, height = 4)

# 2. Top Genes Table (visualized as text/bar? No, just top 20 bar)
top20 <- head(df_res[order(-df_res$concordance_r), ], 20)
# Fix factor order
top20$gene_symbol <- factor(top20$gene_symbol, levels = rev(top20$gene_symbol))

p2 <- ggplot(top20, aes(x = gene_symbol, y = concordance_r)) +
  geom_bar(stat = "identity", fill = "#E34A33") +
  coord_flip() +
  theme_tannin() +
  labs(
    title = "Top 20 Most Concordant Genes",
    x = "",
    y = "Pearson Correlation"
  )

ggsave(file.path(OUTPUT_DIR, "top_concordant_genes.pdf"), p2, width = 6, height = 5)

print("Visualizations complete.")
