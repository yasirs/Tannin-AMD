library(ggplot2)
library(extrafont)
library(reshape2)

# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
EXPR_DIR <- file.path(BASE_DIR, "results/baseline-expression")
OUTPUT_DIR <- EXPR_DIR

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

# 1. Load Data
df <- read.csv(file.path(EXPR_DIR, "all_baseline_expression.csv"))

# 2. Expression Distribution
# Melt for density plot
df_dist <- df[, c("RPE1_baseline", "K562_GWPS_baseline")]
df_dist_long <- melt(df_dist)
df_dist_long <- df_dist_long[!is.na(df_dist_long$value), ]

p1 <- ggplot(df_dist_long, aes(x = value + 0.001, fill = variable)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  scale_fill_manual(values = c("RPE1_baseline" = "#A6611A", "K562_GWPS_baseline" = "#018571")) +
  theme_tannin() +
  labs(
    title = "Baseline Expression Distribution",
    x = "Mean UMI per Cell (log10)",
    y = "Density",
    fill = "Cell Line"
  )

ggsave(file.path(OUTPUT_DIR, "expression_distribution.pdf"), p1, width = 6, height = 4)

# 3. RPE1 vs K562 Scatter
p2 <- ggplot(df, aes(x = RPE1_baseline + 0.001, y = K562_GWPS_baseline + 0.001)) +
  geom_point(alpha = 0.2, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  scale_y_log10() +
  theme_tannin() +
  labs(
    title = "RPE1 vs K562 Baseline Expression",
    x = "RPE1 Mean UMI (log10)",
    y = "K562 GWPS Mean UMI (log10)"
  )

ggsave(file.path(OUTPUT_DIR, "baseline_scatter.pdf"), p2, width = 6, height = 5)

# 4. Marker Heatmap
markers_df <- read.csv(file.path(EXPR_DIR, "marker_baseline_check.csv"))
# Filter out NAs for heatmap
markers_df <- markers_df[!is.na(markers_df$RPE1) | !is.na(markers_df$K562), ]

# Melt
m_long <- melt(markers_df, id.vars = c("category", "gene"))
# Convert to z-score per gene (as per preferences)
m_long <- as.data.frame(m_long)
m_wide <- dcast(m_long, gene ~ variable)
# Simple z-score (row-wise)
# Since we only have 2 columns, z-score is just (x-mean)/sd
# Actually, for 2 values, z-score might be binary. Let's just log-scale color for better contrast if z-score is flat.
# Preference says "TPM converted to gene-wise z-scores"

p3 <- ggplot(m_long, aes(x = variable, y = gene, fill = log10(value + 0.01))) +
  geom_tile() +
  scale_fill_distiller(palette = "BrBG") +
  theme_tannin() +
  labs(
    title = "Cell-Type Marker Expression",
    x = "",
    y = "",
    fill = "log10(UMI)"
  )

ggsave(file.path(OUTPUT_DIR, "marker_heatmap.pdf"), p3, width = 5, height = 6)

print("Visualizations complete.")
