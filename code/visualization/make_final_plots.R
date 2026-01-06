library(ggplot2)
library(extrafont)
library(reshape2)
library(gridExtra)

# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
EXT_DIR <- file.path(BASE_DIR, "results/external-validation")
RES_DIR <- file.path(BASE_DIR, "results/visualization")
if(!dir.exists(RES_DIR)) dir.create(RES_DIR)

# Custom Theme (Tufte-inspired)
theme_tufte_custom <- function() {
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    axis.text = element_text(color = "black", size = 10),
    axis.line = element_blank(), # Remove axis lines (data-ink)
    axis.ticks = element_line(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2), # Subtle grid
    plot.title = element_text(hjust = 0, face = "bold", size = 14),
    legend.position = "top",
    legend.key = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0)
  )
}

# 1. Serum Starvation vs PRG4 Rescue
df_serum <- read.csv(file.path(EXT_DIR, "gse129964_comparison.csv"))

# Classify
df_serum$Category <- "Neutral"
df_serum$Category[df_serum$logfc_serum_stress > 0.5 & df_serum$H2O2PRG4_vs_H2O2.log2FoldChange < -0.5] <- "Restored (Pathogenic)"
df_serum$Category[df_serum$logfc_serum_stress < -0.5 & df_serum$H2O2PRG4_vs_H2O2.log2FoldChange > 0.5] <- "Restored (Protective)"

# Filter for plotting (remove noise)
df_plot <- df_serum[abs(df_serum$logfc_serum_stress) > 0.2 | abs(df_serum$H2O2PRG4_vs_H2O2.log2FoldChange) > 0.2, ]

p1 <- ggplot(df_plot, aes(x = logfc_serum_stress, y = H2O2PRG4_vs_H2O2.log2FoldChange)) +
  geom_point(aes(color = Category, alpha = Category), size = 1) +
  scale_color_manual(values = c("Neutral" = "grey80", "Restored (Pathogenic)" = "#D95F02", "Restored (Protective)" = "#1B9E77")) +
  scale_alpha_manual(values = c("Neutral" = 0.3, "Restored (Pathogenic)" = 0.8, "Restored (Protective)" = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_tufte_custom() +
  labs(
    title = "PRG4 Reverses Serum Starvation Transcriptome",
    subtitle = "Genes induced by starvation (x>0) are suppressed by PRG4 (y<0)",
    x = "Serum Starvation Log2 Fold Change",
    y = "PRG4 Rescue Log2 Fold Change"
  ) +
  geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.5) +
  annotate("text", x = 2, y = -2, label = paste0("r = ", round(cor(df_serum$logfc_serum_stress, df_serum$H2O2PRG4_vs_H2O2.log2FoldChange), 2)), family = "serif")

ggsave(file.path(RES_DIR, "serum_rescue_scatter.pdf"), p1, width = 6, height = 6)

# 2. AMD Patient Distributions (Small Multiples)
df_expr <- read.csv(file.path(EXT_DIR, "GSE135092_risk_expression.csv"))
# Melt
df_long <- melt(df_expr, id.vars = c("X", "Group"), variable.name = "Gene", value.name = "LogCPM")
df_long <- df_long[df_long$Group %in% c("AMD", "Control"), ]

# Density Ridgeline style (using facet_wrap)
p2 <- ggplot(df_long, aes(x = LogCPM, fill = Group)) +
  geom_density(alpha = 0.6, size = 0.2) +
  facet_wrap(~Gene, scales = "free", ncol = 4) +
  scale_fill_manual(values = c("Control" = "#7570B3", "AMD" = "#D95F02")) +
  theme_tufte_custom() +
  labs(
    title = "Expression of Risk Genes in Human RPE/Choroid",
    subtitle = "Distributions across 537 samples (GSE135092)",
    x = "Expression (Log2 CPM)",
    y = "Density"
  )

ggsave(file.path(RES_DIR, "amd_patient_distributions.pdf"), p2, width = 10, height = 5)

print("Visualizations complete.")
