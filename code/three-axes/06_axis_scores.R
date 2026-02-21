
# 06_axis_scores.R
# Purpose: Summarize Axis Activity using GSEA NES (Statistics) and Z-Scores (Visualization).

source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(fgsea)

out_dir <- "results/three-axes/axis-scores"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. GSEA NES Summary Table
# Load DE Results
de_res <- read_csv("data/RPE_DE_results.csv", show_col_types = FALSE)

# We need GSEA for 3 contrasts: H2O2 vs CTRL, PRG4 vs CTRL, H2O2PRG4 vs H2O2
# Helper to run GSEA
run_gsea_contrast <- function(lfc_col, pval_col, label) {
    # Rank
    min_p <- min(de_res[[pval_col]][de_res[[pval_col]] > 0], na.rm=TRUE)
    stats <- sign(de_res[[lfc_col]]) * -log10(ifelse(de_res[[pval_col]]==0, min_p, de_res[[pval_col]]))
    names(stats) <- de_res$hgnc_symbol
    stats <- stats[!is.na(names(stats))]
    stats <- sort(stats, decreasing = TRUE)
    
    # Pathways
    res <- fgseaMultilevel(pathways = pathways, stats = stats, minSize=15, maxSize=500)
    res$Contrast <- label
    res %>% select(pathway, Contrast, NES, pval, padj)
}

# Load Pathways
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

# Run for all contrasts
gsea_h2o2 <- run_gsea_contrast("H2O2_vs_CTRL.log2FoldChange", "H2O2_vs_CTRL.pvalue", "H2O2_vs_CTRL")
gsea_prg4 <- run_gsea_contrast("PRG4_vs_CTRL.log2FoldChange", "PRG4_vs_CTRL.pvalue", "PRG4_vs_CTRL")
gsea_rescue <- run_gsea_contrast("H2O2PRG4_vs_H2O2.log2FoldChange", "H2O2PRG4_vs_H2O2.pvalue", "H2O2PRG4_vs_H2O2")

all_nes <- rbind(gsea_h2o2, gsea_prg4, gsea_rescue)
write_csv(all_nes, file.path(out_dir, "axis_nes_summary.csv"))
print("GSEA NES Summary:")
print(all_nes)

# 2. Per-Sample Visualization (Boxplots/Radar) using Z-scores
# (We keep this for visual intuition, but rely on NES for stats)
tpms <- read_csv("data/RPE_TPMS.csv", show_col_types = FALSE)
gene_col <- "hgnc_symbol"
tpm_subset <- tpms %>% select(-any_of("ensembl_gene_id"))

long_tpms <- tpm_subset %>%
  pivot_longer(cols = -!!sym(gene_col), names_to = "Sample", values_to = "TPM")

long_tpms$Condition <- case_when(
    long_tpms$Sample %in% c("TS1", "TS2", "TS3") ~ "CTRL",
    long_tpms$Sample %in% c("TS5", "TS6", "TS7") ~ "H2O2",
    long_tpms$Sample %in% c("TS9", "TS11", "TS12") ~ "PRG4",
    long_tpms$Sample %in% c("TS13", "TS15", "TS16") ~ "H2O2+PRG4",
    TRUE ~ "Other"
)
long_tpms <- long_tpms %>% filter(Condition != "Other")

# Z-score per gene
long_tpms <- long_tpms %>%
  group_by(!!sym(gene_col)) %>%
  mutate(Z = (TPM - mean(TPM)) / sd(TPM)) %>%
  ungroup()

# Calculate Score = Mean(Pro-Disease Z)
scores <- data.frame()
samples <- unique(long_tpms$Sample)
for (s in samples) {
  s_dat <- long_tpms %>% filter(Sample == s)
  cond <- unique(s_dat$Condition)
  for (ax in names(pathways)) {
    pro_g <- pathways[[ax]]
    z_vals <- s_dat %>% filter(!!sym(gene_col) %in% pro_g) %>% pull(Z)
    s_score <- mean(z_vals, na.rm=TRUE)
    scores <- rbind(scores, data.frame(Sample=s, Condition=cond, Axis=ax, Score=s_score))
  }
}
scores$Condition <- factor(scores$Condition, levels = c("CTRL", "H2O2", "PRG4", "H2O2+PRG4"))

# Boxplot
p_box <- ggplot(scores, aes(x = Condition, y = Score, fill = Condition)) +
  geom_boxplot() +
  facet_wrap(~Axis, scales = "free_y") +
  theme_pub() +
  labs(title = "Axis Activity (Mean Z-score)", subtitle = "For visual comparison only. See NES for stats.") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "Fig_Axis_scores_boxplot.pdf"), p_box, width = 8, height = 4)

# Radar
radar_dat <- scores %>% group_by(Condition, Axis) %>% summarize(Mean=mean(Score))
min_val <- min(radar_dat$Mean)
radar_dat$Shifted <- radar_dat$Mean - min_val + 0.1 # Make positive

p_radar <- ggplot(radar_dat, aes(x = Axis, y = Shifted, color = Condition, group = Condition)) +
  geom_polygon(fill = NA, size = 1) +
  geom_point() + 
  coord_polar() +
  theme_pub() +
  labs(title = "Axis Radar (Relative Activity)")
ggsave(file.path(out_dir, "Fig_Axis_scores_radar.pdf"), p_radar, width = 5, height = 5)

write_csv(scores, file.path(out_dir, "sample_z_scores_viz_only.csv"))
