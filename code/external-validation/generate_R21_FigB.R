
library(data.table)
library(ggplot2)
library(ggrepel)

# Set Paths
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE29801/age_nrf2_analysis")
axes_dir <- file.path(base_dir, "results/perturbseq-analysis/Axes")
output_file_tiff <- file.path(results_dir, "R21_FigB.tiff")
output_file_pdf <- file.path(results_dir, "R21_FigB.pdf")

# Load Data
cat("Loading scores...\n")
aging_df <- fread(file.path(results_dir, "backward_aging_mimetics.csv"))
prg4_df <- fread(file.path(results_dir, "backward_prg4_mimetics.csv"))
redox_df <- fread(file.path(axes_dir, "Redox/RPE1/backward/backward_Redox_mimetics.csv"), 
                   select = c("perturbation_id", "mimetic_score"))
setnames(redox_df, "mimetic_score", "redox_mimetic_score")

# Merge All
cat("Merging...\n")
merged <- merge(aging_df[, .(perturbation_id, target_gene, mimetic_score)], 
                prg4_df[, .(perturbation_id, prg4_mimetic_score)], 
                by = "perturbation_id")
merged <- merge(merged, redox_df, by = "perturbation_id")

# Calculate Antagonism
merged[, `:=`(
  aging_antagonism = -mimetic_score,
  redox_antagonism = -redox_mimetic_score
)]

# Identify Convergent Genes (Using SAME logic as Fig A)
# Logic: Top 10 unique genes by (Aging Antag + PRG4 Mimetic), strictly in Q1 (Positive/Positive)
merged[, combined_score_A := aging_antagonism + prg4_mimetic_score]
q1_candidates_A <- merged[aging_antagonism > 0 & prg4_mimetic_score > 0]
top_10_genes <- head(unique(q1_candidates_A[order(-combined_score_A)]$target_gene), 10)
convergent_perturbations <- q1_candidates_A[target_gene %in% top_10_genes]$perturbation_id

cat(sprintf("Highlighting %d convergent perturbations (identified from Fig A logic)\n", length(convergent_perturbations)))

# Define Groups
merged[, plot_group := "Other"]
merged[perturbation_id %in% convergent_perturbations, plot_group := "Convergent"]

# Set Factor levels
merged[, plot_group := factor(plot_group, levels = c("Other", "Convergent"))]
merged <- merged[order(plot_group)]

# Colors
color_map <- c(
  "Other" = "#BDBDBD",      # Grey
  "Convergent" = "#8E44AD"  # Purple
)

# Alpha map
alpha_map <- c(
  "Other" = 0.4,
  "Convergent" = 1.0
)

# Size map
size_map <- c(
  "Other" = 0.1,
  "Convergent" = 0.6
)

# Labels - Top 5 genes (best perturbation per gene)
top_5_genes <- head(top_10_genes, 5)
label_data <- merged[perturbation_id %in% convergent_perturbations & target_gene %in% top_5_genes]
# Pick best by combined score A (Aging+PRG4) to use the SAME labels as Fig A
label_data <- label_data[order(-combined_score_A)]
label_data <- label_data[!duplicated(target_gene)]

cat("Generating plot...\n")
p <- ggplot(merged, aes(x = aging_antagonism, y = redox_antagonism, color = plot_group, size = plot_group, alpha = plot_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(stroke = 0) +
  scale_color_manual(values = color_map) +
  scale_alpha_manual(values = alpha_map) +
  scale_size_manual(values = size_map) +
  # Highlight circle
  geom_point(data = merged[plot_group == "Convergent"], shape = 1, color = "black", size = 0.6, stroke = 0.2, alpha = 1) + 
  labs(
    x = "Aging Antagonism Score",
    y = "Redox Antagonism Score",
    title = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 7),
    axis.title = element_text(size = 8),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

# Add labels
p <- p + geom_text_repel(data = label_data, aes(label = target_gene), 
                         size = 2.5, color = "black", fontface = "bold", segment.size = 0.2,
                         min.segment.length = 0, show.legend = FALSE)

# Save
cat("Saving to", output_file_tiff, "and PDF...\n")
ggsave(output_file_pdf, p, width = 2.5, height = 2.5, units = "in")

# For TIFF
tiff(output_file_tiff, width = 2.5, height = 2.5, units = "in", res = 300, compression = "lzw")
print(p)
dev.off()

cat("Done.\n")
