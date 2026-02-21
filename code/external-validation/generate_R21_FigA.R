
library(ggplot2)
library(dplyr)
library(data.table)

# Set Paths
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE29801/age_nrf2_analysis")
output_file_tiff <- file.path(results_dir, "R21_FigA.tiff")
output_file_pdf <- file.path(results_dir, "R21_FigA.pdf")

# Load Data
cat("Loading data...\n")
aging_df <- fread(file.path(results_dir, "backward_aging_mimetics.csv"))
prg4_df <- fread(file.path(results_dir, "backward_prg4_mimetics.csv"))
are_df <- fread(file.path(results_dir, "ARE_gene_set.csv"))
are_genes <- are_df$gene

# Merge Data
# We want target_gene, mimetic_score (aging), prg4_mimetic_score
merged_df <- merge(aging_df[, .(perturbation_id, target_gene, mimetic_score)], 
                   prg4_df[, .(perturbation_id, prg4_mimetic_score)], 
                   by = "perturbation_id")

# Calculate Attributes
merged_df[, `:=`(
  aging_antagonism = -mimetic_score,
  is_ARE = target_gene %in% are_genes
)]



# Define Top/Bottom Thresholds (Fix: strict quadrant filtering)
# User asked "why is one point in a different quadrant?"
# Cause: We selected by GENE, so if a gene had 2 perturbations (one good, one bad), both were colored.
# Fix: 1. Filter for Q1 (Positive/Positive). 2. Select top 10. 3. Color only those specific points.

merged_df[, combined_score := aging_antagonism + prg4_mimetic_score]

# Filter for strictly positive quadrant candidates
q1_candidates <- merged_df[aging_antagonism > 0 & prg4_mimetic_score > 0]

# Take top 10 unique genes from this high-quality set
top_10_genes <- head(unique(q1_candidates[order(-combined_score)]$target_gene), 10)

# Identify the specific perturbations for these genes that are ALSO in Q1
# This implies we might color fewer than all perturbations for a gene, which is correct for visual consistency
convergent_perturbations <- q1_candidates[target_gene %in% top_10_genes]$perturbation_id

cat(sprintf("Highlighting %d specific perturbations for top 10 convergent genes (Q1 only)\n", length(convergent_perturbations)))

# Define Groups
merged_df[, plot_group := "Other"]
merged_df[perturbation_id %in% convergent_perturbations, plot_group := "Convergent"]

# Set Factor levels (draw Other first)
merged_df[, plot_group := factor(plot_group, levels = c("Other", "Convergent"))]
merged_df <- merged_df[order(plot_group)]

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

# Labels - Top 5 genes (attach label to the best scoring perturbation for that gene)
top_5_genes <- head(top_10_genes, 5)
# Select best perturbation per gene for labeling to avoid double labels
label_data <- merged_df[perturbation_id %in% convergent_perturbations & target_gene %in% top_5_genes]
label_data <- label_data[order(-combined_score)]
label_data <- label_data[!duplicated(target_gene)]

# Try to load ggrepel for better labeling
use_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if(use_ggrepel) library(ggrepel)

cat("Generating plot...\n")
p <- ggplot(merged_df, aes(x = aging_antagonism, y = prg4_mimetic_score, color = plot_group, size = plot_group, alpha = plot_group)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(stroke = 0) +
  scale_color_manual(values = color_map) +
  scale_alpha_manual(values = alpha_map) +
  scale_size_manual(values = size_map) +
  # Highlight circle (matched to new smaller size)
  geom_point(data = merged_df[plot_group == "Convergent"], shape = 1, color = "black", size = 0.6, stroke = 0.2, alpha = 1) + 
  labs(
    x = "Aging Antagonism Score",
    y = "PRG4 Mimetic Score",
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
if(use_ggrepel) {
  p <- p + geom_text_repel(data = label_data, aes(label = target_gene), 
                           size = 2.5, color = "black", fontface = "bold", segment.size = 0.2,
                           min.segment.length = 0, show.legend = FALSE)
} else {
  p <- p + geom_text(data = label_data, aes(label = target_gene), 
                     size = 2.5, color = "black", fontface = "bold", vjust = -1, show.legend = FALSE)
}

# Save
cat("Saving to", output_file_tiff, "and PDF...\n")
ggsave(output_file_pdf, p, width = 2.5, height = 2.5, units = "in")

# For TIFF
tiff(output_file_tiff, width = 2.5, height = 2.5, units = "in", res = 300, compression = "lzw")
print(p)
dev.off()

cat("Done.\n")
