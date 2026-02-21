
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(tidyr)

# Set Paths
base_dir <- "/home/ysuhail/work/Tannin-AMD"
tpm_file <- file.path(base_dir, "data/RPE_TPMS.csv")
results_dir <- file.path(base_dir, "results/cohort-GSE29801/age_nrf2_analysis")
output_file_tiff <- file.path(results_dir, "R21_FigD.tiff")
output_file_pdf <- file.path(results_dir, "R21_FigD.pdf")

# Load Data
cat("Loading TPMs...\n")
tpm_df <- fread(tpm_file)

# Filter for Genes of Interest
target_genes <- c("HMOX1", "NQO1", "SQSTM1", "GPX3", "GCLC", "TXNIP", "GDF15")
subset_df <- tpm_df[hgnc_symbol %in% target_genes]

# Define Samples
samples_ctrl <- c("TS1", "TS2", "TS3")
samples_stress <- c("TS5", "TS6", "TS7")
samples_rescue <- c("TS13", "TS15", "TS16")
all_samples <- c(samples_ctrl, samples_stress, samples_rescue)

# Extract Matrix
mat_data <- as.matrix(subset_df[, ..all_samples])
rownames(mat_data) <- subset_df$hgnc_symbol

# Compute Row-wise Z-scores
cat("Computing Z-scores...\n")
z_mat <- t(apply(mat_data, 1, scale))
colnames(z_mat) <- colnames(mat_data)

# Convert to Long Format for ggplot
z_df <- as.data.frame(z_mat)
z_df$Gene <- rownames(z_df)
long_df <- pivot_longer(z_df, cols = -Gene, names_to = "Sample", values_to = "Z")

# Add Condition
long_df$Condition <- NA
long_df$Condition[long_df$Sample %in% samples_ctrl] <- "Cntrl"
long_df$Condition[long_df$Sample %in% samples_stress] <- "H2O2"
long_df$Condition[long_df$Sample %in% samples_rescue] <- "H2O2+PRG4"

# Set Factor Levels
long_df$Condition <- factor(long_df$Condition, levels = c("Cntrl", "H2O2", "H2O2+PRG4"))

# Order Genes: Manual ordering or Clustering order?
# User said "no dendrogram", so manual or simple sort is fine. 
# Let's cluster once just to get a nice visual order, then fix the factor levels.
hc <- hclust(dist(z_mat))
gene_order <- rownames(z_mat)[hc$order]
long_df$Gene <- factor(long_df$Gene, levels = gene_order)

# Colors: BrBG reversed (Brown=High)
# brewer.pal(11, "BrBG") -> Brown...Green
# rev -> Green...Brown
colors <- rev(brewer.pal(11, "BrBG"))

cat("Generating ggplot heatmap...\n")
p <- ggplot(long_df, aes(x = Sample, y = Gene, fill = Z)) +
  geom_tile(color = "white", size = 0.2) +
  facet_grid(~Condition, scales = "free_x", space = "free_x") +
  scale_fill_gradientn(colours = colors, name = "Z-score") +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(), # Remove sample names if redundant with facet
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    strip.text = element_text(size = 9, face = "bold", color = "black"),
    strip.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.4, "cm")
  )

# Save
cat("Saving to", output_file_tiff, "and PDF...\n")
ggsave(output_file_pdf, p, width = 4, height = 2.5, units = "in")

# For TIFF
tiff(output_file_tiff, width = 4, height = 2.5, units = "in", res = 300, compression = "lzw")
print(p)
dev.off()

cat("Done.\n")
