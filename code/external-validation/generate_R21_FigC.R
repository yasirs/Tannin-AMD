
library(data.table)
library(pheatmap)
library(ggplot2)

# Set Paths
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE29801/age_nrf2_analysis")
axes_dir <- file.path(base_dir, "results/perturbseq-analysis/Axes")
output_file_tiff <- file.path(results_dir, "R21_FigC.tiff")
output_file_pdf <- file.path(results_dir, "R21_FigC.pdf")

# Load Data
cat("Loading axis mimetics...\n")
aging_df <- fread(file.path(results_dir, "backward_aging_mimetics.csv"), select = c("perturbation_id", "target_gene", "mimetic_score"))
setnames(aging_df, "mimetic_score", "Aging")

prg4_df <- fread(file.path(results_dir, "backward_prg4_mimetics.csv"), select = c("perturbation_id", "prg4_mimetic_score"))
setnames(prg4_df, "prg4_mimetic_score", "PRG4")

# Helper to load axis
load_axis <- function(name, path) {
  df <- fread(path, select = c("perturbation_id", "mimetic_score"))
  setnames(df, "mimetic_score", name)
  return(df)
}

redox_df <- load_axis("Redox", file.path(axes_dir, "Redox/RPE1/backward/backward_Redox_mimetics.csv"))
sen_df <- load_axis("Senescence", file.path(axes_dir, "Senescence/RPE1/backward/backward_Senescence_mimetics.csv"))
inf_df <- load_axis("Inflammation", file.path(axes_dir, "Inflammation/RPE1/backward/backward_Inflammation_mimetics.csv"))

# Merge All
cat("Merging scores...\n")
merged <- merge(aging_df, prg4_df, by = "perturbation_id", all = FALSE)
merged <- merge(merged, redox_df, by = "perturbation_id", all = FALSE)
merged <- merge(merged, sen_df, by = "perturbation_id", all = FALSE)
merged <- merge(merged, inf_df, by = "perturbation_id", all = FALSE)

# Select Representative Perturbations
# 1. Convergent: XPO5, FGFR1OP, RBBP6 (from Fig A)
# Use grep to find the perturbation IDs for these genes. We want the ones that performed well.
# Ideally, we sort by optimal combined score first to pick the "good" perturbation if multiple exist
merged[, combined_score := -Aging + PRG4] # Antagonism + Mimetic
merged <- merged[order(-combined_score)]

get_best_id <- function(gene) {
  return(merged[target_gene == gene, perturbation_id][1])
}

conv_genes <- c("XPO5", "FGFR1OP", "RBBP6", "DICER1")
conv_ids <- sapply(conv_genes, get_best_id)
conv_ids <- conv_ids[!is.na(conv_ids)]

# 2. Exclusive Aging Antagonists
# Low Aging score, but NOT high PRG4 score (to be distinct)
# Sort by Aging (ascending), filtering for low PRG4
antag_candidates <- merged[PRG4 < 0.05][order(Aging)]
antag_ids <- head(antag_candidates$perturbation_id, 3)

# 3. Exclusive PRG4 Mimetics
# High PRG4, but 'normal' Aging (or at least not super low)
# Sort by PRG4 (descending), filtering for Aging > -0.1
prg4_candidates <- merged[Aging > -0.1][order(-PRG4)]
prg4_ids <- head(prg4_candidates$perturbation_id, 3)

# 4. Neutral
# Scores near 0 on all axes
neutral_candidates <- merged[abs(Aging) < 0.05 & abs(PRG4) < 0.05 & abs(Redox) < 0.05]
neutral_ids <- head(neutral_candidates$perturbation_id, 3)


# Combine selection
selected_ids <- c(conv_ids, antag_ids, prg4_ids, neutral_ids)
selected_df <- merged[perturbation_id %in% selected_ids]

# Define Annotations (Categories)
selected_df[, Category := "Neutral"]
selected_df[perturbation_id %in% conv_ids, Category := "Convergent"]
selected_df[perturbation_id %in% antag_ids, Category := "Aging Antagonist"]
selected_df[perturbation_id %in% prg4_ids, Category := "PRG4 Mimetic"]

# Set factor levels for ordering
selected_df[, Category := factor(Category, levels = c("Convergent", "Aging Antagonist", "PRG4 Mimetic", "Neutral"))]
# Sort by Category so rows are grouped in heatmap
selected_df <- selected_df[order(Category)]

# Prepare for Heatmap
# Rows: Target Gene
# Columns: Aging, Redox, Senescence, Inflammation (and PRG4 for context)
mat_data <- as.matrix(selected_df[, .(Aging, Redox, Senescence, Inflammation, PRG4)])
rownames(mat_data) <- selected_df$target_gene

# Define Groups for Annotation
annot_row <- data.frame(row.names = selected_df$target_gene, Category = selected_df$Category)

# Annot Colors
ann_colors <- list(
  Category = c(
    "Convergent" = "#8E44AD",
    "Aging Antagonist" = "#27AE60",
    "PRG4 Mimetic" = "#3498DB",
    "Neutral" = "#BDBDBD"
  )
)

# Colors: BrBG (Brown for high values)
# Standard BrBG is Brown(Low) -> Green(High). We want High=Brown.
# So we want Green(Low) -> Brown(High).
# Note: RColorBrewer::brewer.pal(11, "BrBG") gives [Brown, ..., Green]
# So rev() gives [Green, ..., Brown]
heatmap_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG")))(100)

# Plot
# 3x3 inches is small, adjust logical size
cat("Generating heatmap...\n")

# Save
# TIFF
tiff(output_file_tiff, width = 4, height = 4, units = "in", res = 300, compression = "lzw")
pheatmap(mat_data,
         color = heatmap_colors,
         scale = "none", # Scores are already comparable (correlation-based/z-score-based)
         cluster_rows = FALSE, # User requested no clustering, just grouping
         cluster_cols = FALSE, # Keep axes ordered
         annotation_row = annot_row,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         border_color = NA,
         treeheight_row = 10,
         fontsize = 8,
         main = "Multi-axis Fingerprint"
)
dev.off()

# PDF
pdf(output_file_pdf, width = 4, height = 4)
pheatmap(mat_data,
         color = heatmap_colors,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annot_row,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         border_color = NA,
         treeheight_row = 10,
         fontsize = 8,
         main = "Multi-axis Fingerprint"
)
dev.off()

cat("Done.\n")
