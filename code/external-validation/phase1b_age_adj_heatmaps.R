# Phase 1b: Create Age-Adjusted Variance Heatmaps
# Hierarchical clustering using age-adjusted high-variance and bimodal genes

#%%
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#%%
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- file.path(BASE_DIR, "results/external-validation/gse29801_prg4_comparison")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")

dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

FONT_FAMILY <- "Palatino Linotype"

#%%
# Load data
message("Loading data...")
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

# Load age-adjusted gene sets
macular_high <- read.csv(file.path(OUTPUT_DIR, "gene_sets/macular_high_variance_age_adj.csv"))
macular_bimod <- read.csv(file.path(OUTPUT_DIR, "gene_sets/macular_bimodality_age_adjusted.csv"))
extramacular_high <- read.csv(file.path(OUTPUT_DIR, "gene_sets/extramacular_high_variance_age_adj.csv"))
extramacular_bimod <- read.csv(file.path(OUTPUT_DIR, "gene_sets/extramacular_bimodality_age_adjusted.csv"))

#%%
# Function to create heatmap
create_heatmap <- function(expr_data, metadata, gene_subset, title, filename) {
  
  message(sprintf("\nCreating: %s", title))
  
  # Filter to AMD samples
  amd_idx <- which(metadata$disease_status == "AMD")
  
  if (length(amd_idx) == 0 || nrow(gene_subset) == 0) {
    message("  Skipping - no data")
    return(NULL)
  }
  
  # Get gene IDs
  gene_ids <- as.character(gene_subset$gene)
  expr_subset <- expr_data[amd_idx, gene_ids, drop = FALSE]
  
  # Transpose and z-score
  expr_mat <- t(as.matrix(expr_subset))
  expr_mat_scaled <- t(scale(t(expr_mat)))
  
  # Gene labels
  gene_symbols <- gene_subset$gene_symbol
  gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- 
    paste0("Gene_", gene_ids[is.na(gene_symbols) | gene_symbols == ""])
  rownames(expr_mat_scaled) <- gene_symbols
  
  # Sample annotation
  annotation_col <- data.frame(
    Stage = metadata$dry_stage[amd_idx],
    Age = metadata$age[amd_idx],
    row.names = rownames(expr_subset)
  )
  
  # Colors
  stage_colors <- c("MD1" = "#fee391", "MD2" = "#fec44f", 
                   "dry_AMD" = "#fe9929", "GA" = "#cc4c02")
  annotation_colors <- list(
    Stage = stage_colors,
    Age = colorRampPalette(c("lightblue", "darkblue"))(100)
  )
  
  heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "BrBG")))(100)
  
  # Determine if show gene names
  show_rownames <- nrow(expr_mat_scaled) <= 100
  
  # Create heatmap
  message(sprintf("  Samples: %d, Genes: %d", ncol(expr_mat_scaled), nrow(expr_mat_scaled)))
  
  pheatmap(
    expr_mat_scaled,
    color = heatmap_colors,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    show_rownames = show_rownames,
    show_colnames = FALSE,
    fontsize_row = 8,
    main = title,
    filename = file.path(FIGURE_DIR, filename),
    width = 10,
    height = ifelse(nrow(expr_mat_scaled) > 50, 14, 10)
  )
  
  # TIFF version
  tiff_file <- sub("\\.pdf$", ".tiff", filename)
  pheatmap(
    expr_mat_scaled,
    color = heatmap_colors,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "complete",
    show_rownames = show_rownames,
    show_colnames = FALSE,
    fontsize_row = 8,
    main = title,
    filename = file.path(FIGURE_DIR, tiff_file),
    width = 10,
    height = ifelse(nrow(expr_mat_scaled) > 50, 14, 10),
    compression = "lzw",
    res = 300
  )
  
  message(sprintf("  Saved: %s", filename))
}

#%%
# Create 4 heatmaps

# Macular - High variance
create_heatmap(
  macular_data$expr,
  macular_data$metadata,
  macular_high,
  "Macular AMD: High-Variance Genes (Age-Adjusted)",
  "macular_amd_clustering_highvar_age_adj.pdf"
)

# Macular - Bimodal
macular_bimod_filtered <- macular_bimod %>% filter(is_bimodal == TRUE)
create_heatmap(
  macular_data$expr,
  macular_data$metadata,
  macular_bimod_filtered,
  "Macular AMD: Bimodal Genes (Age-Adjusted)",
  "macular_amd_clustering_bimodal_age_adj.pdf"
)

# Extramacular - High variance
create_heatmap(
  extramacular_data$expr,
  extramacular_data$metadata,
  extramacular_high,
  "Extramacular AMD: High-Variance Genes (Age-Adjusted)",
  "extramacular_amd_clustering_highvar_age_adj.pdf"
)

# Extramacular - Bimodal
extramacular_bimod_filtered <- extramacular_bimod %>% filter(is_bimodal == TRUE)
create_heatmap(
  extramacular_data$expr,
  extramacular_data$metadata,
  extramacular_bimod_filtered,
  "Extramacular AMD: Bimodal Genes (Age-Adjusted)",
  "extramacular_amd_clustering_bimodal_age_adj.pdf"
)

message("\n=== Heatmaps Complete ===")
message(sprintf("Saved to: %s", FIGURE_DIR))
