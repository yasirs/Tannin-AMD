# Hierarchical Clustering of AMD Samples Using High-Variance Genes
# Creates heatmaps to identify patient subtypes

#%%
# Load libraries
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#%%
# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- DATA_DIR

# Coding preferences
FONT_FAMILY <- "Palatino Linotype"

#%%
# Load data
message("Loading data...")

# Expression data
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

# Variance results (annotated)
macular_var <- read.csv(file.path(DATA_DIR, "macular_variance_analysis_annotated.csv"))
extramacular_var <- read.csv(file.path(DATA_DIR, "extramacular_variance_analysis_annotated.csv"))

# Bimodality results
macular_bimod <- read.csv(file.path(DATA_DIR, "macular_variance_bimodality.csv"))
extramacular_bimod <- read.csv(file.path(DATA_DIR, "extramacular_variance_bimodality.csv"))

#%%
# Function to create heatmap with hierarchical clustering
create_variance_heatmap <- function(expr_data, metadata, var_results, gene_subset, 
                                   title, filename, show_gene_names = TRUE) {
  
  message(sprintf("\nCreating heatmap: %s", title))
  
  # Filter to AMD samples only
  amd_idx <- which(metadata$disease_status == "AMD")
  
  if (length(amd_idx) == 0) {
    message("  No AMD samples found!")
    return(NULL)
  }
  
  # Get gene IDs for subset
  gene_ids <- as.character(gene_subset$gene)
  
  # Filter expression data
  expr_subset <- expr_data[amd_idx, gene_ids, drop = FALSE]
  
  # Convert to matrix and transpose (genes as rows, samples as columns)
  expr_mat <- t(as.matrix(expr_subset))
  
  # Get gene symbols for row labels
  gene_symbols <- gene_subset$gene_symbol
  gene_symbols[is.na(gene_symbols) | gene_symbols == ""] <- 
    paste0("Gene_", gene_ids[is.na(gene_symbols) | gene_symbols == ""])
  
  rownames(expr_mat) <- gene_symbols
  
  # Z-score normalize genes (row-wise)
  expr_mat_scaled <- t(scale(t(expr_mat)))
  
  # Create annotation for samples
  annotation_col <- data.frame(
    Stage = metadata$dry_stage[amd_idx],
    Age = metadata$age[amd_idx],
    row.names = rownames(expr_subset)
  )
  
  # Annotation colors
  stage_colors <- c(
    "MD1" = "#fee391",
    "MD2" = "#fec44f", 
    "dry_AMD" = "#fe9929",
    "GA" = "#cc4c02"
  )
  
  annotation_colors <- list(
    Stage = stage_colors,
    Age = colorRampPalette(c("lightblue", "darkblue"))(100)
  )
  
  # Heatmap color palette (BrBG diverging)
  heatmap_colors <- colorRampPalette(rev(brewer.pal(11, "BrBG")))(100)
  
  # Create heatmap
  message(sprintf("  Samples: %d, Genes: %d", ncol(expr_mat_scaled), nrow(expr_mat_scaled)))
  
  # Determine if gene names should be shown
  show_rownames <- show_gene_names && nrow(expr_mat_scaled) <= 100
  
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
    fontsize_col = 8,
    main = title,
    filename = file.path(OUTPUT_DIR, filename),
    width = 10,
    height = ifelse(nrow(expr_mat_scaled) > 50, 14, 10)
  )
  
  # Also save as TIFF with LZW compression
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
    fontsize_col = 8,
    main = title,
    filename = file.path(OUTPUT_DIR, tiff_file),
    width = 10,
    height = ifelse(nrow(expr_mat_scaled) > 50, 14, 10),
    compression = "lzw",
    res = 300
  )
  
  message(sprintf("  Saved: %s", filename))
  message(sprintf("  Saved: %s", tiff_file))
  
  return(expr_mat_scaled)
}

#%%
# ========================================
# MACULAR TISSUE
# ========================================

message("\n=== MACULAR RPE-CHOROID ===")

# 1. All high-variance genes (FDR < 0.05, increased variance)
macular_highvar <- macular_var %>%
  filter(f_fdr < 0.05 & var_fc > 0) %>%
  arrange(f_fdr)

message(sprintf("\nHigh-variance genes (FDR<0.05, var_fc>0): %d", nrow(macular_highvar)))

if (nrow(macular_highvar) > 0) {
  create_variance_heatmap(
    expr_data = macular_data$expr,
    metadata = macular_data$metadata,
    var_results = macular_var,
    gene_subset = macular_highvar,
    title = "Macular AMD: High-Variance Genes (FDR<0.05)",
    filename = "macular_amd_clustering_highvar.pdf",
    show_gene_names = TRUE
  )
}

# 2. Bimodal genes only
macular_bimodal_genes <- macular_bimod %>%
  filter(is_bimodal == TRUE) %>%
  arrange(desc(bic_improvement))

message(sprintf("\nBimodal genes: %d", nrow(macular_bimodal_genes)))

if (nrow(macular_bimodal_genes) > 0) {
  create_variance_heatmap(
    expr_data = macular_data$expr,
    metadata = macular_data$metadata,
    var_results = macular_var,
    gene_subset = macular_bimodal_genes,
    title = "Macular AMD: Bimodal Genes (Patient Subtypes)",
    filename = "macular_amd_clustering_bimodal.pdf",
    show_gene_names = TRUE
  )
}

#%%
# ========================================
# EXTRAMACULAR TISSUE
# ========================================

message("\n=== EXTRAMACULAR RPE-CHOROID ===")

# 1. All high-variance genes (FDR < 0.05, increased variance)
extramacular_highvar <- extramacular_var %>%
  filter(f_fdr < 0.05 & var_fc > 0) %>%
  arrange(f_fdr)

message(sprintf("\nHigh-variance genes (FDR<0.05, var_fc>0): %d", nrow(extramacular_highvar)))

if (nrow(extramacular_highvar) > 0) {
  create_variance_heatmap(
    expr_data = extramacular_data$expr,
    metadata = extramacular_data$metadata,
    var_results = extramacular_var,
    gene_subset = extramacular_highvar,
    title = "Extramacular AMD: High-Variance Genes (FDR<0.05)",
    filename = "extramacular_amd_clustering_highvar.pdf",
    show_gene_names = FALSE  # Too many genes
  )
}

# 2. Bimodal genes only
extramacular_bimodal_genes <- extramacular_bimod %>%
  filter(is_bimodal == TRUE) %>%
  arrange(desc(bic_improvement))

message(sprintf("\nBimodal genes: %d", nrow(extramacular_bimodal_genes)))

if (nrow(extramacular_bimodal_genes) > 0) {
  create_variance_heatmap(
    expr_data = extramacular_data$expr,
    metadata = extramacular_data$metadata,
    var_results = extramacular_var,
    gene_subset = extramacular_bimodal_genes,
    title = "Extramacular AMD: Bimodal Genes (Patient Subtypes)",
    filename = "extramacular_amd_clustering_bimodal.pdf",
    show_gene_names = TRUE
  )
}

#%%
message("\n=== CLUSTERING COMPLETE ===")
message("Generated 4 heatmaps (PDF + TIFF):")
message("  - macular_amd_clustering_highvar.pdf/.tiff")
message("  - macular_amd_clustering_bimodal.pdf/.tiff")
message("  - extramacular_amd_clustering_highvar.pdf/.tiff")
message("  - extramacular_amd_clustering_bimodal.pdf/.tiff")
