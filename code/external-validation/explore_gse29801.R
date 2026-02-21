# GSE29801 Exploratory Visualization Script
# Following coding_preferences.md: R, ggplot2, theme_classic, Palatino Linotype, BrBG

#%%
# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

#%%
# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "data/external/geo/GSE29801")
OUTPUT_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801")

# Visualization settings (from coding_preferences.md)
# Note: Using 'serif' as fallback since Palatino Linotype may not be installed
FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

#%%
# Parse series matrix to extract metadata
parse_series_matrix <- function(path) {
  con <- gzfile(path, "rt")
  lines <- readLines(con)
  close(con)
  
  # Find sample characteristics lines
  char_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
  geo_line <- lines[grepl("^!Sample_geo_accession", lines)]
  
  # Parse sample IDs
  sample_ids <- strsplit(geo_line, "\t")[[1]][-1]
  sample_ids <- gsub("\"", "", sample_ids)
  
  # Initialize metadata data frame
  metadata <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  
  # Parse each characteristic line
  for (line in char_lines) {
    parts <- strsplit(line, "\t")[[1]][-1]
    parts <- gsub("\"", "", parts)
    
    # Extract key-value pairs
    if (length(parts) > 0 && grepl(":", parts[1])) {
      key <- trimws(sub(":.*", "", parts[1]))
      values <- sapply(parts, function(x) trimws(sub("^[^:]+:\\s*", "", x)))
      metadata[[key]] <- values
    }
  }
  
  return(metadata)
}

#%%
# Load metadata
message("Loading metadata from series matrix...")
metadata <- parse_series_matrix(file.path(DATA_DIR, "GSE29801_series_matrix.txt.gz"))

# Clean up column names
colnames(metadata) <- gsub(" ", "_", colnames(metadata))
colnames(metadata) <- gsub("[()]", "", colnames(metadata))

# Create simplified disease status
metadata$disease_status <- ifelse(metadata$ocular_disease == "normal", "Control", "AMD")

# Create staging categories
metadata$staging <- case_when(
  metadata$amd_classification == "normal" ~ "Control",
  metadata$amd_classification %in% c("MD1", "MD2") ~ "Early",
  metadata$amd_classification == "dry AMD" ~ "Intermediate",
  metadata$amd_classification %in% c("GA", "CNV", "GA/CNV") ~ "Advanced",
  TRUE ~ "Other AMD"
)
metadata$staging <- factor(metadata$staging, 
                           levels = c("Control", "Early", "Intermediate", "Advanced", "Other AMD"))

# Parse age as numeric
metadata$age <- as.numeric(metadata$age_years)

message(sprintf("Loaded metadata for %d samples", nrow(metadata)))
print(table(metadata$disease_status))

#%%
# Load expression matrix
message("Loading expression matrix...")
expr_raw <- read.csv(gzfile(file.path(DATA_DIR, "GSE29801_expression_matrix.csv.gz")), 
                     row.names = 1, check.names = FALSE)

# Separate expression data from Group column
if ("Group" %in% colnames(expr_raw)) {
  expr_data <- expr_raw[, colnames(expr_raw) != "Group"]
} else {
  expr_data <- expr_raw
}

# Ensure metadata order matches expression data
metadata <- metadata[match(rownames(expr_data), metadata$sample_id), ]

message(sprintf("Expression matrix: %d samples × %d probes", nrow(expr_data), ncol(expr_data)))

#%%
# Color definitions
staging_colors <- c(
  "Control" = "navy",
  "Early" = "#fee391",
  "Intermediate" = "#fe9929",
  "Advanced" = "#cc4c02",
  "Other AMD" = "#662506"
)

# BrBG color palette for heatmap (from coding preferences)
heatmap_palette <- colorRampPalette(rev(brewer.pal(11, "BrBG")))(100)

#%%
# Contrastive PCA implementation
# Identifies directions of variation enriched in foreground (disease) vs background (control)
contrastive_pca <- function(foreground, background, n_components = 2, alpha = 1.0) {
  # Center both datasets
  fg_centered <- scale(foreground, center = TRUE, scale = FALSE)
  bg_centered <- scale(background, center = TRUE, scale = FALSE)
  
  # Compute covariance matrices
  cov_fg <- cov(fg_centered)
  cov_bg <- cov(bg_centered)
  
  # Contrastive covariance: C_fg - alpha * C_bg
  cov_contrastive <- cov_fg - alpha * cov_bg
  
  # Eigen decomposition
  eigen_result <- eigen(cov_contrastive)
  
  # Get top components
  loadings <- eigen_result$vectors[, 1:n_components, drop = FALSE]
  eigenvalues <- eigen_result$values[1:n_components]
  
  # Project both datasets
  fg_projected <- fg_centered %*% loadings
  bg_projected <- bg_centered %*% loadings
  
  # Combine projections
  all_projected <- rbind(bg_projected, fg_projected)
  
  # Calculate variance explained (as proportion of contrastive variance)
  total_var <- sum(abs(eigenvalues))
  var_explained <- abs(eigenvalues) / total_var
  
  return(list(
    projection = all_projected,
    loadings = loadings,
    eigenvalues = eigenvalues,
    var_explained = var_explained
  ))
}

#%%
# Function to generate PCA and heatmap for a subset of data
generate_tissue_plots <- function(expr_subset, meta_subset, tissue_name, output_dir) {
  
  # Skip if too few samples
  if (nrow(expr_subset) < 10) {
    message(sprintf("  Skipping %s: only %d samples", tissue_name, nrow(expr_subset)))
    return(NULL)
  }
  
  # Clean tissue name for filenames
  tissue_slug <- gsub(" ", "_", tolower(tissue_name))
  tissue_slug <- gsub("-", "_", tissue_slug)
  
  message(sprintf("  Processing %s (%d samples)...", tissue_name, nrow(expr_subset)))
  
  # Get top variable genes for this tissue
  gene_vars <- apply(expr_subset, 2, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(2000, ncol(expr_subset))])
  
  # PCA computation
  expr_scaled <- scale(expr_subset[, top_genes])
  pca_result <- prcomp(expr_scaled, center = FALSE, scale. = FALSE)
  
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    disease_status = meta_subset$disease_status,
    staging = meta_subset$staging,
    age = meta_subset$age
  )
  
  var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)
  
  # PCA Plot 1: Disease Status
  p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = disease_status)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = c("Control" = "navy", "AMD" = "firebrick")) +
    labs(
      title = sprintf("%s: PCA by Disease", tissue_name),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      color = "Disease"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(file.path(output_dir, sprintf("pca_%s_disease.pdf", tissue_slug)), p1, width = 6, height = 5)
  ggsave(file.path(output_dir, sprintf("pca_%s_disease.tiff", tissue_slug)), p1, width = 6, height = 5,
         dpi = 300, compression = "lzw")
  
  # PCA Plot 2: AMD Staging
  p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = staging)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = staging_colors) +
    labs(
      title = sprintf("%s: PCA by Staging", tissue_name),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      color = "Stage"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(file.path(output_dir, sprintf("pca_%s_staging.pdf", tissue_slug)), p2, width = 6, height = 5)
  ggsave(file.path(output_dir, sprintf("pca_%s_staging.tiff", tissue_slug)), p2, width = 6, height = 5,
         dpi = 300, compression = "lzw")
  
  # Contrastive PCA: Disease vs Control
  control_idx <- which(meta_subset$disease_status == "Control")
  amd_idx <- which(meta_subset$disease_status == "AMD")
  
  if (length(control_idx) >= 5 && length(amd_idx) >= 5) {
    cpca_result <- contrastive_pca(
      foreground = expr_subset[amd_idx, top_genes],
      background = expr_subset[control_idx, top_genes],
      n_components = 2,
      alpha = 1.0
    )
    
    cpca_df <- data.frame(
      cPC1 = cpca_result$projection[, 1],
      cPC2 = cpca_result$projection[, 2],
      disease_status = c(
        rep("Control", length(control_idx)),
        rep("AMD", length(amd_idx))
      ),
      staging = c(
        meta_subset$staging[control_idx],
        meta_subset$staging[amd_idx]
      )
    )
    
    cpca_var <- round(100 * cpca_result$var_explained, 1)
    
    # cPCA Plot 1: Disease
    p_cpca1 <- ggplot(cpca_df, aes(x = cPC1, y = cPC2, color = disease_status)) +
      geom_point(size = 2.5, alpha = 0.7) +
      scale_color_manual(values = c("Control" = "navy", "AMD" = "firebrick")) +
      labs(
        title = sprintf("%s: Contrastive PCA (Disease)", tissue_name),
        x = sprintf("cPC1 (%.1f%%)", cpca_var[1]),
        y = sprintf("cPC2 (%.1f%%)", cpca_var[2]),
        color = "Disease"
      ) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    ggsave(file.path(output_dir, sprintf("cpca_%s_disease.pdf", tissue_slug)), p_cpca1, width = 6, height = 5)
    ggsave(file.path(output_dir, sprintf("cpca_%s_disease.tiff", tissue_slug)), p_cpca1, width = 6, height = 5,
           dpi = 300, compression = "lzw")
    
    # cPCA Plot 2: Staging
    p_cpca2 <- ggplot(cpca_df, aes(x = cPC1, y = cPC2, color = staging)) +
      geom_point(size = 2.5, alpha = 0.7) +
      scale_color_manual(values = staging_colors) +
      labs(
        title = sprintf("%s: Contrastive PCA (Staging)", tissue_name),
        x = sprintf("cPC1 (%.1f%%)", cpca_var[1]),
        y = sprintf("cPC2 (%.1f%%)", cpca_var[2]),
        color = "Stage"
      ) +
      theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    ggsave(file.path(output_dir, sprintf("cpca_%s_staging.pdf", tissue_slug)), p_cpca2, width = 6, height = 5)
    ggsave(file.path(output_dir, sprintf("cpca_%s_staging.tiff", tissue_slug)), p_cpca2, width = 6, height = 5,
           dpi = 300, compression = "lzw")
  }
  
  # Heatmap: Top 500 variable genes
  top_hm_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(500, ncol(expr_subset))])
  heatmap_data <- expr_subset[, top_hm_genes]
  heatmap_zscore <- scale(heatmap_data)
  
  # Hierarchical clustering
  row_order <- hclust(dist(heatmap_zscore))$order
  
  # Annotation
  annotation_row <- data.frame(
    Disease = meta_subset$disease_status,
    Stage = meta_subset$staging,
    row.names = rownames(heatmap_data)
  )
  
  annotation_colors <- list(
    Disease = c("Control" = "navy", "AMD" = "firebrick"),
    Stage = staging_colors
  )
  
  # Generate heatmap PDF
  pdf(file.path(output_dir, sprintf("heatmap_%s.pdf", tissue_slug)), width = 10, height = 7)
  pheatmap(
    t(heatmap_zscore[row_order, ]),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = annotation_row[row_order, ],
    annotation_colors = annotation_colors,
    color = heatmap_palette,
    breaks = seq(-3, 3, length.out = 101),
    main = sprintf("%s: Top 500 Variable Genes", tissue_name),
    fontsize = 10,
    treeheight_row = 0,
    treeheight_col = 0
  )
  dev.off()
  
  # Generate heatmap TIFF
  tiff(file.path(output_dir, sprintf("heatmap_%s.tiff", tissue_slug)), width = 10, height = 7,
       units = "in", res = 300, compression = "lzw")
  pheatmap(
    t(heatmap_zscore[row_order, ]),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    annotation_col = annotation_row[row_order, ],
    annotation_colors = annotation_colors,
    color = heatmap_palette,
    breaks = seq(-3, 3, length.out = 101),
    main = sprintf("%s: Top 500 Variable Genes", tissue_name),
    fontsize = 10,
    treeheight_row = 0,
    treeheight_col = 0
  )
  dev.off()
  
  return(TRUE)
}

#%%
# Generate plots for each tissue type
message("\n=== Generating Tissue-Specific Visualizations ===")

unique_tissues <- unique(metadata$tissue)
message(sprintf("Found %d tissue types: %s", length(unique_tissues), paste(unique_tissues, collapse = ", ")))

for (tiss in unique_tissues) {
  # Subset data
  idx <- which(metadata$tissue == tiss)
  expr_subset <- expr_data[idx, ]
  meta_subset <- metadata[idx, ]
  
  generate_tissue_plots(expr_subset, meta_subset, tiss, OUTPUT_DIR)
}

#%%
# Summary
message("\n=== Visualization Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
message("Generated files per tissue:")
for (tiss in unique_tissues) {
  slug <- gsub("-", "_", gsub(" ", "_", tolower(tiss)))
  message(sprintf("  %s:", tiss))
  message(sprintf("    - pca_%s_disease.pdf / .tiff", slug))
  message(sprintf("    - pca_%s_staging.pdf / .tiff", slug))
  message(sprintf("    - cpca_%s_disease.pdf / .tiff (contrastive)", slug))
  message(sprintf("    - cpca_%s_staging.pdf / .tiff (contrastive)", slug))
  message(sprintf("    - heatmap_%s.pdf / .tiff", slug))
}


