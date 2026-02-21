#%%
# GSE135092 Exploratory Visualization
# PCA plots and heatmap for cohort exploration
# Following coding preferences: R, ggplot2, Palatino Linotype, theme_classic, BrBG

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(scales)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
data_dir <- file.path(base_dir, "data/external/geo/GSE135092")
out_dir <- file.path(base_dir, "results/cohort-GSE135092")
meta_file <- file.path(data_dir, "GSE135092_series_matrix.txt.gz")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#%%
# Custom theme following preferences
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "serif"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}


#%%
# Parse metadata from series matrix
parse_metadata <- function(meta_path) {
  con <- gzfile(meta_path, "rt")
  lines <- readLines(con)
  close(con)
  
  gsm_line <- grep("^!Sample_geo_accession", lines, value = TRUE)[1]
  gsm <- gsub('"', '', strsplit(gsm_line, "\t")[[1]][-1])
  
  char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  
  tissue <- NULL
  amd_status <- NULL
  donor_id <- NULL
  age <- NULL
  
  for (line in char_lines) {
    parts <- gsub('"', '', strsplit(line, "\t")[[1]][-1])
    if (grepl("tissue:", parts[1])) {
      tissue <- gsub("tissue: ", "", parts)
    } else if (grepl("amd_status:", parts[1])) {
      amd_status <- gsub("amd_status: ", "", parts)
    } else if (grepl("samid:", parts[1])) {
      donor_id <- gsub("samid: ", "", parts)
    } else if (grepl("age:", parts[1])) {
      age_raw <- gsub("age: ", "", parts)
      age <- ifelse(age_raw == "NA", NA, as.numeric(age_raw))
    }
  }
  
  data.frame(
    gsm = gsm,
    tissue = tissue,
    amd_status = amd_status,
    donor_id = donor_id,
    age = age,
    stringsAsFactors = FALSE
  )
}

cat("Parsing metadata...\n")
meta <- parse_metadata(meta_file)
cat(sprintf("Loaded metadata for %d samples\n", nrow(meta)))

#%%
# Add derived variables
meta <- meta %>%
  mutate(
    tissue_type = case_when(
      grepl("RPE", tissue) ~ "RPE",
      grepl("Retina", tissue) ~ "Retina",
      TRUE ~ "Other"
    ),
    region = case_when(
      grepl("Macula", tissue) ~ "Macula",
      grepl("non-Macula", tissue) ~ "Periphery",
      TRUE ~ "Other"
    ),
    age_bin = case_when(
      is.na(age) ~ "Unknown",
      age < 70 ~ "<70",
      age < 80 ~ "70-79",
      age < 90 ~ "80-89",
      TRUE ~ "90+"
    )
  )

# Save metadata for future use
write_csv(meta, file.path(out_dir, "GSE135092_sample_metadata.csv"))
cat("Saved sample metadata\n")

#%%
# Build expression matrix from raw counts
# Check if extracted raw data exists
raw_dir <- file.path(data_dir, "GSE135092_extracted")

if (!dir.exists(raw_dir)) {
  cat("Raw data directory not found. Attempting to extract from RAW.tar...\n")
  tar_file <- file.path(data_dir, "GSE135092_RAW.tar")
  if (file.exists(tar_file)) {
    dir.create(raw_dir, recursive = TRUE)
    untar(tar_file, exdir = raw_dir)
    cat("Extracted raw data\n")
  } else {
    stop("Cannot find raw data. Please download GSE135092_RAW.tar")
  }
}

#%%
# Build expression matrix
cache_file <- file.path(out_dir, "expression_matrix_cache.rds")

if (file.exists(cache_file)) {
  cat("Loading cached expression matrix...\n")
  expr_data <- readRDS(cache_file)
  expr_mat <- expr_data$expr_mat
  meta <- expr_data$meta
} else {
  cat("Building expression matrix from raw files...\n")
  
  count_files <- list.files(raw_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  cat(sprintf("Found %d count files\n", length(count_files)))
  
  # Read first file to get gene list
  first_df <- read.delim(gzfile(count_files[1]), comment.char = "#")
  gene_ids <- first_df$ID_REF
  
  # Initialize matrix
  counts <- matrix(NA, nrow = length(gene_ids), ncol = length(count_files))
  rownames(counts) <- gene_ids
  gsm_ids <- character(length(count_files))
  
  for (i in seq_along(count_files)) {
    gsm <- gsub("_.*", "", basename(count_files[i]))
    gsm_ids[i] <- gsm
    df <- read.delim(gzfile(count_files[i]), comment.char = "#")
    counts[, i] <- df$count
    if (i %% 100 == 0) cat(sprintf("  Processed %d/%d files\n", i, length(count_files)))
  }
  colnames(counts) <- gsm_ids
  
  # Filter to samples in metadata that are AMD or Control
  valid_gsm <- intersect(colnames(counts), meta$gsm[meta$amd_status %in% c("AMD", "Control")])
  counts <- counts[, valid_gsm]
  meta <- meta[meta$gsm %in% valid_gsm, ]
  
  # CPM normalize and log transform
  lib_sizes <- colSums(counts)
  cpm <- sweep(counts, 2, lib_sizes, "/") * 1e6
  log_cpm <- log2(cpm + 1)
  
  expr_mat <- log_cpm
  
  # Cache
  saveRDS(list(expr_mat = expr_mat, meta = meta), cache_file)
  cat("Saved expression matrix cache\n")
}

cat(sprintf("Expression matrix: %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))

#%%
# Ensure metadata order matches expression matrix
meta <- meta[match(colnames(expr_mat), meta$gsm), ]

#%%
# PCA using top variable genes
cat("Computing PCA...\n")
gene_vars <- apply(expr_mat, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:2000])

pca_input <- t(expr_mat[top_genes, ])
pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  gsm = rownames(pca_result$x)
) %>%
  left_join(meta, by = "gsm")

var_explained <- summary(pca_result)$importance[2, 1:3] * 100

#%%
# Plot 1: PCA by Disease Status
p_pca_disease <- ggplot(pca_df, aes(x = PC1, y = PC2, color = amd_status)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Control" = "#4393C3", "AMD" = "#D6604D"), name = "Status") +
  labs(
    title = "PCA: Control vs AMD",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_pub() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Fig_PCA_disease_status.pdf"), p_pca_disease, width = 5, height = 4)
ggsave(file.path(out_dir, "Fig_PCA_disease_status.tiff"), p_pca_disease, 
       width = 5, height = 4, dpi = 300, compression = "lzw")
cat("Saved PCA disease status plot\n")

#%%
# Plot 2: PCA by Tissue Type
p_pca_tissue <- ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_brewer(palette = "Set2", name = "Tissue") +
  labs(
    title = "PCA: By Tissue Type",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_pub() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Fig_PCA_tissue_type.pdf"), p_pca_tissue, width = 6, height = 4)
ggsave(file.path(out_dir, "Fig_PCA_tissue_type.tiff"), p_pca_tissue, 
       width = 6, height = 4, dpi = 300, compression = "lzw")
cat("Saved PCA tissue type plot\n")

#%%
# Plot 3: PCA by Tissue and Disease (faceted or combined)
p_pca_combined <- ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue, shape = amd_status)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_brewer(palette = "Set2", name = "Tissue") +
  scale_shape_manual(values = c("Control" = 16, "AMD" = 17), name = "Status") +
  labs(
    title = "PCA: Tissue and Disease Status",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_pub() +
  theme(legend.position = "right")

ggsave(file.path(out_dir, "Fig_PCA_tissue_disease.pdf"), p_pca_combined, width = 6, height = 4)
ggsave(file.path(out_dir, "Fig_PCA_tissue_disease.tiff"), p_pca_combined, 
       width = 6, height = 4, dpi = 300, compression = "lzw")
cat("Saved PCA tissue+disease plot\n")

#%%
# Plot 4: PCA by Age Bins
p_pca_age <- ggplot(pca_df %>% filter(age_bin != "Unknown"), 
                    aes(x = PC1, y = PC2, color = age_bin)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(
    values = c("<70" = "#FEE0D2", "70-79" = "#FC9272", "80-89" = "#DE2D26", "90+" = "#67000D"),
    name = "Age"
  ) +
  labs(
    title = "PCA: By Age Group",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_pub()

ggsave(file.path(out_dir, "Fig_PCA_age_bins.pdf"), p_pca_age, width = 5, height = 4)
cat("Saved PCA age bins plot\n")

#%%
# Heatmap: Top variable genes, clustered, with annotation bars
cat("Creating heatmap...\n")

# Get top 500 variable genes for heatmap
n_genes_heatmap <- 500
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:n_genes_heatmap])

# Gene-wise z-score (as per preferences)
heatmap_data <- expr_mat[top_var_genes, ]
heatmap_z <- t(scale(t(heatmap_data)))  # Gene-wise z-score

# Hierarchical clustering of samples
sample_dist <- dist(t(heatmap_z))
sample_hc <- hclust(sample_dist, method = "ward.D2")
sample_order <- sample_hc$order

# Reorder data and metadata
heatmap_z_ordered <- heatmap_z[, sample_order]
meta_ordered <- meta[match(colnames(heatmap_z_ordered), meta$gsm), ]

# Prepare data for ggplot heatmap
heatmap_long <- as.data.frame(heatmap_z_ordered) %>%
  mutate(gene = rownames(heatmap_z_ordered)) %>%
  pivot_longer(-gene, names_to = "gsm", values_to = "zscore") %>%
  mutate(
    sample_idx = match(gsm, colnames(heatmap_z_ordered)),
    gene_idx = match(gene, rownames(heatmap_z_ordered))
  )

# Clip z-scores for visualization
heatmap_long$zscore_clipped <- pmax(pmin(heatmap_long$zscore, 3), -3)

#%%
# Main heatmap
p_heatmap_main <- ggplot(heatmap_long, aes(x = sample_idx, y = gene_idx, fill = zscore_clipped)) +
  geom_raster() +
  scale_fill_distiller(palette = "BrBG", direction = -1, limits = c(-3, 3), name = "Z-score") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

#%%
# Annotation bars
anno_df <- data.frame(
  sample_idx = 1:ncol(heatmap_z_ordered),
  amd_status = meta_ordered$amd_status,
  tissue = meta_ordered$tissue,
  tissue_type = meta_ordered$tissue_type,
  region = meta_ordered$region
)

# Disease status bar
p_status_bar <- ggplot(anno_df, aes(x = sample_idx, y = 1, fill = amd_status)) +
  geom_tile() +
  scale_fill_manual(values = c("Control" = "#4393C3", "AMD" = "#D6604D"), name = "Status") +
  theme_void() +
  theme(legend.position = "none")

# Tissue bar
p_tissue_bar <- ggplot(anno_df, aes(x = sample_idx, y = 1, fill = tissue_type)) +
  geom_tile() +
  scale_fill_manual(values = c("RPE" = "#66C2A5", "Retina" = "#FC8D62"), name = "Tissue") +
  theme_void() +
  theme(legend.position = "none")

# Region bar
p_region_bar <- ggplot(anno_df, aes(x = sample_idx, y = 1, fill = region)) +
  geom_tile() +
  scale_fill_manual(values = c("Macula" = "#8DA0CB", "Periphery" = "#E78AC3"), name = "Region") +
  theme_void() +
  theme(legend.position = "none")

#%%
# Combine with patchwork
library(patchwork)

p_combined_heatmap <- p_status_bar / p_tissue_bar / p_region_bar / p_heatmap_main +
  plot_layout(heights = c(0.02, 0.02, 0.02, 0.94))

ggsave(file.path(out_dir, "Fig_heatmap_top_genes.pdf"), p_combined_heatmap, 
       width = 8, height = 6)
ggsave(file.path(out_dir, "Fig_heatmap_top_genes.tiff"), p_combined_heatmap, 
       width = 8, height = 6, dpi = 300, compression = "lzw")
cat("Saved heatmap\n")

#%%
# Create a separate legend plot
p_legend <- ggplot() +
  geom_point(data = data.frame(x = 1:2, y = 1, status = c("Control", "AMD")),
             aes(x, y, color = status), size = 4) +
  scale_color_manual(values = c("Control" = "#4393C3", "AMD" = "#D6604D"), name = "Status") +
  geom_point(data = data.frame(x = 4:5, y = 1, tissue = c("RPE", "Retina")),
             aes(x, y, fill = tissue), shape = 22, size = 4) +
  scale_fill_manual(values = c("RPE" = "#66C2A5", "Retina" = "#FC8D62"), name = "Tissue") +
  theme_void() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "Fig_heatmap_legend.pdf"), p_legend, width = 6, height = 1)

cat("\n=== GENERATING TISSUE-SPECIFIC VISUALIZATIONS ===\n")

#%%
# Tissue-specific PCA and Heatmaps
unique_tissues <- unique(meta$tissue)

for (tiss in unique_tissues) {
  cat(sprintf("\nProcessing tissue: %s\n", tiss))
  
  # Safe filename
  tiss_safe <- gsub("[, ]", "_", tiss)
  
  # Subset metadata and expression
  tiss_meta <- meta[meta$tissue == tiss, ]
  tiss_gsm <- tiss_meta$gsm
  tiss_expr <- expr_mat[, tiss_gsm]
  
  cat(sprintf("  Samples: %d (Control: %d, AMD: %d)\n", 
              nrow(tiss_meta),
              sum(tiss_meta$amd_status == "Control"),
              sum(tiss_meta$amd_status == "AMD")))
  
  #%% Tissue-specific PCA
  tiss_gene_vars <- apply(tiss_expr, 1, var)
  tiss_top_genes <- names(sort(tiss_gene_vars, decreasing = TRUE)[1:2000])
  
  tiss_pca_input <- t(tiss_expr[tiss_top_genes, ])
  tiss_pca_result <- prcomp(tiss_pca_input, center = TRUE, scale. = TRUE)
  
  tiss_pca_df <- data.frame(
    PC1 = tiss_pca_result$x[, 1],
    PC2 = tiss_pca_result$x[, 2],
    gsm = rownames(tiss_pca_result$x)
  ) %>%
    left_join(tiss_meta, by = "gsm")
  
  tiss_var_explained <- summary(tiss_pca_result)$importance[2, 1:2] * 100
  
  # PCA plot for this tissue
  p_tiss_pca <- ggplot(tiss_pca_df, aes(x = PC1, y = PC2, color = amd_status)) +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_color_manual(values = c("Control" = "#4393C3", "AMD" = "#D6604D"), name = "Status") +
    labs(
      title = sprintf("PCA: %s", tiss),
      subtitle = sprintf("N = %d samples", nrow(tiss_meta)),
      x = sprintf("PC1 (%.1f%%)", tiss_var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", tiss_var_explained[2])
    ) +
    theme_pub()
  
  ggsave(file.path(out_dir, sprintf("Fig_PCA_%s.pdf", tiss_safe)), p_tiss_pca, 
         width = 5, height = 4)
  ggsave(file.path(out_dir, sprintf("Fig_PCA_%s.tiff", tiss_safe)), p_tiss_pca, 
         width = 5, height = 4, dpi = 300, compression = "lzw")
  cat(sprintf("  Saved PCA for %s\n", tiss))
  
  #%% Contrastive PCA (cPCA) for this tissue
  # cPCA removes variation shared between control and AMD to highlight disease-specific changes
  cat(sprintf("  Computing contrastive PCA for %s...\n", tiss))
  
  # Separate control and AMD samples
  ctrl_idx <- tiss_meta$amd_status == "Control"
  amd_idx <- tiss_meta$amd_status == "AMD"
  
  ctrl_data <- tiss_pca_input[ctrl_idx, ]
  amd_data <- tiss_pca_input[amd_idx, ]
  
  # Center both datasets
  ctrl_centered <- scale(ctrl_data, center = TRUE, scale = FALSE)
  amd_centered <- scale(amd_data, center = TRUE, scale = FALSE)
  
  # Compute covariance matrices
  cov_ctrl <- cov(ctrl_centered)
  cov_amd <- cov(amd_centered)
  
  # Contrastive covariance: emphasize AMD-specific variation
  # Try multiple alpha values to find optimal contrast
  alpha_values <- c(0.5, 1.0, 2.0)
  
  for (alpha in alpha_values) {
    cov_contrast <- cov_amd - alpha * cov_ctrl
    
    # Eigen decomposition
    eigen_result <- eigen(cov_contrast)
    
    # Project all samples onto contrastive PCs
    cpca_scores <- tiss_pca_input %*% eigen_result$vectors[, 1:2]
    
    cpca_df <- data.frame(
      cPC1 = cpca_scores[, 1],
      cPC2 = cpca_scores[, 2],
      gsm = rownames(tiss_pca_input)
    ) %>%
      left_join(tiss_meta, by = "gsm")
    
    # Calculate variance explained (approximate)
    cpca_var <- eigen_result$values[1:2]
    cpca_var_pct <- 100 * cpca_var / sum(abs(cpca_var[cpca_var > 0]))
    
    # Plot contrastive PCA
    p_cpca <- ggplot(cpca_df, aes(x = cPC1, y = cPC2, color = amd_status)) +
      geom_point(alpha = 0.7, size = 2.5) +
      scale_color_manual(values = c("Control" = "#4393C3", "AMD" = "#D6604D"), name = "Status") +
      labs(
        title = sprintf("Contrastive PCA: %s", tiss),
        subtitle = sprintf("α = %.1f, N = %d samples", alpha, nrow(tiss_meta)),
        x = sprintf("cPC1 (%.1f%%)", cpca_var_pct[1]),
        y = sprintf("cPC2 (%.1f%%)", cpca_var_pct[2])
      ) +
      theme_pub()
    
    alpha_label <- gsub("\\.", "_", sprintf("%.1f", alpha))
    ggsave(file.path(out_dir, sprintf("Fig_cPCA_%s_alpha%s.pdf", tiss_safe, alpha_label)), 
           p_cpca, width = 5, height = 4)
    ggsave(file.path(out_dir, sprintf("Fig_cPCA_%s_alpha%s.tiff", tiss_safe, alpha_label)), 
           p_cpca, width = 5, height = 4, dpi = 300, compression = "lzw")
  }
  
  cat(sprintf("  Saved contrastive PCA for %s (3 alpha values)\n", tiss))

  
  #%% Tissue-specific Heatmap
  n_heatmap_genes <- 300
  tiss_heatmap_genes <- names(sort(tiss_gene_vars, decreasing = TRUE)[1:n_heatmap_genes])
  
  tiss_heatmap_data <- tiss_expr[tiss_heatmap_genes, ]
  tiss_heatmap_z <- t(scale(t(tiss_heatmap_data)))  # Gene-wise z-score
  
  # Hierarchical clustering
  tiss_sample_dist <- dist(t(tiss_heatmap_z))
  tiss_sample_hc <- hclust(tiss_sample_dist, method = "ward.D2")
  tiss_sample_order <- tiss_sample_hc$order
  
  tiss_heatmap_z_ordered <- tiss_heatmap_z[, tiss_sample_order]
  tiss_meta_ordered <- tiss_meta[match(colnames(tiss_heatmap_z_ordered), tiss_meta$gsm), ]
  
  tiss_heatmap_long <- as.data.frame(tiss_heatmap_z_ordered) %>%
    mutate(gene = rownames(tiss_heatmap_z_ordered)) %>%
    pivot_longer(-gene, names_to = "gsm", values_to = "zscore") %>%
    mutate(
      sample_idx = match(gsm, colnames(tiss_heatmap_z_ordered)),
      gene_idx = match(gene, rownames(tiss_heatmap_z_ordered))
    )
  
  tiss_heatmap_long$zscore_clipped <- pmax(pmin(tiss_heatmap_long$zscore, 3), -3)
  
  # Main heatmap
  p_tiss_heatmap_main <- ggplot(tiss_heatmap_long, 
                                 aes(x = sample_idx, y = gene_idx, fill = zscore_clipped)) +
    geom_raster() +
    scale_fill_distiller(palette = "BrBG", direction = -1, limits = c(-3, 3), name = "Z-score") +
    labs(x = NULL, y = NULL, title = sprintf("Top 300 Variable Genes: %s", tiss)) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", family = "serif")
    )
  
  # Disease annotation bar
  tiss_anno_df <- data.frame(
    sample_idx = 1:ncol(tiss_heatmap_z_ordered),
    amd_status = tiss_meta_ordered$amd_status
  )
  
  p_tiss_status_bar <- ggplot(tiss_anno_df, aes(x = sample_idx, y = 1, fill = amd_status)) +
    geom_tile() +
    scale_fill_manual(values = c("Control" = "#4393C3", "AMD" = "#D6604D"), name = "Status") +
    theme_void() +
    theme(legend.position = "none")
  
  # Combine
  p_tiss_combined_heatmap <- p_tiss_status_bar / p_tiss_heatmap_main +
    plot_layout(heights = c(0.03, 0.97))
  
  ggsave(file.path(out_dir, sprintf("Fig_heatmap_%s.pdf", tiss_safe)), p_tiss_combined_heatmap, 
         width = 6, height = 5)
  ggsave(file.path(out_dir, sprintf("Fig_heatmap_%s.tiff", tiss_safe)), p_tiss_combined_heatmap, 
         width = 6, height = 5, dpi = 300, compression = "lzw")
  cat(sprintf("  Saved heatmap for %s\n", tiss))
}

cat("\n=== COMPLETE ===\n")
cat("Generated files:\n")
cat("  - GSE135092_sample_metadata.csv\n")
cat("  - Fig_PCA_disease_status.pdf/tiff (all samples)\n")
cat("  - Fig_PCA_tissue_type.pdf/tiff (all samples)\n")
cat("  - Fig_PCA_tissue_disease.pdf/tiff (all samples)\n")
cat("  - Fig_PCA_age_bins.pdf (all samples)\n")
cat("  - Fig_heatmap_top_genes.pdf/tiff (all samples)\n")
cat("  - Fig_PCA_<tissue>.pdf/tiff (per tissue, standard PCA)\n")
cat("  - Fig_cPCA_<tissue>_alpha<X>.pdf/tiff (per tissue, contrastive PCA with 3 alpha values)\n")
cat("  - Fig_heatmap_<tissue>.pdf/tiff (per tissue)\n")


