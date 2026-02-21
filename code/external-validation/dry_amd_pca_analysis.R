# GSE29801 Dry AMD PCA/cPCA Visualization Script
# Generates PCA and contrastive PCA plots for dry AMD stages

#%%
# Load libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

#%%
# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- DATA_DIR

# Visualization settings
FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

#%%
# Color definitions
stage_colors <- c(
  "Control" = "navy",
  "MD1" = "#fee391",
  "MD2" = "#fe9929",
  "dry_AMD" = "#cc4c02",
  "GA" = "#662506"
)

#%%
# Contrastive PCA implementation
contrastive_pca <- function(foreground, background, n_components = 2, alpha = 1.0) {
  fg_centered <- scale(foreground, center = TRUE, scale = FALSE)
  bg_centered <- scale(background, center = TRUE, scale = FALSE)
  
  cov_fg <- cov(fg_centered)
  cov_bg <- cov(bg_centered)
  cov_contrastive <- cov_fg - alpha * cov_bg
  
  eigen_result <- eigen(cov_contrastive)
  loadings <- eigen_result$vectors[, 1:n_components, drop = FALSE]
  eigenvalues <- eigen_result$values[1:n_components]
  
  fg_projected <- fg_centered %*% loadings
  bg_projected <- bg_centered %*% loadings
  all_projected <- rbind(bg_projected, fg_projected)
  
  total_var <- sum(abs(eigenvalues))
  var_explained <- abs(eigenvalues) / total_var
  
  return(list(
    projection = all_projected,
    var_explained = var_explained
  ))
}

#%%
# Function to generate PCA/cPCA for a tissue
generate_pca_plots <- function(data_list, tissue_name) {
  
  expr <- data_list$expr
  metadata <- data_list$metadata
  
  tissue_slug <- gsub(" ", "_", tolower(tissue_name))
  
  message(sprintf("\nProcessing %s (%d samples)...", tissue_name, nrow(expr)))
  
  # Get top variable genes
  gene_vars <- apply(expr, 2, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(2000, ncol(expr))])
  
  # ===== 1. PCA/cPCA with all dry stages =====
  message("  Generating all-stages PCA/cPCA...")
  
  # Standard PCA
  expr_scaled <- scale(expr[, top_genes])
  pca_result <- prcomp(expr_scaled, center = FALSE, scale. = FALSE)
  
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    stage = metadata$dry_stage
  )
  
  var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)
  
  p_pca_all <- ggplot(pca_df, aes(x = PC1, y = PC2, color = stage)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = stage_colors) +
    labs(
      title = sprintf("%s: PCA (All Dry Stages)", tissue_name),
      x = sprintf("PC1 (%.1f%%)", var_explained[1]),
      y = sprintf("PC2 (%.1f%%)", var_explained[2]),
      color = "Stage"
    ) +
    theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(file.path(OUTPUT_DIR, sprintf("pca_%s_all_stages.pdf", tissue_slug)), 
         p_pca_all, width = 6, height = 5)
  
  # Contrastive PCA (all AMD vs Control)
  control_idx <- which(metadata$disease_status == "Control")
  amd_idx <- which(metadata$disease_status == "AMD")
  
  if (length(control_idx) >= 5 && length(amd_idx) >= 5) {
    cpca_result <- contrastive_pca(
      foreground = expr[amd_idx, top_genes],
      background = expr[control_idx, top_genes],
      n_components = 2
    )
    
    cpca_df <- data.frame(
      cPC1 = cpca_result$projection[, 1],
      cPC2 = cpca_result$projection[, 2],
      stage = c(metadata$dry_stage[control_idx], metadata$dry_stage[amd_idx])
    )
    
    cpca_var <- round(100 * cpca_result$var_explained, 1)
    
    p_cpca_all <- ggplot(cpca_df, aes(x = cPC1, y = cPC2, color = stage)) +
      geom_point(size = 2.5, alpha = 0.7) +
      scale_color_manual(values = stage_colors) +
      labs(
        title = sprintf("%s: cPCA (All Dry Stages)", tissue_name),
        x = sprintf("cPC1 (%.1f%%)", cpca_var[1]),
        y = sprintf("cPC2 (%.1f%%)", cpca_var[2]),
        color = "Stage"
      ) +
      theme(legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(file.path(OUTPUT_DIR, sprintf("cpca_%s_all_stages.pdf", tissue_slug)), 
           p_cpca_all, width = 6, height = 5)
  }
  
  # ===== 2. Individual stage analyses (panel figure) =====
  message("  Generating individual stage PCA/cPCA panels...")
  
  stages_to_test <- c("MD1", "MD2", "dry_AMD", "GA")
  panel_plots <- list()
  
  for (stage in stages_to_test) {
    stage_idx <- which(metadata$dry_stage == stage)
    
    if (length(stage_idx) < 3) {
      message(sprintf("    Skipping %s: only %d samples", stage, length(stage_idx)))
      next
    }
    
    # Subset to control + this stage only
    subset_idx <- c(control_idx, stage_idx)
    subset_expr <- expr[subset_idx, top_genes]
    subset_meta <- metadata[subset_idx, ]
    
    # PCA
    subset_scaled <- scale(subset_expr)
    subset_pca <- prcomp(subset_scaled, center = FALSE, scale. = FALSE)
    
    pca_df_sub <- data.frame(
      PC1 = subset_pca$x[, 1],
      PC2 = subset_pca$x[, 2],
      stage = subset_meta$dry_stage
    )
    
    var_exp_sub <- round(100 * summary(subset_pca)$importance[2, 1:2], 1)
    
    p_pca_sub <- ggplot(pca_df_sub, aes(x = PC1, y = PC2, color = stage)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = stage_colors[c("Control", stage)]) +
      labs(
        title = sprintf("PCA: Control vs %s", stage),
        x = sprintf("PC1 (%.1f%%)", var_exp_sub[1]),
        y = sprintf("PC2 (%.1f%%)", var_exp_sub[2])
      ) +
      theme(legend.position = "bottom", plot.title = element_text(size = 10, face = "bold"))
    
    # cPCA
    cpca_sub <- contrastive_pca(
      foreground = subset_expr[subset_meta$dry_stage == stage, ],
      background = subset_expr[subset_meta$dry_stage == "Control", ],
      n_components = 2
    )
    
    cpca_df_sub <- data.frame(
      cPC1 = cpca_sub$projection[, 1],
      cPC2 = cpca_sub$projection[, 2],
      stage = subset_meta$dry_stage
    )
    
    cpca_var_sub <- round(100 * cpca_sub$var_explained, 1)
    
    p_cpca_sub <- ggplot(cpca_df_sub, aes(x = cPC1, y = cPC2, color = stage)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = stage_colors[c("Control", stage)]) +
      labs(
        title = sprintf("cPCA: Control vs %s", stage),
        x = sprintf("cPC1 (%.1f%%)", cpca_var_sub[1]),
        y = sprintf("cPC2 (%.1f%%)", cpca_var_sub[2])
      ) +
      theme(legend.position = "bottom", plot.title = element_text(size = 10, face = "bold"))
    
    panel_plots[[stage]] <- list(pca = p_pca_sub, cpca = p_cpca_sub)
  }
  
  # Create panel figure
  if (length(panel_plots) > 0) {
    plot_list <- unlist(lapply(panel_plots, function(x) list(x$pca, x$cpca)), recursive = FALSE)
    
    pdf(file.path(OUTPUT_DIR, sprintf("pca_cpca_panel_%s.pdf", tissue_slug)), 
        width = 12, height = 10)
    grid.arrange(grobs = plot_list, ncol = 2)
    dev.off()
  }
  
  return(TRUE)
}

#%%
# Load and process macular data
message("=== Generating PCA/cPCA Visualizations ===")

macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
generate_pca_plots(macular_data, "Macular RPE-choroid")

#%%
# Load and process extramacular data
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))
generate_pca_plots(extramacular_data, "Extramacular RPE-choroid")

#%%
# Summary
message("\n=== PCA/cPCA Visualization Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
message("\nGenerated files:")
message("  Macular:")
message("    - pca_macular_rpe_choroid_all_stages.pdf")
message("    - cpca_macular_rpe_choroid_all_stages.pdf")
message("    - pca_cpca_panel_macular_rpe_choroid.pdf")
message("  Extramacular:")
message("    - pca_extramacular_rpe_choroid_all_stages.pdf")
message("    - cpca_extramacular_rpe_choroid_all_stages.pdf")
message("    - pca_cpca_panel_extramacular_rpe_choroid.pdf")
