#%%
# RPE-Specific Differential Expression Analysis with Age Covariate and Variance Testing
# GSE135092 - Focus on RPE samples only
# Compare: Mean DE, Variance testing, different thresholds

library(limma)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(car)  # For Levene test
library(VennDiagram)
library(patchwork)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092")
out_dir <- file.path(results_dir, "rpe_covariate_de")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load cached expression and metadata
cache_file <- file.path(results_dir, "expression_matrix_cache.rds")
expr_data <- readRDS(cache_file)
expr_mat <- expr_data$expr_mat
meta_all <- read_csv(file.path(results_dir, "GSE135092_sample_metadata.csv"))

cat("Loaded data:\n")
cat(sprintf("  Expression: %d genes x %d samples\n", nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("  Metadata: %d samples\n", nrow(meta_all)))

#%%
# Custom theme
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
# Function to run DE analysis with age covariate using limma
run_limma_de <- function(expr, meta, tissue_name) {
  cat(sprintf("\n=== Running limma DE for %s ===\n", tissue_name))
  
  # Filter samples
  meta_filt <- meta %>%
    filter(!is.na(age), amd_status %in% c("Control", "AMD"))
  
  expr_filt <- expr[, meta_filt$gsm]
  
  cat(sprintf("  Samples: %d total (%d Control, %d AMD)\n",
              nrow(meta_filt),
              sum(meta_filt$amd_status == "Control"),
              sum(meta_filt$amd_status == "AMD")))
  cat(sprintf("  Age range: %.0f-%.0f years\n", min(meta_filt$age), max(meta_filt$age)))
  
  # Design matrix with age covariate
  design <- model.matrix(~ age + amd_status, data = meta_filt)
  
  # Fit linear model
  fit <- lmFit(expr_filt, design)
  fit <- eBayes(fit)
  
  # Extract results for AMD vs Control (amd_statusControl)
  res_amd <- topTable(fit, coef = "amd_statusControl", number = Inf, sort.by = "none")
  res_amd$gene <- rownames(res_amd)
  
  # Extract results for Age (age)
  res_age <- topTable(fit, coef = "age", number = Inf, sort.by = "none")
  res_age$gene <- rownames(res_age)
  
  return(list(amd = res_amd, age = res_age))
}

#%%
# Function to test for differential variance (Levene/Bartlett approach)
run_variance_test <- function(expr, meta, tissue_name) {
  cat(sprintf("\n=== Running variance test for %s ===\n", tissue_name))
  
  # Filter samples
  meta_filt <- meta %>%
    filter(!is.na(age), amd_status %in% c("Control", "AMD"))
  
  expr_filt <- expr[, meta_filt$gsm]
  
  cat(sprintf("  Testing variance differences for %d genes\n", nrow(expr_filt)))
  
  # For each gene, test if variance differs between AMD and Control
  # Using Levene's test (more robust than Bartlett for non-normal data)
  results <- data.frame(
    gene = rownames(expr_filt),
    levene_pval = NA_real_,
    ctrl_var = NA_real_,
    amd_var = NA_real_,
    var_ratio = NA_real_
  )
  
  ctrl_idx <- meta_filt$amd_status == "Control"
  amd_idx <- meta_filt$amd_status == "AMD"
  
  for (i in 1:nrow(expr_filt)) {
    gene_expr <- as.numeric(expr_filt[i, ])
    group <- factor(meta_filt$amd_status)
    
    # Levene test
    tryCatch({
      levene_result <- leveneTest(gene_expr ~ group)
      results$levene_pval[i] <- levene_result$`Pr(>F)`[1]
    }, error = function(e) {
      results$levene_pval[i] <- NA
    })
    
    # Calculate variances
    results$ctrl_var[i] <- var(gene_expr[ctrl_idx])
    results$amd_var[i] <- var(gene_expr[amd_idx])
    results$var_ratio[i] <- results$amd_var[i] / results$ctrl_var[i]
    
    if (i %% 5000 == 0) cat(sprintf("  Processed %d/%d genes\n", i, nrow(expr_filt)))
  }
  
  # Calculate FDR
  results$levene_fdr <- p.adjust(results$levene_pval, method = "BH")
  
  # Add direction (increased vs decreased variance in AMD)
  results$variance_direction <- ifelse(results$var_ratio > 1, "Increased", "Decreased")
  
  return(results)
}

#%%
# Function to run combined variance-mean test (residual variance after accounting for mean)
run_residual_variance_test <- function(expr, meta, tissue_name) {
  cat(sprintf("\n=== Running residual variance test for %s ===\n", tissue_name))
  
  # Filter samples
  meta_filt <- meta %>%
    filter(!is.na(age), amd_status %in% c("Control", "AMD"))
  
  expr_filt <- expr[, meta_filt$gsm]
  
  # For each gene, fit model with age, extract residuals, test variance of residuals
  results <- data.frame(
    gene = rownames(expr_filt),
    resid_var_pval = NA_real_,
    ctrl_resid_var = NA_real_,
    amd_resid_var = NA_real_
  )
  
  ctrl_idx <- meta_filt$amd_status == "Control"
  amd_idx <- meta_filt$amd_status == "AMD"
  
  for (i in 1:nrow(expr_filt)) {
    gene_expr <- as.numeric(expr_filt[i, ])
    
    # Fit model: expression ~ age
    fit <- lm(gene_expr ~ meta_filt$age)
    residuals <- residuals(fit)
    
    # Test variance of residuals between groups
    group <- factor(meta_filt$amd_status)
    tryCatch({
      levene_result <- leveneTest(residuals ~ group)
      results$resid_var_pval[i] <- levene_result$`Pr(>F)`[1]
    }, error = function(e) {
      results$resid_var_pval[i] <- NA
    })
    
    results$ctrl_resid_var[i] <- var(residuals[ctrl_idx])
    results$amd_resid_var[i] <- var(residuals[amd_idx])
    
    if (i %% 5000 == 0) cat(sprintf("  Processed %d/%d genes\n", i, nrow(expr_filt)))
  }
  
  results$resid_var_fdr <- p.adjust(results$resid_var_pval, method = "BH")
  
  return(results)
}

#%%
# Analyze RPE Macula
cat("\n" , rep("=", 70), "\n", sep = "")
cat("ANALYZING RPE MACULA\n")
cat(rep("=", 70), "\n", sep = "")

meta_rpe_mac <- meta_all %>% filter(tissue == "RPE, Macula")
expr_rpe_mac <- expr_mat[, meta_rpe_mac$gsm]

# Method 1: Standard DE with age covariate & Age DE extraction
de_results_mac <- run_limma_de(expr_rpe_mac, meta_rpe_mac, "RPE Macula")
de_rpe_mac <- de_results_mac$amd
age_rpe_mac <- de_results_mac$age

# Save Age Results
write_csv(age_rpe_mac, file.path(out_dir, "RPE_Macula_Age_DE_results.csv"))
cat("Saved: RPE_Macula_Age_DE_results.csv\n")

# Method 2: Variance testing
var_rpe_mac <- run_variance_test(expr_rpe_mac, meta_rpe_mac, "RPE Macula")

# Method 3: Residual variance testing
resid_var_rpe_mac <- run_residual_variance_test(expr_rpe_mac, meta_rpe_mac, "RPE Macula")

# Combine results (AMD focus for main file)
results_rpe_mac <- de_rpe_mac %>%
  left_join(var_rpe_mac, by = "gene") %>%
  left_join(resid_var_rpe_mac, by = "gene")

write_csv(results_rpe_mac, file.path(out_dir, "RPE_Macula_all_results.csv"))
cat("Saved: RPE_Macula_all_results.csv\n")

#%%
# Analyze RPE non-Macula
cat("\n", rep("=", 70), "\n", sep = "")
cat("ANALYZING RPE NON-MACULA\n")
cat(rep("=", 70), "\n", sep = "")

meta_rpe_nonmac <- meta_all %>% filter(tissue == "RPE, non-Macula")
expr_rpe_nonmac <- expr_mat[, meta_rpe_nonmac$gsm]

# Method 1: Standard DE with age covariate
de_results_nonmac <- run_limma_de(expr_rpe_nonmac, meta_rpe_nonmac, "RPE non-Macula")
de_rpe_nonmac <- de_results_nonmac$amd
age_rpe_nonmac <- de_results_nonmac$age

# Save Age Results
write_csv(age_rpe_nonmac, file.path(out_dir, "RPE_nonMacula_Age_DE_results.csv"))
cat("Saved: RPE_nonMacula_Age_DE_results.csv\n")

# Method 2: Variance testing
var_rpe_nonmac <- run_variance_test(expr_rpe_nonmac, meta_rpe_nonmac, "RPE non-Macula")

# Method 3: Residual variance testing
resid_var_rpe_nonmac <- run_residual_variance_test(expr_rpe_nonmac, meta_rpe_nonmac, "RPE non-Macula")

# Combine results
results_rpe_nonmac <- de_rpe_nonmac %>%
  left_join(var_rpe_nonmac, by = "gene") %>%
  left_join(resid_var_rpe_nonmac, by = "gene")

write_csv(results_rpe_nonmac, file.path(out_dir, "RPE_nonMacula_all_results.csv"))
cat("Saved: RPE_nonMacula_all_results.csv\n")

#%%
# Count significant genes for each method and threshold
count_sig_genes <- function(results_df, tissue_label) {
  # Initialize with all rows at once
  counts <- data.frame(
    Tissue = rep(tissue_label, 6),
    Method = c(
      "Mean DE (limma)", "Mean DE (limma)",
      "Variance (Levene)", "Variance (Levene)",
      "Residual Variance", "Residual Variance"
    ),
    Threshold = rep(c("p < 0.01", "FDR < 0.05"), 3),
    Count = c(
      sum(results_df$P.Value < 0.01, na.rm = TRUE),
      sum(results_df$adj.P.Val < 0.05, na.rm = TRUE),
      sum(results_df$levene_pval < 0.01, na.rm = TRUE),
      sum(results_df$levene_fdr < 0.05, na.rm = TRUE),
      sum(results_df$resid_var_pval < 0.01, na.rm = TRUE),
      sum(results_df$resid_var_fdr < 0.05, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  
  return(counts)
}


counts_mac <- count_sig_genes(results_rpe_mac, "RPE Macula")
counts_nonmac <- count_sig_genes(results_rpe_nonmac, "RPE non-Macula")
all_counts <- rbind(counts_mac, counts_nonmac)

write_csv(all_counts, file.path(out_dir, "gene_counts_summary.csv"))

#%%
# Create comparison plots
p_counts <- ggplot(all_counts, aes(x = Method, y = Count, fill = Threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Tissue, scales = "free_y") +
  scale_fill_manual(values = c("p < 0.01" = "#FC8D62", "FDR < 0.05" = "#8DA0CB")) +
  labs(
    title = "Significant Genes Detected by Method and Threshold",
    subtitle = "RPE samples with age covariate",
    x = NULL,
    y = "Number of Significant Genes"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "Fig_gene_counts_comparison.pdf"), p_counts, width = 8, height = 5)
ggsave(file.path(out_dir, "Fig_gene_counts_comparison.tiff"), p_counts, 
       width = 8, height = 5, dpi = 300, compression = "lzw")
cat("Saved comparison plot\n")

#%%
# Print summary
cat("\n", rep("=", 70), "\n", sep = "")
cat("SUMMARY OF SIGNIFICANT GENES\n")
cat(rep("=", 70), "\n\n")
print(all_counts)

cat("\n", rep("=", 70), "\n", sep = "")
cat("FILES SAVED:\n")
cat(rep("=", 70), "\n")
cat("  - RPE_Macula_all_results.csv\n")
cat("  - RPE_nonMacula_all_results.csv\n")
cat("  - gene_counts_summary.csv\n")
cat("  - Fig_gene_counts_comparison.pdf/tiff\n")
cat("\nAll files in:", out_dir, "\n")
