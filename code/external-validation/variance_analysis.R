# Extended Dry AMD - Differential Variability Analysis
# Tests variance differences, not just mean differences

#%%
# Load libraries
library(limma)
library(dplyr)
library(ggplot2)

#%%
# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- DATA_DIR

FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

FDR_THRESH <- 0.05
P_THRESH <- 0.01

#%%
# Load data
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

#%%
# Function to test differential variability
test_variance <- function(expr, metadata, tissue_name) {
  
  message(sprintf("\nTesting %s...", tissue_name))
  
  control_idx <- which(metadata$disease_status == "Control")
  amd_idx <- which(metadata$disease_status == "AMD")
  
  results <- data.frame(
    gene = colnames(expr),
    mean_control = NA,
    mean_amd = NA,
    var_control = NA,
    var_amd = NA,
    cv_control = NA,
    cv_amd = NA,
    f_stat = NA,
    f_pvalue = NA,
    levene_pvalue = NA,
    extreme_freq_control = NA,
    extreme_freq_amd = NA
  )
  
  for (i in 1:ncol(expr)) {
    gene_expr <- expr[, i]
    
    control_expr <- gene_expr[control_idx]
    amd_expr <- gene_expr[amd_idx]
    
    # Remove NAs
    control_expr <- control_expr[!is.na(control_expr)]
    amd_expr <- amd_expr[!is.na(amd_expr)]
    
    if (length(control_expr) < 3 || length(amd_expr) < 3) next
    
    # Mean and variance
    results$mean_control[i] <- mean(control_expr)
    results$mean_amd[i] <- mean(amd_expr)
    results$var_control[i] <- var(control_expr)
    results$var_amd[i] <- var(amd_expr)
    
    # Coefficient of variation
    results$cv_control[i] <- sd(control_expr) / abs(mean(control_expr) + 1e-6)
    results$cv_amd[i] <- sd(amd_expr) / abs(mean(amd_expr) + 1e-6)
    
    # F-test for equal variances
    f_test <- tryCatch({
      var.test(amd_expr, control_expr)
    }, error = function(e) list(statistic = NA, p.value = NA))
    
    results$f_stat[i] <- as.numeric(f_test$statistic)
    results$f_pvalue[i] <- f_test$p.value
    
    # Levene's test (using median)
    levene_p <- tryCatch({
      group <- c(rep("Control", length(control_expr)), rep("AMD", length(amd_expr)))
      values <- c(control_expr, amd_expr)
      median_vals <- tapply(values, group, median)
      abs_dev <- abs(values - median_vals[group])
      anova(lm(abs_dev ~ group))$`Pr(>F)`[1]
    }, error = function(e) NA)
    
    results$levene_pvalue[i] <- levene_p
    
    # Extreme values (>2 SD from mean)
    pooled_sd <- sd(c(control_expr, amd_expr))
    pooled_mean <- mean(c(control_expr, amd_expr))
    
    results$extreme_freq_control[i] <- sum(abs(control_expr - pooled_mean) > 2 * pooled_sd) / length(control_expr)
    results$extreme_freq_amd[i] <- sum(abs(amd_expr - pooled_mean) > 2 * pooled_sd) / length(amd_expr)
  }
  
  # FDR correction
  results$f_fdr <- p.adjust(results$f_pvalue, method = "BH")
  results$levene_fdr <- p.adjust(results$levene_pvalue, method = "BH")
  
  # Variance fold change
  results$var_fc <- log2(results$var_amd / (results$var_control + 1e-6))
  
  # CV difference
  results$cv_diff <- results$cv_amd - results$cv_control
  
  return(results)
}

#%%
# Run variance analysis
macular_var <- test_variance(macular_data$expr, macular_data$metadata, "Macular")
extramacular_var <- test_variance(extramacular_data$expr, extramacular_data$metadata, "Extramacular")

# Save results
write.csv(macular_var, file.path(OUTPUT_DIR, "macular_variance_analysis.csv"), row.names = FALSE)
write.csv(extramacular_var, file.path(OUTPUT_DIR, "extramacular_variance_analysis.csv"), row.names = FALSE)

#%%
# Count differentially variable genes
count_var_genes <- function(var_results, tissue) {
  data.frame(
    tissue = tissue,
    f_test_fdr05 = sum(var_results$f_fdr < 0.05, na.rm = TRUE),
    f_test_p001 = sum(var_results$f_pvalue < 0.01, na.rm = TRUE),
    levene_fdr05 = sum(var_results$levene_fdr < 0.05, na.rm = TRUE),
    levene_p001 = sum(var_results$levene_pvalue < 0.01, na.rm = TRUE),
    high_cv_diff = sum(abs(var_results$cv_diff) > 0.5, na.rm = TRUE),
    high_var_fc = sum(abs(var_results$var_fc) > 1, na.rm = TRUE)
  )
}

var_counts <- rbind(
  count_var_genes(macular_var, "Macular"),
  count_var_genes(extramacular_var, "Extramacular")
)

write.csv(var_counts, file.path(OUTPUT_DIR, "variance_analysis_summary.csv"), row.names = FALSE)

message("\n=== Variance Analysis Summary ===")
print(var_counts)

#%%
# Visualization: Variance vs Mean plot
plot_variance_analysis <- function(var_results, tissue_name) {
  
  # Filter to genes with sufficient data
  plot_data <- var_results %>%
    filter(!is.na(f_pvalue) & !is.na(mean_control))
  
  plot_data$signif <- ifelse(plot_data$f_fdr < 0.05, "FDR < 0.05",
                             ifelse(plot_data$f_pvalue < 0.01, "p < 0.01", "NS"))
  
  p <- ggplot(plot_data, aes(x = (mean_control + mean_amd) / 2, y = var_fc, color = signif)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("FDR < 0.05" = "red", "p < 0.01" = "orange", "NS" = "gray70")) +
    labs(
      title = sprintf("%s: Differential Variability", tissue_name),
      x = "Average Expression",
      y = "Variance Fold Change (log2)",
      color = "Significance"
    ) +
    theme(legend.position = "right")
  
  return(p)
}

p_mac <- plot_variance_analysis(macular_var, "Macular")
p_extra <- plot_variance_analysis(extramacular_var, "Extramacular")

ggsave(file.path(OUTPUT_DIR, "macular_variance_plot.pdf"), p_mac, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "extramacular_variance_plot.pdf"), p_extra, width = 7, height = 5)

#%%
# Extreme value comparison
plot_extreme_values <- function(var_results, tissue_name) {
  
  plot_data <- var_results %>%
    filter(!is.na(extreme_freq_control) & !is.na(extreme_freq_amd)) %>%
    mutate(extreme_diff = extreme_freq_amd - extreme_freq_control)
  
  p <- ggplot(plot_data, aes(x = extreme_diff)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = sprintf("%s: Extreme Value Frequency Difference", tissue_name),
      x = "Extreme Frequency (AMD - Control)",
      y = "Number of Genes"
    )
  
  return(p)
}

p_extreme_mac <- plot_extreme_values(macular_var, "Macular")
p_extreme_extra <- plot_extreme_values(extramacular_var, "Extramacular")

ggsave(file.path(OUTPUT_DIR, "macular_extreme_values.pdf"), p_extreme_mac, width = 6, height = 4)
ggsave(file.path(OUTPUT_DIR, "extramacular_extreme_values.pdf"), p_extreme_extra, width = 6, height = 4)

message("\n=== Variance Analysis Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
