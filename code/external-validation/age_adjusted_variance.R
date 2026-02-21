# Age-Adjusted Variance Analysis
# Tests disease variance after removing age effect
# Also tests if age itself drives variance

#%%
library(dplyr)
library(ggplot2)

#%%
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- DATA_DIR

#%%
# Load data
message("Loading data...")
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))
macular_var <- read.csv(file.path(DATA_DIR, "macular_variance_analysis_annotated.csv"))
extramacular_var <- read.csv(file.path(DATA_DIR, "extramacular_variance_analysis_annotated.csv"))

#%%
# Function: Age-Adjusted Disease Variance Test
# Removes age effect first, then tests disease variance
test_age_adjusted_variance <- function(expr, metadata, tissue_name) {
  
  message(sprintf("\n=== Age-Adjusted Disease Variance: %s ===", tissue_name))
  
  results <- data.frame(
    gene = colnames(expr),
    f_stat_adj = NA,
    f_pvalue_adj = NA
  )
  
  for (i in 1:ncol(expr)) {
    gene_expr <- expr[, i]
    
    # Remove NA
    valid_idx <- !is.na(gene_expr) & !is.na(metadata$age)
    if (sum(valid_idx) < 10) next
    
    gene_expr_clean <- gene_expr[valid_idx]
    age_clean <- metadata$age[valid_idx]
    disease_clean <- metadata$disease_status[valid_idx]
    
    # Step 1: Remove age effect (regress out age)
    age_model <- lm(gene_expr_clean ~ age_clean)
    residuals <- residuals(age_model)
    
    # Step 2: Test variance on residuals
    control_resid <- residuals[disease_clean == "Control"]
    amd_resid <- residuals[disease_clean == "AMD"]
    
    if (length(control_resid) < 3 || length(amd_resid) < 3) next
    
    # F-test on residuals
    f_test <- tryCatch({
      var.test(amd_resid, control_resid)
    }, error = function(e) list(statistic = NA, p.value = NA))
    
    results$f_stat_adj[i] <- as.numeric(f_test$statistic)
    results$f_pvalue_adj[i] <- f_test$p.value
  }
  
  # FDR correction
  results$f_fdr_adj <- p.adjust(results$f_pvalue_adj, method = "BH")
  
  return(results)
}

#%%
# Function: Age-Driven Variance Test
# Tests if variance increases with age (continuous)
test_age_variance <- function(expr, metadata, tissue_name) {
  
  message(sprintf("\n=== Age-Driven Variance: %s ===", tissue_name))
  
  results <- data.frame(
    gene = colnames(expr),
    age_var_pvalue = NA
  )
  
  for (i in 1:ncol(expr)) {
    gene_expr <- expr[, i]
    
    # Remove NA
    valid_idx <- !is.na(gene_expr) & !is.na(metadata$age)
    if (sum(valid_idx) < 10) next
    
    gene_expr_clean <- gene_expr[valid_idx]
    age_clean <- metadata$age[valid_idx]
    
    # Calculate squared residuals from overall mean
    mean_expr <- mean(gene_expr_clean)
    sq_residuals <- (gene_expr_clean - mean_expr)^2
    
    # Test if variance depends on age (Breusch-Pagan test concept)
    # Model: squared_residuals ~ age
    tryCatch({
      age_var_model <- lm(sq_residuals ~ age_clean)
      age_coef_pvalue <- summary(age_var_model)$coefficients[2, 4]  # p-value for age coefficient
      results$age_var_pvalue[i] <- age_coef_pvalue
    }, error = function(e) {
      results$age_var_pvalue[i] <- NA
    })
  }
  
  # FDR correction
  results$age_var_fdr <- p.adjust(results$age_var_pvalue, method = "BH")
  
  return(results)
}

#%%
# Run age-adjusted disease variance
macular_var_adj <- test_age_adjusted_variance(macular_data$expr, macular_data$metadata, "Macular")
extramacular_var_adj <- test_age_adjusted_variance(extramacular_data$expr, extramacular_data$metadata, "Extramacular")

# Count significant
mac_adj_fdr05 <- sum(macular_var_adj$f_fdr_adj < 0.05, na.rm = TRUE)
mac_adj_p001 <- sum(macular_var_adj$f_pvalue_adj < 0.01, na.rm = TRUE)
extra_adj_fdr05 <- sum(extramacular_var_adj$f_fdr_adj < 0.05, na.rm = TRUE)
extra_adj_p001 <- sum(extramacular_var_adj$f_pvalue_adj < 0.01, na.rm = TRUE)

message(sprintf("\nAge-Adjusted Disease Variance Results:"))
message(sprintf("  Macular: %d (FDR<0.05), %d (p<0.01)", mac_adj_fdr05, mac_adj_p001))
message(sprintf("  Extramacular: %d (FDR<0.05), %d (p<0.01)", extra_adj_fdr05, extra_adj_p001))

# Save
write.csv(macular_var_adj, 
          file.path(OUTPUT_DIR, "macular_variance_age_adjusted.csv"),
          row.names = FALSE)
write.csv(extramacular_var_adj,
          file.path(OUTPUT_DIR, "extramacular_variance_age_adjusted.csv"),
          row.names = FALSE)

#%%
# Run age-driven variance
macular_age_var <- test_age_variance(macular_data$expr, macular_data$metadata, "Macular")
extramacular_age_var <- test_age_variance(extramacular_data$expr, extramacular_data$metadata, "Extramacular")

# Count significant
mac_agevar_fdr05 <- sum(macular_age_var$age_var_fdr < 0.05, na.rm = TRUE)
mac_agevar_p001 <- sum(macular_age_var$age_var_pvalue < 0.01, na.rm = TRUE)
extra_agevar_fdr05 <- sum(extramacular_age_var$age_var_fdr < 0.05, na.rm = TRUE)
extra_agevar_p001 <- sum(extramacular_age_var$age_var_pvalue < 0.01, na.rm = TRUE)

message(sprintf("\nAge-Driven Variance Results:"))
message(sprintf("  Macular: %d (FDR<0.05), %d (p<0.01)", mac_agevar_fdr05, mac_agevar_p001))
message(sprintf("  Extramacular: %d (FDR<0.05), %d (p<0.01)", extra_agevar_fdr05, extra_agevar_p001))

# Save
write.csv(macular_age_var,
          file.path(OUTPUT_DIR, "macular_age_variance.csv"),
          row.names = FALSE)
write.csv(extramacular_age_var,
          file.path(OUTPUT_DIR, "extramacular_age_variance.csv"),
          row.names = FALSE)

#%%
# Update comprehensive comparison
message("\n=== Updating Comprehensive Comparison ===")

# Load existing comparison
comparison_data <- read.csv(file.path(OUTPUT_DIR, "comprehensive_comparison.csv"))

# Add new rows for age-adjusted variance
new_rows <- rbind(
  data.frame(
    Method = "Disease Variance (Age-Adjusted)",
    Tissue = "Macular",
    Features_FDR05 = mac_adj_fdr05,
    Features_P001 = mac_adj_p001,
    Type = "Genes"
  ),
  data.frame(
    Method = "Disease Variance (Age-Adjusted)",
    Tissue = "Extramacular",
    Features_FDR05 = extra_adj_fdr05,
    Features_P001 = extra_adj_p001,
    Type = "Genes"
  ),
  data.frame(
    Method = "Age Variance",
    Tissue = "Macular",
    Features_FDR05 = mac_agevar_fdr05,
    Features_P001 = mac_agevar_p001,
    Type = "Genes"
  ),
  data.frame(
    Method = "Age Variance",
    Tissue = "Extramacular",
    Features_FDR05 = extra_agevar_fdr05,
    Features_P001 = extra_agevar_p001,
    Type = "Genes"
  )
)

comparison_data <- rbind(comparison_data, new_rows)

# Save updated comparison
write.csv(comparison_data, 
          file.path(OUTPUT_DIR, "comprehensive_comparison.csv"),
          row.names = FALSE)

message("\nUpdated comprehensive_comparison.csv")

#%%
# Create updated comparison plot
library(tidyr)

comparison_long <- comparison_data %>%
  pivot_longer(cols = c(Features_FDR05, Features_P001),
               names_to = "Threshold",
               values_to = "Count") %>%
  filter(!is.na(Count))

comparison_long$Threshold <- ifelse(comparison_long$Threshold == "Features_FDR05",
                                   "FDR < 0.05", "p < 0.01")

# Order methods logically
method_order <- c(
  "Original DE (individual stages)",
  "Improved DE: Combined Early",
  "Improved DE: Combined Intermediate", 
  "Improved DE: All Dry AMD",
  "Differential Variability (F-test)",
  "Disease Variance (Age-Adjusted)",
  "Age Variance",
  "Age Main Effect",
  "Age × Disease Interaction"
)

comparison_long$Method <- factor(comparison_long$Method, 
                                levels = method_order)

p <- ggplot(comparison_long, aes(x = Method, y = Count, fill = Threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Tissue, scales = "free_y", ncol = 1) +
  labs(
    title = "Comprehensive Signal Detection Comparison (Updated)",
    x = "Analysis Method",
    y = "Number of Genes Detected",
    fill = "Significance"
  ) +
  theme_classic(base_family = "serif") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("FDR < 0.05" = "#1f78b4", "p < 0.01" = "#e31a1c"))

ggsave(file.path(OUTPUT_DIR, "comprehensive_comparison_plot_updated.pdf"), 
       p, width = 12, height = 10)

message("\nSaved: comprehensive_comparison_plot_updated.pdf")

message("\n=== Analysis Complete ===")
