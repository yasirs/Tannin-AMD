# GSE29801 Dry AMD Differential Expression Analysis with Covariate Evaluation
# Uses limma-voom for microarray data

#%%
# Load libraries
library(limma)
library(edgeR)
library(car)
library(ggplot2)
library(dplyr)
library(pheatmap)

#%%
# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- DATA_DIR

# Visualization settings
FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

# DE thresholds
FDR_THRESH <- 0.05
LOGFC_THRESH <- 0.5

#%%
# Function to run DE analysis with different covariate models
run_de_analysis <- function(expr, metadata, tissue_name, stage, covariates_to_test) {
  
  tissue_slug <- gsub(" ", "_", tolower(tissue_name))
  
  message(sprintf("\n  Analyzing %s vs Control...", stage))
  
  # Subset to control + this stage
  subset_idx <- which(metadata$dry_stage %in% c("Control", stage))
  
  if (length(subset_idx) < 10) {
    message(sprintf("    Skipping: only %d samples", length(subset_idx)))
    return(NULL)
  }
  
  expr_subset <- expr[subset_idx, ]
  meta_subset <- metadata[subset_idx, ]
  
  # Create binary disease variable
  meta_subset$is_disease <- ifelse(meta_subset$dry_stage == stage, 1, 0)
  
  message(sprintf("    Samples: %d control, %d %s", 
                  sum(meta_subset$is_disease == 0), 
                  sum(meta_subset$is_disease == 1), stage))
  
  # Test different covariate models
  results_list <- list()
  de_counts <- data.frame(model = character(), n_de_genes = integer(), stringsAsFactors = FALSE)
  
  # Model 0: No covariates
  message("    Model 0: ~ disease (no covariates)")
  design0 <- model.matrix(~ is_disease, data = meta_subset)
  fit0 <- lmFit(t(expr_subset), design0)
  fit0 <- eBayes(fit0)
  results0 <- topTable(fit0, coef = "is_disease", number = Inf, sort.by = "none")
  results0$gene <- rownames(results0)
  
  n_de0 <- sum(results0$adj.P.Val < FDR_THRESH & abs(results0$logFC) > LOGFC_THRESH, na.rm = TRUE)
  de_counts <- rbind(de_counts, data.frame(model = "No_covariates", n_de_genes = n_de0))
  results_list[["No_covariates"]] <- results0
  
  # Model 1: + age
  if ("age" %in% covariates_to_test && !all(is.na(meta_subset$age))) {
    message("    Model 1: ~ disease + age")
    design1 <- model.matrix(~ is_disease + age, data = meta_subset)
    fit1 <- lmFit(t(expr_subset), design1)
    fit1 <- eBayes(fit1)
    results1 <- topTable(fit1, coef = "is_disease", number = Inf, sort.by = "none")
    results1$gene <- rownames(results1)
    
    n_de1 <- sum(results1$adj.P.Val < FDR_THRESH & abs(results1$logFC) > LOGFC_THRESH, na.rm = TRUE)
    de_counts <- rbind(de_counts, data.frame(model = "Age", n_de_genes = n_de1))
    results_list[["Age"]] <- results1
  }
  
  # Model 2: + sex
  if ("sex" %in% covariates_to_test && length(unique(meta_subset$sex)) > 1) {
    message("    Model 2: ~ disease + sex")
    design2 <- model.matrix(~ is_disease + sex, data = meta_subset)
    fit2 <- lmFit(t(expr_subset), design2)
    fit2 <- eBayes(fit2)
    results2 <- topTable(fit2, coef = "is_disease", number = Inf, sort.by = "none")
    results2$gene <- rownames(results2)
    
    n_de2 <- sum(results2$adj.P.Val < FDR_THRESH & abs(results2$logFC) > LOGFC_THRESH, na.rm = TRUE)
    de_counts <- rbind(de_counts, data.frame(model = "Sex", n_de_genes = n_de2))
    results_list[["Sex"]] <- results2
  }
  
  # Model 3: + age + sex
  if ("age" %in% covariates_to_test && "sex" %in% covariates_to_test && 
      !all(is.na(meta_subset$age)) && length(unique(meta_subset$sex)) > 1) {
    message("    Model 3: ~ disease + age + sex")
    design3 <- model.matrix(~ is_disease + age + sex, data = meta_subset)
    fit3 <- lmFit(t(expr_subset), design3)
    fit3 <- eBayes(fit3)
    results3 <- topTable(fit3, coef = "is_disease", number = Inf, sort.by = "none")
    results3$gene <- rownames(results3)
    
    n_de3 <- sum(results3$adj.P.Val < FDR_THRESH & abs(results3$logFC) > LOGFC_THRESH, na.rm = TRUE)
    de_counts <- rbind(de_counts, data.frame(model = "Age_Sex", n_de_genes = n_de3))
    results_list[["Age_Sex"]] <- results3
  }
  
  # Model 4: + age + sex + RIN (if available)
  if ("rin" %in% covariates_to_test && "age" %in% covariates_to_test && "sex" %in% covariates_to_test) {
    # Check if we have enough non-NA RIN values
    complete_idx <- complete.cases(meta_subset[, c("is_disease", "age", "sex", "rin")])
    
    if (sum(complete_idx) >= 10 && sum(meta_subset$is_disease[complete_idx]) >= 3) {
      message(sprintf("    Model 4: ~ disease + age + sex + RIN (n=%d with complete data)", sum(complete_idx)))
      
      # Subset to complete cases
      expr_complete <- expr_subset[complete_idx, ]
      meta_complete <- meta_subset[complete_idx, ]
      
      design4 <- model.matrix(~ is_disease + age + sex + rin, data = meta_complete)
      fit4 <- lmFit(t(expr_complete), design4)
      fit4 <- eBayes(fit4)
      results4 <- topTable(fit4, coef = "is_disease", number = Inf, sort.by = "none")
      results4$gene <- rownames(results4)
      
      n_de4 <- sum(results4$adj.P.Val < FDR_THRESH & abs(results4$logFC) > LOGFC_THRESH, na.rm = TRUE)
      de_counts <- rbind(de_counts, data.frame(model = "Age_Sex_RIN", n_de_genes = n_de4))
      results_list[["Age_Sex_RIN"]] <- results4
    } else {
      message("    Model 4: Skipped (insufficient complete RIN data)")
    }
  }
  
  # Save results
  for (model_name in names(results_list)) {
    filename <- sprintf("%s_%s_vs_control_%s.csv", tissue_slug, stage, model_name)
    write.csv(results_list[[model_name]], file.path(OUTPUT_DIR, filename), row.names = FALSE)
  }
  
  return(list(
    de_counts = de_counts,
    results = results_list,
    stage = stage,
    tissue = tissue_name
  ))
}

#%%
# Function to check VIF for covariates
check_vif <- function(metadata) {
  message("\nChecking multicollinearity (VIF)...")
  
  # Create design matrix with all covariates
  complete_cases <- complete.cases(metadata[, c("age", "sex")])
  meta_complete <- metadata[complete_cases, ]
  
  if (nrow(meta_complete) < 10) {
    message("  Not enough complete cases for VIF analysis")
    return(NULL)
  }
  
  # Fit model with all covariates
  model <- lm(age ~ sex, data = meta_complete)
  
  tryCatch({
    vif_values <- vif(model)
    message("  VIF values:")
    print(vif_values)
    
    if (any(vif_values > 5)) {
      message("  ⚠️  High multicollinearity detected (VIF > 5)")
    } else {
      message("  ✅ No problematic multicollinearity (all VIF < 5)")
    }
    
    return(vif_values)
  }, error = function(e) {
    message("  Could not calculate VIF (model may be singular)")
    return(NULL)
  })
}

#%%
# Main analysis
message("=== Dry AMD Differential Expression Analysis ===")

# Load data
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

# Check VIF for each tissue
message("\n--- Macular RPE-choroid ---")
check_vif(macular_data$metadata)

message("\n--- Extramacular RPE-choroid ---")
check_vif(extramacular_data$metadata)

#%%
# Run DE analyses
covariates <- c("age", "sex", "rin")
stages <- c("MD1", "MD2", "dry_AMD", "GA")

all_results <- list()

# Macular
message("\n\n=== MACULAR RPE-CHOROID ===")
macular_results <- list()
for (stage in stages) {
  result <- run_de_analysis(
    macular_data$expr,
    macular_data$metadata,
    "Macular RPE-choroid",
    stage,
    covariates
  )
  if (!is.null(result)) {
    macular_results[[stage]] <- result
  }
}

# Extramacular
message("\n\n=== EXTRAMACULAR RPE-CHOROID ===")
extramacular_results <- list()
for (stage in stages) {
  result <- run_de_analysis(
    extramacular_data$expr,
    extramacular_data$metadata,
    "Extramacular RPE-choroid",
    stage,
    covariates
  )
  if (!is.null(result)) {
    extramacular_results[[stage]] <- result
  }
}

#%%
# Combine DE counts for visualization
message("\n\n=== Creating Summary Visualizations ===")

all_de_counts <- data.frame()

for (stage in names(macular_results)) {
  df <- macular_results[[stage]]$de_counts
  df$stage <- stage
  df$tissue <- "Macular"
  all_de_counts <- rbind(all_de_counts, df)
}

for (stage in names(extramacular_results)) {
  df <- extramacular_results[[stage]]$de_counts
  df$stage <- stage
  df$tissue <- "Extramacular"
  all_de_counts <- rbind(all_de_counts, df)
}

# Save DE counts summary
write.csv(all_de_counts, file.path(OUTPUT_DIR, "de_gene_counts_summary.csv"), row.names = FALSE)

#%%
# Plot 1: Bar plot of DE genes per model per comparison
p1 <- ggplot(all_de_counts, aes(x = stage, y = n_de_genes, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ tissue) +
  labs(
    title = "DE Genes Detected by Model and Stage",
    x = "AMD Stage",
    y = "Number of DE Genes",
    fill = "Model"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  scale_fill_brewer(palette = "Set2")

ggsave(file.path(OUTPUT_DIR, "de_counts_by_model.pdf"), p1, width = 10, height = 6)

#%%
# Plot 2: Heatmap of DE gene counts
de_matrix <- all_de_counts %>%
  tidyr::pivot_wider(names_from = model, values_from = n_de_genes, values_fill = 0) %>%
  as.data.frame()

rownames(de_matrix) <- paste(de_matrix$tissue, de_matrix$stage, sep = "_")
de_matrix <- de_matrix[, !(colnames(de_matrix) %in% c("tissue", "stage"))]

pdf(file.path(OUTPUT_DIR, "de_counts_heatmap.pdf"), width = 8, height = 6)
pheatmap(
  as.matrix(de_matrix),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.0f",
  main = "DE Gene Counts: Rows=Comparisons, Cols=Models",
  color = colorRampPalette(c("white", "orange", "red"))(50)
)
dev.off()

#%%
# Summary
message("\n=== Analysis Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
message("\nGenerated files:")
message("  - Individual DE results: [tissue]_[stage]_vs_control_[model].csv")
message("  - de_gene_counts_summary.csv")
message("  - de_counts_by_model.pdf")
message("  - de_counts_heatmap.pdf")

message("\nSummary of DE genes detected:")
print(all_de_counts)
