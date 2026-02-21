# Extended Dry AMD - Improved DE Detection Analysis
# Combines stages, relaxes thresholds, performs meta-analysis

#%%
# Load libraries
library(limma)
library(edgeR)
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

# DE thresholds
FDR_STRICT <- 0.05
FDR_RELAXED <- 0.10
P_NOMINAL <- 0.01
LOGFC_THRESH <- 0.5

#%%
# Load data
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

#%%
# Function to run DE with different thresholds
test_thresholds <- function(expr, metadata, comparison_name, tissue_name) {
  
  # Create binary disease variable
  is_disease <- ifelse(metadata$disease_status == "AMD", 1, 0)
  
  # Design matrix with age and sex
  design <- model.matrix(~ is_disease + age + sex, data = metadata)
  
  # Fit model
  fit <- lmFit(t(expr), design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "is_disease", number = Inf, sort.by = "none")
  results$gene <- rownames(results)
  
  # Count DE genes at different thresholds
  counts <- data.frame(
    comparison = comparison_name,
    tissue = tissue_name,
    fdr_0.05 = sum(results$adj.P.Val < FDR_STRICT & abs(results$logFC) > LOGFC_THRESH, na.rm = TRUE),
    fdr_0.10 = sum(results$adj.P.Val < FDR_RELAXED & abs(results$logFC) > LOGFC_THRESH, na.rm = TRUE),
    p_0.01 = sum(results$P.Value < P_NOMINAL & abs(results$logFC) > LOGFC_THRESH, na.rm = TRUE)
  )
  
  return(list(results = results, counts = counts))
}

#%%
# 1. COMBINED STAGES ANALYSIS
message("=== 1. Combined Stages Analysis ===\n")

# Create combined stage labels
combine_stages <- function(data_list) {
  metadata <- data_list$metadata
  metadata$combined_stage <- case_when(
    metadata$dry_stage == "Control" ~ "Control",
    metadata$dry_stage %in% c("MD1", "MD2") ~ "Early",
    metadata$dry_stage == "dry_AMD" ~ "Intermediate",
    metadata$dry_stage == "GA" ~ "Advanced"
  )
  metadata$combined_stage <- factor(metadata$combined_stage, 
                                    levels = c("Control", "Early", "Intermediate", "Advanced"))
  return(list(expr = data_list$expr, metadata = metadata))
}

macular_combined <- combine_stages(macular_data)
extramacular_combined <- combine_stages(extramacular_data)

# Test each combined stage
all_counts <- data.frame()

for (stage in c("Early", "Intermediate", "Advanced")) {
  # Macular
  idx_mac <- which(macular_combined$metadata$combined_stage %in% c("Control", stage))
  if (sum(macular_combined$metadata$combined_stage[idx_mac] == stage) >= 5) {
    message(sprintf("Testing Macular: %s vs Control...", stage))
    result_mac <- test_thresholds(
      macular_combined$expr[idx_mac, ],
      macular_combined$metadata[idx_mac, ],
      paste("Combined", stage),
      "Macular"
    )
    all_counts <- rbind(all_counts, result_mac$counts)
    
    # Save results
    filename <- sprintf("macular_combined_%s_vs_control.csv", tolower(stage))
    write.csv(result_mac$results, file.path(OUTPUT_DIR, filename), row.names = FALSE)
  }
  
  # Extramacular
  idx_extra <- which(extramacular_combined$metadata$combined_stage %in% c("Control", stage))
  if (sum(extramacular_combined$metadata$combined_stage[idx_extra] == stage) >= 5) {
    message(sprintf("Testing Extramacular: %s vs Control...", stage))
    result_extra <- test_thresholds(
      extramacular_combined$expr[idx_extra, ],
      extramacular_combined$metadata[idx_extra, ],
      paste("Combined", stage),
      "Extramacular"
    )
    all_counts <- rbind(all_counts, result_extra$counts)
    
    # Save results
    filename <- sprintf("extramacular_combined_%s_vs_control.csv", tolower(stage))
    write.csv(result_extra$results, file.path(OUTPUT_DIR, filename), row.names = FALSE)
  }
}

#%%
# 2. ALL DRY AMD POOLED VS CONTROL
message("\n=== 2. All Dry AMD Pooled Analysis ===\n")

# Macular: All AMD vs Control
message("Testing Macular: All Dry AMD vs Control...")
result_mac_all <- test_thresholds(
  macular_data$expr,
  macular_data$metadata,
  "All Dry AMD",
  "Macular"
)
all_counts <- rbind(all_counts, result_mac_all$counts)
write.csv(result_mac_all$results, file.path(OUTPUT_DIR, "macular_all_dry_amd_vs_control.csv"), row.names = FALSE)

# Extramacular: All AMD vs Control
message("Testing Extramacular: All Dry AMD vs Control...")
result_extra_all <- test_thresholds(
  extramacular_data$expr,
  extramacular_data$metadata,
  "All Dry AMD",
  "Extramacular"
)
all_counts <- rbind(all_counts, result_extra_all$counts)
write.csv(result_extra_all$results, file.path(OUTPUT_DIR, "extramacular_all_dry_amd_vs_control.csv"), row.names = FALSE)

#%%
# 3. META-ANALYSIS (Fisher's method to combine p-values)
message("\n=== 3. Meta-Analysis (Pooling Tissues) ===\n")

# Function to combine p-values using Fisher's method
fisher_combine <- function(p_values) {
  # Remove NA values
  p_values <- p_values[!is.na(p_values)]
  
  # Fisher's combined test statistic
  chi_sq <- -2 * sum(log(p_values))
  
  # Combined p-value (chi-squared distribution with 2k degrees of freedom)
  combined_p <- pchisq(chi_sq, df = 2 * length(p_values), lower.tail = FALSE)
  
  return(combined_p)
}

# For each comparison, combine macular and extramacular
meta_comparisons <- list(
  "Combined Early" = list(
    mac = "macular_combined_early_vs_control.csv",
    extra = "extramacular_combined_early_vs_control.csv"
  ),
  "Combined Intermediate" = list(
    mac = "macular_combined_intermediate_vs_control.csv",
    extra = "extramacular_combined_intermediate_vs_control.csv"
  ),
  "All Dry AMD" = list(
    mac = "macular_all_dry_amd_vs_control.csv",
    extra = "extramacular_all_dry_amd_vs_control.csv"
  )
)

meta_results <- list()

for (comp_name in names(meta_comparisons)) {
  files <- meta_comparisons[[comp_name]]
  
  # Check if both files exist
  mac_file <- file.path(OUTPUT_DIR, files$mac)
  extra_file <- file.path(OUTPUT_DIR, files$extra)
  
  if (file.exists(mac_file) && file.exists(extra_file)) {
    message(sprintf("Meta-analyzing: %s...", comp_name))
    
    mac_res <- read.csv(mac_file)
    extra_res <- read.csv(extra_file)
    
    # Merge by gene
    merged <- merge(mac_res[, c("gene", "P.Value", "logFC", "adj.P.Val")],
                    extra_res[, c("gene", "P.Value", "logFC", "adj.P.Val")],
                    by = "gene", suffixes = c("_mac", "_extra"))
    
    # Combine p-values
    merged$meta_p <- mapply(function(p1, p2) fisher_combine(c(p1, p2)),
                            merged$P.Value_mac, merged$P.Value_extra)
    
    # Meta logFC (average weighted by inverse variance - simplified as mean)
    merged$meta_logFC <- (merged$logFC_mac + merged$logFC_extra) / 2
    
    # FDR correction
    merged$meta_fdr <- p.adjust(merged$meta_p, method = "BH")
    
    # Count DE genes
    meta_count <- data.frame(
      comparison = comp_name,
      tissue = "Meta-analysis",
      fdr_0.05 = sum(merged$meta_fdr < FDR_STRICT & abs(merged$meta_logFC) > LOGFC_THRESH, na.rm = TRUE),
      fdr_0.10 = sum(merged$meta_fdr < FDR_RELAXED & abs(merged$meta_logFC) > LOGFC_THRESH, na.rm = TRUE),
      p_0.01 = sum(merged$meta_p < P_NOMINAL & abs(merged$meta_logFC) > LOGFC_THRESH, na.rm = TRUE)
    )
    all_counts <- rbind(all_counts, meta_count)
    
    # Save
    meta_filename <- sprintf("meta_%s.csv", gsub(" ", "_", tolower(comp_name)))
    write.csv(merged, file.path(OUTPUT_DIR, meta_filename), row.names = FALSE)
    
    meta_results[[comp_name]] <- merged
  }
}

#%%
# 4. CROSS-SUBTYPE SIMILARITY ANALYSIS
message("\n=== 4. Cross-Subtype Similarity ===\n")

# Load individual stage results from original analysis
load_original_results <- function(tissue, stage) {
  filename <- sprintf("%s_rpe-choroid_%s_vs_control_Age_Sex.csv", 
                      tolower(tissue), stage)
  path <- file.path(OUTPUT_DIR, filename)
  if (file.exists(path)) {
    return(read.csv(path))
  }
  return(NULL)
}

# Calculate pairwise similarities
stages <- c("MD1", "MD2", "dry_AMD")
tissues <- c("macular", "extramacular")

similarity_results <- list()

for (tissue in tissues) {
  # Load all stage results
  stage_results <- list()
  for (stage in stages) {
    res <- load_original_results(tissue, stage)
    if (!is.null(res)) {
      stage_results[[stage]] <- res
    }
  }
  
  if (length(stage_results) >= 2) {
    # Correlation matrix of logFC
    logfc_matrix <- do.call(cbind, lapply(stage_results, function(x) {
      setNames(x$logFC, x$gene)
    }))
    
    cor_matrix <- cor(logfc_matrix, use = "pairwise.complete.obs", method = "spearman")
    
    message(sprintf("\n%s tissue - logFC correlations:", tools::toTitleCase(tissue)))
    print(round(cor_matrix, 3))
    
    similarity_results[[tissue]] <- list(
      correlation = cor_matrix,
      stage_results = stage_results
    )
  }
}

#%%
# Save all counts summary
write.csv(all_counts, file.path(OUTPUT_DIR, "improved_de_counts_summary.csv"), row.names = FALSE)

# Print summary
message("\n=== Summary: DE Genes Detected ===")
print(all_counts)

#%%
# Visualization: Comparison of detection methods
all_counts_long <- all_counts %>%
  tidyr::pivot_longer(cols = c(fdr_0.05, fdr_0.10, p_0.01),
                     names_to = "threshold", values_to = "n_genes")

p <- ggplot(all_counts_long, aes(x = comparison, y = n_genes, fill = threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ tissue, scales = "free_x") +
  labs(
    title = "DE Genes Detected: Improved Methods",
    x = "Comparison",
    y = "Number of DE Genes",
    fill = "Threshold"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(
    values = c("fdr_0.05" = "#1f78b4", "fdr_0.10" = "#33a02c", "p_0.01" = "#e31a1c"),
    labels = c("FDR < 0.05", "FDR < 0.10", "p < 0.01")
  )

ggsave(file.path(OUTPUT_DIR, "improved_de_comparison.pdf"), p, width = 10, height = 6)

message("\n=== Improved DE Analysis Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
