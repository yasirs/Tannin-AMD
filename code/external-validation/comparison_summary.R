# Extended Dry AMD - Comprehensive Comparison Summary
# Aggregates results from all analyses and creates comparison visualizations

#%%
library(dplyr)
library(ggplot2)
library(tidyr)

#%%
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- DATA_DIR

FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

#%%
# Load all count summaries
message("=== Loading Analysis Results ===\n")

# Original DE counts
original_counts <- read.csv(file.path(OUTPUT_DIR, "de_gene_counts_summary.csv"))

# Improved DE counts
if (file.exists(file.path(OUTPUT_DIR, "improved_de_counts_summary.csv"))) {
  improved_counts <- read.csv(file.path(OUTPUT_DIR, "improved_de_counts_summary.csv"))
} else {
  improved_counts <- NULL
}

# Variance analysis
if (file.exists(file.path(OUTPUT_DIR, "variance_analysis_summary.csv"))) {
  variance_counts <- read.csv(file.path(OUTPUT_DIR, "variance_analysis_summary.csv"))
} else {
  variance_counts <- NULL
}

# Age analysis
if (file.exists(file.path(OUTPUT_DIR, "age_analysis_summary.csv"))) {
  age_counts <- read.csv(file.path(OUTPUT_DIR, "age_analysis_summary.csv"))
} else {
  age_counts <- NULL
}

# GSVA pathway analysis
if (file.exists(file.path(OUTPUT_DIR, "gsva_pathway_counts.csv"))) {
  gsva_counts <- read.csv(file.path(OUTPUT_DIR, "gsva_pathway_counts.csv"))
} else {
  gsva_counts <- NULL
}

#%%
# Create comprehensive comparison table
message("Creating comprehensive comparison table...\n")

comparison_data <- data.frame(
  Method = character(),
  Tissue = character(),
  Features_FDR05 = integer(),
  Features_P001 = integer(),
  Type = character(),
  stringsAsFactors = FALSE
)

# Add original DE results (just FDR < 0.05 equivalent, which was 0-3 genes)
if (nrow(original_counts) > 0) {
  original_summary <- original_counts %>%
    group_by(tissue, stage) %>%
    summarise(max_genes = max(n_de_genes, na.rm = TRUE), .groups = "drop")
  
  for (i in 1:nrow(original_summary)) {
    comparison_data <- rbind(comparison_data, data.frame(
      Method = "Original DE (individual stages)",
      Tissue = original_summary$tissue[i],
      Features_FDR05 = original_summary$max_genes[i],
      Features_P001 = NA,
      Type = "Genes"
    ))
  }
}

# Add improved DE results
if (!is.null(improved_counts)) {
  for (i in 1:nrow(improved_counts)) {
    comparison_data <- rbind(comparison_data, data.frame(
      Method = paste("Improved DE:", improved_counts$comparison[i]),
      Tissue = improved_counts$tissue[i],
      Features_FDR05 = improved_counts$fdr_0.05[i],
      Features_P001 = improved_counts$p_0.01[i],
      Type = "Genes"
    ))
  }
}

# Add variance results
if (!is.null(variance_counts)) {
  for (i in 1:nrow(variance_counts)) {
    comparison_data <- rbind(comparison_data, data.frame(
      Method = "Differential Variability (F-test)",
      Tissue = variance_counts$tissue[i],
      Features_FDR05 = variance_counts$f_test_fdr05[i],
      Features_P001 = variance_counts$f_test_p001[i],
      Type = "Genes"
    ))
  }
}

# Add age results
if (!is.null(age_counts)) {
  for (i in 1:nrow(age_counts)) {
    comparison_data <- rbind(comparison_data, data.frame(
      Method = "Age Main Effect",
      Tissue = age_counts$tissue[i],
      Features_FDR05 = age_counts$age_main_fdr05[i],
      Features_P001 = age_counts$age_main_p001[i],
      Type = "Genes"
    ))
    
    comparison_data <- rbind(comparison_data, data.frame(
      Method = "Age × Disease Interaction",
      Tissue = age_counts$tissue[i],
      Features_FDR05 = age_counts$interaction_fdr05[i],
      Features_P001 = age_counts$interaction_p001[i],
      Type = "Genes"
    ))
  }
}

# Add GSVA results
if (!is.null(gsva_counts)) {
  for (i in 1:nrow(gsva_counts)) {
    comparison_data <- rbind(comparison_data, data.frame(
      Method = paste("GSVA:", gsva_counts$gene_set[i]),
      Tissue = gsva_counts$tissue[i],
      Features_FDR05 = gsva_counts$pathways_fdr05[i],
      Features_P001 = gsva_counts$pathways_p001[i],
      Type = "Pathways"
    ))
  }
}

# Save comparison table
write.csv(comparison_data, file.path(OUTPUT_DIR, "comprehensive_comparison.csv"), row.names = FALSE)

message("Comprehensive Comparison Table:")
print(comparison_data)

#%%
# Visualization: Bar plot comparing all methods
comparison_long <- comparison_data %>%
  pivot_longer(cols = c(Features_FDR05, Features_P001),
               names_to = "Threshold",
               values_to = "Count") %>%
  filter(!is.na(Count))

comparison_long$Threshold <- ifelse(comparison_long$Threshold == "Features_FDR05",
                                   "FDR < 0.05", "p < 0.01")

p <- ggplot(comparison_long, aes(x = Method, y = Count, fill = Threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Type ~ Tissue, scales = "free") +
  labs(
    title = "Comprehensive Signal Detection Comparison",
    x = "Analysis Method",
    y = "Number of Features Detected",
    fill = "Significance"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(face = "bold")
  ) +
  scale_fill_manual(values = c("FDR < 0.05" = "#1f78b4", "p < 0.01" = "#e31a1c"))

ggsave(file.path(OUTPUT_DIR, "comprehensive_comparison_plot.pdf"), p, width = 14, height = 10)

message("\n=== Comprehensive Comparison Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
message("\nKey findings:")
message("  - Check which methods detected the most features")
message("  - Compare gene-level vs pathway-level signal")
message("  - Assess value of age analysis and variance testing")
