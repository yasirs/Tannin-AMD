# Extended Dry AMD - Age Effects Analysis
# Tests age main effect and age × disease interaction

#%%
library(limma)
library(dplyr)
library(ggplot2)
library(VennDiagram)

#%%
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- DATA_DIR

FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

FDR_THRESH <- 0.05
P_THRESH <- 0.01
LOGFC_THRESH <- 0.5

#%%
# Load data
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

#%%
# Function to test age effects
test_age_effects <- function(expr, metadata, tissue_name) {
  
  message(sprintf("\nAnalyzing %s...", tissue_name))
  
  # Create binary disease variable
  metadata$is_disease <- ifelse(metadata$disease_status == "AMD", 1, 0)
  
  # Model 1: Age main effect (ignore disease)
  message("  Model 1: expr ~ age (main effect)")
  design_age <- model.matrix(~ age, data = metadata)
  fit_age <- lmFit(t(expr), design_age)
  fit_age <- eBayes(fit_age)
  results_age <- topTable(fit_age, coef = "age", number = Inf, sort.by = "none")
  results_age$gene <- rownames(results_age)
  
  # Model 2: Age + Disease (additive)
  message("  Model 2: expr ~ age + disease")
  design_add <- model.matrix(~ age + is_disease, data = metadata)
  fit_add <- lmFit(t(expr), design_add)
  fit_add <- eBayes(fit_add)
  results_age_add <- topTable(fit_add, coef = "age", number = Inf, sort.by = "none")
  results_age_add$gene <- rownames(results_age_add)
  results_disease_add <- topTable(fit_add, coef = "is_disease", number = Inf, sort.by = "none")
  results_disease_add$gene <- rownames(results_disease_add)
  
  # Model 3: Age × Disease interaction
  message("  Model 3: expr ~ age * disease")
  design_int <- model.matrix(~ age * is_disease, data = metadata)
  fit_int <- lmFit(t(expr), design_int)
  fit_int <- eBayes(fit_int)
  results_age_int <- topTable(fit_int, coef = "age", number = Inf, sort.by = "none")
  results_age_int$gene <- rownames(results_age_int)
  results_disease_int <- topTable(fit_int, coef = "is_disease", number = Inf, sort.by = "none")
  results_disease_int$gene <- rownames(results_disease_int)
  results_interaction <- topTable(fit_int, coef = "age:is_disease", number = Inf, sort.by = "none")
  results_interaction$gene <- rownames(results_interaction)
  
  return(list(
    age_main = results_age,
    age_additive = results_age_add,
    disease_additive = results_disease_add,
    age_interaction = results_age_int,
    disease_interaction = results_disease_int,
    interaction = results_interaction
  ))
}

#%%
# Run age analysis
macular_age <- test_age_effects(macular_data$expr, macular_data$metadata, "Macular")
extramacular_age <- test_age_effects(extramacular_data$expr, extramacular_data$metadata, "Extramacular")

# Save results
for (name in names(macular_age)) {
  write.csv(macular_age[[name]], 
            file.path(OUTPUT_DIR, sprintf("macular_age_%s.csv", name)), 
            row.names = FALSE)
  write.csv(extramacular_age[[name]], 
            file.path(OUTPUT_DIR, sprintf("extramacular_age_%s.csv", name)), 
            row.names = FALSE)
}

#%%
# Count significant genes
count_age_genes <- function(age_results, tissue) {
  data.frame(
    tissue = tissue,
    age_main_fdr05 = sum(age_results$age_main$adj.P.Val < FDR_THRESH, na.rm = TRUE),
    age_main_p001 = sum(age_results$age_main$P.Value < P_THRESH, na.rm = TRUE),
    disease_add_fdr05 = sum(age_results$disease_additive$adj.P.Val < FDR_THRESH & 
                           abs(age_results$disease_additive$logFC) > LOGFC_THRESH, na.rm = TRUE),
    disease_add_p001 = sum(age_results$disease_additive$P.Value < P_THRESH & 
                          abs(age_results$disease_additive$logFC) > LOGFC_THRESH, na.rm = TRUE),
    interaction_fdr05 = sum(age_results$interaction$adj.P.Val < FDR_THRESH, na.rm = TRUE),
    interaction_p001 = sum(age_results$interaction$P.Value < P_THRESH, na.rm = TRUE)
  )
}

age_counts <- rbind(
  count_age_genes(macular_age, "Macular"),
  count_age_genes(extramacular_age, "Extramacular")
)

write.csv(age_counts, file.path(OUTPUT_DIR, "age_analysis_summary.csv"), row.names = FALSE)

message("\n=== Age Analysis Summary ===")
print(age_counts)

#%%
# Overlap analysis for macular tissue
message("\n=== Overlap Analysis (Macular) ===")

# Define gene sets (p < 0.01)
age_genes <- macular_age$age_main$gene[macular_age$age_main$P.Value < P_THRESH]
disease_genes <- macular_age$disease_additive$gene[
  macular_age$disease_additive$P.Value < P_THRESH & 
  abs(macular_age$disease_additive$logFC) > LOGFC_THRESH
]
interaction_genes <- macular_age$interaction$gene[macular_age$interaction$P.Value < P_THRESH]

message(sprintf("  Age-associated genes: %d", length(age_genes)))
message(sprintf("  Disease-associated genes: %d", length(disease_genes)))
message(sprintf("  Interaction genes: %d", length(interaction_genes)))

# Overlap counts
overlap_age_disease <- length(intersect(age_genes, disease_genes))
overlap_age_int <- length(intersect(age_genes, interaction_genes))
overlap_disease_int <- length(intersect(disease_genes, interaction_genes))
overlap_all <- length(Reduce(intersect, list(age_genes, disease_genes, interaction_genes)))

message(sprintf("\n  Age ∩ Disease: %d", overlap_age_disease))
message(sprintf("  Age ∩ Interaction: %d", overlap_age_int))
message(sprintf("  Disease ∩ Interaction: %d", overlap_disease_int))
message(sprintf("  All three: %d", overlap_all))

#%%
# Venn diagram
if (length(age_genes) > 0 || length(disease_genes) > 0 || length(interaction_genes) > 0) {
  
  pdf(file.path(OUTPUT_DIR, "macular_age_disease_venn.pdf"), width = 7, height = 7)
  
  grid.newpage()
  venn.plot <- venn.diagram(
    x = list(
      Age = age_genes,
      Disease = disease_genes,
      Interaction = interaction_genes
    ),
    filename = NULL,
    col = "black",
    fill = c("skyblue", "pink", "lightgreen"),
    alpha = 0.5,
    cex = 1.2,
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontfamily = "sans",
    main = "Macular: Age vs Disease vs Interaction Genes (p < 0.01)"
  )
  grid.draw(venn.plot)
  
  dev.off()
}

message("\n=== Age Analysis Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
