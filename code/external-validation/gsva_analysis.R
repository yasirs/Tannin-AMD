# Extended Dry AMD - GSVA Pathway Analysis
# Tests pathway-level differences using GSVA

#%%
library(GSVA)
library(msigdbr)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)

#%%
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
message("Loading data...")
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

#%%
# Load gene sets from MSigDB
message("\nLoading gene sets from MSigDB...")

# H: Hallmark
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
message(sprintf("  Hallmark: %d gene sets", length(hallmark_list)))

# C5: GO BP (limit to avoid memory issues)
gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gobp_list <- split(gobp$gene_symbol, gobp$gs_name)
# Filter to pathways with 15-500 genes
gobp_sizes <- sapply(gobp_list, length)
gobp_list <- gobp_list[gobp_sizes >= 15 & gobp_sizes <= 500]
message(sprintf("  GO BP (filtered 15-500 genes): %d gene sets", length(gobp_list)))

#%%
# Function to run GSVA and test differences
run_gsva_analysis <- function(expr, metadata, gene_sets, set_name, tissue_name) {
  
  message(sprintf("\nRunning GSVA for %s - %s...", tissue_name, set_name))
  
  # Convert to matrix if needed
  expr_mat <- as.matrix(t(expr))
  
  # Run GSVA
  message("  Computing GSVA scores...")
  gsva_scores <- gsva(expr_mat, gene_sets, method = "gsva", verbose = FALSE)
  
  message(sprintf("  GSVA scores: %d pathways × %d samples", 
                  nrow(gsva_scores), ncol(gsva_scores)))
  
  # Test pathway differences
  message("  Testing pathway differences...")
  
  # Create design matrix
  metadata$is_disease <- ifelse(metadata$disease_status == "AMD", 1, 0)
  design <- model.matrix(~ is_disease + age + sex, data = metadata)
  
  # Fit model
  fit <- lmFit(gsva_scores, design)
  fit <- eBayes(fit)
  results <- topTable(fit, coef = "is_disease", number = Inf, sort.by = "P")
  results$pathway <- rownames(results)
  
  # Count significant pathways
  n_fdr05 <- sum(results$adj.P.Val < FDR_THRESH, na.rm = TRUE)
  n_p001 <- sum(results$P.Value < P_THRESH, na.rm = TRUE)
  
  message(sprintf("  Significant pathways: %d (FDR < 0.05), %d (p < 0.01)", n_fdr05, n_p001))
  
  return(list(
    gsva_scores = gsva_scores,
    results = results,
    n_fdr05 = n_fdr05,
    n_p001 = n_p001
  ))
}

#%%
# Run GSVA for all gene set collections
gsva_results <- list()

for (gset_name in c("Hallmark", "GOBP")) {
  
  gene_sets <- switch(gset_name,
    "Hallmark" = hallmark_list,
    "GOBP" = gobp_list
  )
  
  # Macular
  gsva_results[[paste0("Macular_", gset_name)]] <- run_gsva_analysis(
    macular_data$expr,
    macular_data$metadata,
    gene_sets,
    gset_name,
    "Macular"
  )
  
  # Extramacular
  gsva_results[[paste0("Extramacular_", gset_name)]] <- run_gsva_analysis(
    extramacular_data$expr,
    extramacular_data$metadata,
    gene_sets,
    gset_name,
    "Extramacular"
  )
}

#%%
# Save results
message("\nSaving GSVA results...")

for (name in names(gsva_results)) {
  # Save pathway results
  write.csv(gsva_results[[name]]$results,
            file.path(OUTPUT_DIR, sprintf("gsva_%s_results.csv", tolower(name))),
            row.names = FALSE)
  
  # Save GSVA scores
  write.csv(gsva_results[[name]]$gsva_scores,
            file.path(OUTPUT_DIR, sprintf("gsva_%s_scores.csv", tolower(name))))
}

#%%
# Summary of pathway counts
pathway_counts <- data.frame(
  analysis = names(gsva_results),
  pathways_fdr05 = sapply(gsva_results, function(x) x$n_fdr05),
  pathways_p001 = sapply(gsva_results, function(x) x$n_p001)
)

pathway_counts$tissue <- ifelse(grepl("Macular", pathway_counts$analysis), "Macular", "Extramacular")
pathway_counts$gene_set <- gsub(".*_", "", pathway_counts$analysis)

write.csv(pathway_counts, file.path(OUTPUT_DIR, "gsva_pathway_counts.csv"), row.names = FALSE)

message("\n=== GSVA Pathway Counts ===")
print(pathway_counts)

#%%
# Visualization: Top pathways heatmap (Hallmark)
plot_top_pathways <- function(gsva_result, metadata, title, top_n = 20) {
  
  # Get top pathways
  top_pathways <- head(gsva_result$results, top_n)$pathway
  
  if (length(top_pathways) == 0) {
    message(sprintf("  No significant pathways for %s", title))
    return(NULL)
  }
  
  # Subset GSVA scores
  scores_subset <- gsva_result$gsva_scores[top_pathways, , drop = FALSE]
  
  # Annotation
  annotation <- data.frame(
    Disease = metadata$disease_status,
    row.names = colnames(scores_subset)
  )
  
  annotation_colors <- list(
    Disease = c("Control" = "navy", "AMD" = "firebrick")
  )
  
  # Clean pathway names
  rownames(scores_subset) <- gsub("HALLMARK_|KEGG_|GOBP_", "", rownames(scores_subset))
  rownames(scores_subset) <- gsub("_", " ", rownames(scores_subset))
  
  # Plot
  pheatmap(
    scores_subset,
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    show_colnames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "row",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = title,
    fontsize_row = 8
  )
}

# Plot Hallmark for macular
pdf(file.path(OUTPUT_DIR, "gsva_macular_hallmark_top20.pdf"), width = 10, height = 8)
plot_top_pathways(gsva_results$Macular_Hallmark, macular_data$metadata, 
                  "Macular: Top 20 Hallmark Pathways")
dev.off()

message("\n=== GSVA Analysis Complete ===")
message(sprintf("Output directory: %s", OUTPUT_DIR))
