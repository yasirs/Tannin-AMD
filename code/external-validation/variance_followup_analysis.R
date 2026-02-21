# High-Variance Gene Follow-Up Analysis
# Analyzes differentially variable genes from GSE29801

#%%
# Load libraries
library(dplyr)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(VennDiagram)
library(mclust)
library(gridExtra)

#%%
# Configuration
BASE_DIR <- "/home/ysuhail/work/Tannin-AMD"
DATA_DIR <- file.path(BASE_DIR, "results/cohort-GSE29801/dry_amd_de")
OUTPUT_DIR <- DATA_DIR

FONT_FAMILY <- "serif"
theme_set(theme_classic(base_family = FONT_FAMILY) +
  theme(axis.text = element_text(color = "black")))

#%%
# Load variance analysis results (annotated)
message("=== Loading Variance Analysis Results ===\n")

macular_var <- read.csv(file.path(DATA_DIR, "macular_variance_analysis_annotated.csv"))
extramacular_var <- read.csv(file.path(DATA_DIR, "extramacular_variance_analysis_annotated.csv"))

message(sprintf("Macular: %d genes", nrow(macular_var)))
message(sprintf("Extramacular: %d genes", nrow(extramacular_var)))

# Load expression data
macular_data <- readRDS(file.path(DATA_DIR, "macular_rpe_choroid_data.rds"))
extramacular_data <- readRDS(file.path(DATA_DIR, "extramacular_rpe_choroid_data.rds"))

#%%
# ========================================
# ANALYSIS 1: PATHWAY ENRICHMENT
# ========================================

message("\n=== ANALYSIS 1: Pathway Enrichment ===\n")

# ========================================
# ANALYSIS 1: PATHWAY ENRICHMENT (GSEA)
# ========================================

message("\n=== ANALYSIS 1: Pathway Enrichment (GSEA) ===\n")

perform_gsea_variance <- function(var_results, tissue_name) {
  
  message(sprintf("\nTissue: %s", tissue_name))
  
  # Create ranked gene list by variance fold change (or F-statistic)
  # Use -log10(p-value) * sign(var_fc) as ranking metric
  gene_ranks <- var_results %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    mutate(rank_score = -log10(f_pvalue + 1e-300) * sign(var_fc)) %>%
    arrange(desc(rank_score)) %>%
    select(gene_symbol, rank_score)
  
  # Remove duplicates (keep highest ranked)
  gene_ranks <- gene_ranks[!duplicated(gene_ranks$gene_symbol), ]
  
  # Create named vector for fgsea
  gene_list <- setNames(gene_ranks$rank_score, gene_ranks$gene_symbol)
  
  message(sprintf("  Created ranked gene list: %d genes", length(gene_list)))
  message(sprintf("  Top gene: %s (score: %.2f)", names(gene_list)[1], gene_list[1]))
  message(sprintf("  Bottom gene: %s (score: %.2f)", 
                  names(gene_list)[length(gene_list)], 
                  gene_list[length(gene_list)]))
  
  # Load gene sets
  message("\n  Loading gene sets...")
  
  # GO BP
  go_bp_sets <- msigdbr(species = "Homo sapiens", 
                       collection = "C5", 
                       subcollection = "GO:BP")
  go_bp_list <- split(go_bp_sets$gene_symbol, go_bp_sets$gs_name)
  
  # Filter to pathways with 15-500 genes
  go_bp_sizes <- sapply(go_bp_list, length)
  go_bp_list <- go_bp_list[go_bp_sizes >= 15 & go_bp_sizes <= 500]
  
  message(sprintf("  GO BP pathways: %d", length(go_bp_list)))
  
  # Hallmark
  hallmark_sets <- msigdbr(species = "Homo sapiens", collection = "H")
  hallmark_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)
  
  message(sprintf("  Hallmark pathways: %d", length(hallmark_list)))
  
  # Run GSEA
  message("\n  Running GSEA on GO BP...")
  gsea_go <- fgsea(pathways = go_bp_list,
                  stats = gene_list,
                  minSize = 15,
                  maxSize = 500,
                  nproc = 1)
  
  gsea_go <- gsea_go %>%
    arrange(pval) %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(head(x, 5), collapse=", ")))
  
  message("\n  Running GSEA on Hallmark...")
  gsea_hallmark <- fgsea(pathways = hallmark_list,
                        stats = gene_list,
                        minSize = 15,
                        nproc = 1)
  
  gsea_hallmark <- gsea_hallmark %>%
    arrange(pval) %>%
    mutate(leadingEdge = sapply(leadingEdge, function(x) paste(head(x, 5), collapse=", ")))
  
  # Print top results
  n_sig_go <- sum(gsea_go$padj < 0.05, na.rm = TRUE)
  n_sig_hallmark <- sum(gsea_hallmark$padj < 0.05, na.rm = TRUE)
  
  message(sprintf("\n  Significant pathways (FDR < 0.05):"))
  message(sprintf("    GO BP: %d", n_sig_go))
  message(sprintf("    Hallmark: %d", n_sig_hallmark))
  
  if (n_sig_go > 0) {
    message("\n  Top GO BP pathways:")
    top_go <- gsea_go %>%
      filter(padj < 0.05) %>%
      head(10) %>%
      select(pathway, NES, pval, padj, size)
    print(top_go)
  }
  
  if (n_sig_hallmark > 0) {
    message("\n  Top Hallmark pathways:")
    top_hallmark <- gsea_hallmark %>%
      filter(padj < 0.05) %>%
      head(10) %>%
      select(pathway, NES, pval, padj, size)
    print(top_hallmark)
  }
  
  return(list(
    go_bp = gsea_go,
    hallmark = gsea_hallmark,
    gene_ranks = gene_list,
    tissue = tissue_name
  ))
}

# Run GSEA for both tissues
library(fgsea)
library(msigdbr)

gsea_macular <- perform_gsea_variance(macular_var, "Macular")
gsea_extramacular <- perform_gsea_variance(extramacular_var, "Extramacular")

# Save GSEA results
write.csv(gsea_macular$go_bp,
          file.path(OUTPUT_DIR, "macular_variance_GSEA_GO_BP.csv"),
          row.names = FALSE)
write.csv(gsea_macular$hallmark,
          file.path(OUTPUT_DIR, "macular_variance_GSEA_Hallmark.csv"),
          row.names = FALSE)
write.csv(gsea_extramacular$go_bp,
          file.path(OUTPUT_DIR, "extramacular_variance_GSEA_GO_BP.csv"),
          row.names = FALSE)
write.csv(gsea_extramacular$hallmark,
          file.path(OUTPUT_DIR, "extramacular_variance_GSEA_Hallmark.csv"),
          row.names = FALSE)

# Visualize top GSEA results
if (!is.null(gsea_extramacular$hallmark)) {
  
  top_pathways <- gsea_extramacular$hallmark %>%
    filter(padj < 0.1) %>%
    arrange(pval) %>%
    head(20)
  
  if (nrow(top_pathways) > 0) {
    p_gsea <- ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
      geom_col() +
      coord_flip() +
      scale_fill_gradient(low = "red", high = "blue") +
      labs(
        title = "Extramacular: Variance-Associated Pathways (GSEA)",
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)",
        fill = "FDR"
      ) +
      theme(axis.text.y = element_text(size = 8))
    
    ggsave(file.path(OUTPUT_DIR, "extramacular_variance_GSEA_barplot.pdf"),
           p_gsea, width = 10, height = 8)
  }
}

#%%
# ========================================
# ANALYSIS 2: CHECK KNOWN AMD GENES
# ========================================

message("\n\n=== ANALYSIS 2: Known AMD Genes ===\n")

# Known AMD-related genes
amd_genes <- c(
  # Complement pathway
  "CFH", "C3", "C2", "CFB", "CFI", "C9",
  # ECM/drusen
  "HTRA1", "ARMS2", "EFEMP1", "TIMP3", "COL8A1",
  # Lipid metabolism
  "APOE", "CETP", "LIPC", "ABCA1",
  # Oxidative stress
  "SOD2", "NRF2", "NFE2L2", "HMOX1", "GPX1",
  # RPE-specific
  "RPE65", "BEST1", "CTNNA1",
  # Angiogenesis
  "VEGFA", "FLT1", "KDR"
)

check_amd_genes <- function(var_results, tissue_name) {
  
  message(sprintf("\nTissue: %s", tissue_name))
  
  amd_in_data <- var_results %>%
    filter(gene_symbol %in% amd_genes) %>%
    arrange(f_pvalue) %>%
    select(gene_symbol, mean_control, mean_amd, var_control, var_amd, 
           var_fc, f_pvalue, f_fdr, cv_control, cv_amd)
  
  if (nrow(amd_in_data) > 0) {
    message(sprintf("  Found %d / %d known AMD genes in data", nrow(amd_in_data), length(amd_genes)))
    message("\n  Top by variance difference:")
    print(head(amd_in_data, 10))
    
    # Count significant
    n_sig <- sum(amd_in_data$f_fdr < 0.05, na.rm = TRUE)
    message(sprintf("\n  Significant at FDR < 0.05: %d", n_sig))
    
  } else {
    message("  No known AMD genes found")
  }
  
  return(amd_in_data)
}

amd_macular <- check_amd_genes(macular_var, "Macular")
amd_extramacular <- check_amd_genes(extramacular_var, "Extramacular")

# Save AMD gene variance results
write.csv(amd_macular,
          file.path(OUTPUT_DIR, "macular_variance_AMD_genes.csv"),
          row.names = FALSE)
write.csv(amd_extramacular,
          file.path(OUTPUT_DIR, "extramacular_variance_AMD_genes.csv"),
          row.names = FALSE)

#%%
# ========================================
# ANALYSIS 3: OVERLAP WITH DE GENES
# ========================================

message("\n\n=== ANALYSIS 3: Overlap with DE Genes ===\n")

# Load DE results
de_results <- read.csv(file.path(DATA_DIR, "meta_combined_intermediate_annotated.csv"))

# High-variance genes (FDR < 0.05)
var_genes_mac <- macular_var %>%
  filter(f_fdr < 0.05 & !is.na(gene_symbol) & gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

var_genes_extra <- extramacular_var %>%
  filter(f_fdr < 0.05 & !is.na(gene_symbol) & gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

# DE genes (different thresholds for comparison)
de_genes_strict <- de_results %>%
  filter(meta_fdr < 0.10 & !is.na(gene_symbol) & gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

de_genes_relaxed <- de_results %>%
  filter(meta_p < 0.01 & abs(meta_logFC) > 0.5 & !is.na(gene_symbol) & gene_symbol != "") %>%
  pull(gene_symbol) %>%
  unique()

message(sprintf("Variance genes (macular, FDR<0.05): %d", length(var_genes_mac)))
message(sprintf("Variance genes (extramacular, FDR<0.05): %d", length(var_genes_extra)))
message(sprintf("DE genes (meta, FDR<0.10): %d", length(de_genes_strict)))
message(sprintf("DE genes (meta, p<0.01, |FC|>0.5): %d", length(de_genes_relaxed)))

# Calculate overlaps
overlap_extra_de_strict <- intersect(var_genes_extra, de_genes_strict)
overlap_extra_de_relaxed <- intersect(var_genes_extra, de_genes_relaxed)

message(sprintf("\nExtramacular Variance ∩ DE (FDR<0.10): %d genes", length(overlap_extra_de_strict)))
message(sprintf("Extramacular Variance ∩ DE (p<0.01): %d genes", length(overlap_extra_de_relaxed)))

if (length(overlap_extra_de_relaxed) > 0) {
  message("\nOverlapping genes:")
  print(overlap_extra_de_relaxed)
}

# Venn diagram
pdf(file.path(OUTPUT_DIR, "variance_vs_DE_venn.pdf"), width = 8, height = 8)

grid.newpage()
venn.plot <- venn.diagram(
  x = list(
    "High Variance\n(Extramacular)" = var_genes_extra,
    "DE Genes\n(p<0.01)" = de_genes_relaxed
  ),
  filename = NULL,
  col = "black",
  fill = c("steelblue", "coral"),
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontfamily = "sans",
  cat.pos = c(-20, 20),
  main = "High-Variance vs DE Genes",
  main.cex = 1.8
)
grid.draw(venn.plot)

dev.off()

# Statistical test: Is overlap more than expected by chance?
total_genes <- nrow(extramacular_var %>% filter(!is.na(gene_symbol) & gene_symbol != ""))
fisher_test <- fisher.test(matrix(c(
  length(overlap_extra_de_relaxed),  # both
  length(var_genes_extra) - length(overlap_extra_de_relaxed),  # var only
  length(de_genes_relaxed) - length(overlap_extra_de_relaxed),  # DE only
  total_genes - length(var_genes_extra) - length(de_genes_relaxed) + length(overlap_extra_de_relaxed)  # neither
), nrow = 2))

message(sprintf("\nFisher's exact test for overlap enrichment:"))
message(sprintf("  OR = %.2f, p = %.2e", fisher_test$estimate, fisher_test$p.value))

if (fisher_test$p.value < 0.05) {
  message("  Overlap is GREATER than expected by chance")
} else {
  message("  Overlap is NOT greater than expected by chance (orthogonal signals)")
}

#%%
# ========================================
# ANALYSIS 5: BIMODALITY TESTING
# ========================================

message("\n\n=== ANALYSIS 5: Bimodality Testing ===\n")

test_bimodality <- function(gene_expr_control, gene_expr_amd) {
  
  # Test bimodality in AMD samples only
  if (length(gene_expr_amd) < 10) {
    return(list(is_bimodal = FALSE, bic_diff = NA))
  }
  
  # Fit 1-component and 2-component mixture models
  tryCatch({
    fit1 <- Mclust(gene_expr_amd, G = 1, verbose = FALSE)
    fit2 <- Mclust(gene_expr_amd, G = 2, verbose = FALSE)
    
    if (is.null(fit1) || is.null(fit2)) {
      return(list(is_bimodal = FALSE, bic_diff = NA))
    }
    
    bic_diff <- fit2$bic - fit1$bic
    
    # BIC improvement > 10 suggests strong evidence for 2 components
    is_bimodal <- bic_diff > 10
    
    return(list(
      is_bimodal = is_bimodal,
      bic_diff = bic_diff,
      fit2 = fit2
    ))
  }, error = function(e) {
    return(list(is_bimodal = FALSE, bic_diff = NA))
  })
}

# Test bimodality for top high-variance genes
test_genes_bimodality <- function(var_results, expr_data, metadata, tissue_name, top_n = 50) {
  
  message(sprintf("\nTesting %s tissue...", tissue_name))
  
  # Get top high-variance genes
  top_var <- var_results %>%
    filter(f_fdr < 0.05) %>%
    arrange(f_pvalue) %>%
    head(top_n)
  
  if (nrow(top_var) == 0) {
    message("  No significant variance genes to test")
    return(NULL)
  }
  
  message(sprintf("  Testing %d genes", nrow(top_var)))
  
  # Prepare sample indices
  control_idx <- which(metadata$disease_status == "Control")
  amd_idx <- which(metadata$disease_status == "AMD")
  
  # Test each gene
  bimodality_results <- data.frame()
  
  for (i in 1:nrow(top_var)) {
    gene_id <- as.character(top_var$gene[i])
    gene_symbol <- top_var$gene_symbol[i]
    
    if (gene_id %in% colnames(expr_data)) {
      gene_expr <- expr_data[, gene_id]
      
      control_expr <- gene_expr[control_idx]
      amd_expr <- gene_expr[amd_idx]
      
      bimod_result <- test_bimodality(control_expr, amd_expr)
      
      bimodality_results <- rbind(bimodality_results, data.frame(
        gene = gene_id,
        gene_symbol = ifelse(is.na(gene_symbol), "", gene_symbol),
        var_fc = top_var$var_fc[i],
        f_fdr = top_var$f_fdr[i],
        is_bimodal = bimod_result$is_bimodal,
        bic_improvement = bimod_result$bic_diff
      ))
    }
  }
  
  # Count bimodal genes
  n_bimodal <- sum(bimodality_results$is_bimodal, na.rm = TRUE)
  message(sprintf("  Bimodal genes (BIC improvement > 10): %d / %d", n_bimodal, nrow(bimodality_results)))
  
  if (n_bimodal > 0) {
    message("\n  Bimodal genes:")
    print(bimodality_results %>% filter(is_bimodal) %>% arrange(desc(bic_improvement)))
  }
  
  return(bimodality_results)
}

# Test bimodality
bimod_macular <- test_genes_bimodality(macular_var, macular_data$expr, 
                                       macular_data$metadata, "Macular")
bimod_extramacular <- test_genes_bimodality(extramacular_var, extramacular_data$expr,
                                            extramacular_data$metadata, "Extramacular")

# Save bimodality results
if (!is.null(bimod_macular)) {
  write.csv(bimod_macular,
            file.path(OUTPUT_DIR, "macular_variance_bimodality.csv"),
            row.names = FALSE)
}
if (!is.null(bimod_extramacular)) {
  write.csv(bimod_extramacular,
            file.path(OUTPUT_DIR, "extramacular_variance_bimodality.csv"),
            row.names = FALSE)
}

#%%
# Summary
message("\n\n=== ANALYSIS COMPLETE ===")
message("Generated files:")
message("  - *_variance_GSEA_GO_BP.csv (GSEA results)")
message("  - *_variance_GSEA_Hallmark.csv (GSEA results)")
message("  - *_variance_AMD_genes.csv")
message("  - variance_vs_DE_venn.pdf")
message("  - *_variance_bimodality.csv")
message("  - *_variance_GSEA_barplot.pdf")
