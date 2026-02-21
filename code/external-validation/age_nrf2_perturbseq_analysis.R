#%% ============================================================================
# Age-NRF2/ARE Regulatory Connection via Perturb-seq
# Analysis connecting age-related DE genes (GSE29801) to NRF2/ARE pathway
# through Perturb-seq causal mapping
# Using data.table for data manipulation
# ============================================================================

library(data.table)
library(fgsea)
library(ggplot2)

# Set paths
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE29801/age_nrf2_analysis")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Visualization settings (per coding preferences - using default font if Palatino not available)
theme_set(theme_classic())

#%% ============================================================================
# PHASE 1: DATA PREPARATION
# ============================================================================

cat("=== PHASE 1: DATA PREPARATION ===\n\n")

#%% Load Age-DE genes from GSE29801
cat("Loading age-DE genes from GSE29801...\n")

# Extramacular
age_extra <- fread(file.path(base_dir, "results/cohort-GSE29801/dry_amd_de/extramacular_age_age_main.csv"))
# The 'gene' column contains probe IDs, not gene symbols
# Load probe annotation to map to gene symbols
probe_annot <- fread(file.path(base_dir, "results/cohort-GSE29801/dry_amd_de/probe_annotation_complete.csv"))
# Get gene symbol mapping (use variance annotated file as it has gene_symbol)
var_annot <- fread(file.path(base_dir, "results/cohort-GSE29801/dry_amd_de/extramacular_variance_analysis_annotated.csv"),
                   select = c("gene", "gene_symbol"))
var_annot <- unique(var_annot[!is.na(gene_symbol) & gene_symbol != ""])

# Merge gene symbols into age-DE data
setnames(age_extra, ncol(age_extra), "probe_id")
age_extra <- merge(age_extra, var_annot, by.x = "probe_id", by.y = "gene", all.x = TRUE)
setnames(age_extra, "gene_symbol", "gene")

# Macular  
age_mac <- fread(file.path(base_dir, "results/cohort-GSE29801/dry_amd_de/macular_age_age_main.csv"))
setnames(age_mac, ncol(age_mac), "probe_id")
age_mac <- merge(age_mac, var_annot, by.x = "probe_id", by.y = "gene", all.x = TRUE)
setnames(age_mac, "gene_symbol", "gene")

cat(sprintf("  Extramacular: %d probes/genes\n", nrow(age_extra)))
cat(sprintf("  Macular: %d probes/genes\n", nrow(age_mac)))

# Filter significant age-DE genes and separate by direction
age_extra_sig <- age_extra[adj.P.Val < 0.05 & !is.na(gene) & gene != ""]
age_extra_sig[, `:=`(
  direction = ifelse(logFC > 0, "UP", "DOWN"),
  tissue = "extramacular"
)]

age_mac_sig <- age_mac[adj.P.Val < 0.05 & !is.na(gene) & gene != ""]
age_mac_sig[, `:=`(
  direction = ifelse(logFC > 0, "UP", "DOWN"),
  tissue = "macular"
)]

cat(sprintf("  Extramacular significant (FDR<0.05): %d genes (%d UP, %d DOWN)\n", 
            nrow(age_extra_sig), 
            sum(age_extra_sig$direction == "UP"), 
            sum(age_extra_sig$direction == "DOWN")))
cat(sprintf("  Macular significant (FDR<0.05): %d genes (%d UP, %d DOWN)\n", 
            nrow(age_mac_sig),
            sum(age_mac_sig$direction == "UP"),
            sum(age_mac_sig$direction == "DOWN")))

# Combine for complete age-DE gene list
age_de_all <- rbind(age_extra_sig, age_mac_sig)

#%% Load PRG4 rescue signature
cat("\nLoading PRG4 rescue signature...\n")

prg4_de <- fread(file.path(base_dir, "data/RPE_DE_results.csv"))

# Clean up - get gene symbols
prg4_de <- prg4_de[!is.na(hgnc_symbol) & hgnc_symbol != ""]
setnames(prg4_de, "hgnc_symbol", "gene")

# Create rescue signature (H2O2+PRG4 vs H2O2 = PRG4 rescue effect)
prg4_rescue <- prg4_de[!is.na(H2O2PRG4_vs_H2O2.log2FoldChange), 
                        .(gene, 
                          rescue_lfc = H2O2PRG4_vs_H2O2.log2FoldChange, 
                          rescue_pval = H2O2PRG4_vs_H2O2.pvalue,
                          h2o2_lfc = H2O2_vs_CTRL.log2FoldChange,
                          h2o2_pval = H2O2_vs_CTRL.pvalue)]
prg4_rescue[, abs_rescue := abs(rescue_lfc)]
setorder(prg4_rescue, -abs_rescue)
prg4_rescue[, abs_rescue := NULL]

cat(sprintf("  PRG4 rescue signature: %d genes\n", nrow(prg4_rescue)))

#%% Define ARE/NRF2 gene set
cat("\nDefining ARE/NRF2 gene set...\n")

# Load MSigDB gene sets
gmt_file <- file.path(base_dir, "data/c2.all.v2023.2.Hs.symbols.gmt")
if(!file.exists(gmt_file)) {
  stop("MSigDB gmt file not found: ", gmt_file)
}

msigdb <- gmtPathways(gmt_file)

# Get ROS/NRF2-related gene sets
ros_related <- names(msigdb)[grepl("REACTIVE_OXYGEN|NRF2|NFE2L2|OXIDATIVE_STRESS|ANTIOXIDANT", 
                                    names(msigdb), ignore.case = TRUE)]
cat("  Found ROS/NRF2-related gene sets:\n")
for(gs in ros_related) {
  cat(sprintf("    - %s (%d genes)\n", gs, length(msigdb[[gs]])))
}

# Define core ARE gene set from literature
known_ARE_genes <- c(
  # Core NRF2 targets
  "NFE2L2", "KEAP1", "NQO1", "HMOX1", "GCLC", "GCLM", "GSR", "GPX1", "GPX2", 
  "TXNRD1", "SRXN1", "PRDX1", "PRDX6", "SOD1", "SOD2", "CAT",
  # Phase II enzymes
  "GSTA1", "GSTA2", "GSTM1", "GSTP1", "UGT1A1", "ABCC1", "ABCC2",
  # Thioredoxin system
  "TXN", "TXNRD2", "PRDX2", "PRDX3", "PRDX4", "PRDX5",
  # Other ARE targets
  "FTH1", "FTL", "SLC7A11", "ME1", "PGD", "G6PD", "MAFG", "MAFK"
)

# Combine with MSigDB if available
if(length(ros_related) > 0) {
  ARE_genes <- unique(c(known_ARE_genes, unlist(msigdb[ros_related[1:min(3, length(ros_related))]])))
} else {
  ARE_genes <- known_ARE_genes
}

cat(sprintf("  ARE gene set: %d genes total\n", length(ARE_genes)))

# Save ARE gene set
fwrite(data.table(gene = ARE_genes), file.path(results_dir, "ARE_gene_set.csv"))

#%% Check overlap with age-DE genes
cat("\nChecking ARE overlap with age-DE genes...\n")

are_in_age_extra <- intersect(ARE_genes, age_extra_sig$gene)
are_in_age_mac <- intersect(ARE_genes, age_mac_sig$gene)

cat(sprintf("  ARE genes in extramacular age-DE: %d / %d\n", 
            length(are_in_age_extra), length(ARE_genes)))
cat(sprintf("  ARE genes in macular age-DE: %d / %d\n", 
            length(are_in_age_mac), length(ARE_genes)))

if(length(are_in_age_extra) > 0) {
  cat("  Extramacular ARE genes:\n")
  for(g in are_in_age_extra) {
    row <- age_extra_sig[gene == g][1]
    cat(sprintf("    %s: %s (logFC=%.3f)\n", g, row$direction, row$logFC))
  }
}

#%% Create age signature vectors (for GSEA)
cat("\nCreating age signature vectors...\n")

# Use t-statistic as ranking metric
age_extra_clean <- age_extra[!is.na(gene) & gene != "" & !is.na(t)]
age_extra_vec <- setNames(age_extra_clean$t, age_extra_clean$gene)
age_extra_vec <- sort(age_extra_vec, decreasing = TRUE)

age_mac_clean <- age_mac[!is.na(gene) & gene != "" & !is.na(t)]
age_mac_vec <- setNames(age_mac_clean$t, age_mac_clean$gene)
age_mac_vec <- sort(age_mac_vec, decreasing = TRUE)

cat(sprintf("  Extramacular signature: %d genes\n", length(age_extra_vec)))
cat(sprintf("  Macular signature: %d genes\n", length(age_mac_vec)))

#%% Save prepared data
cat("\nSaving prepared data...\n")

fwrite(age_extra_sig, file.path(results_dir, "age_de_extramacular.csv"))
fwrite(age_mac_sig, file.path(results_dir, "age_de_macular.csv"))
fwrite(age_de_all, file.path(results_dir, "age_de_all.csv"))
fwrite(prg4_rescue, file.path(results_dir, "prg4_rescue_signature.csv"))

saveRDS(list(
  extramacular = age_extra_vec,
  macular = age_mac_vec
), file.path(results_dir, "age_signatures.rds"))

cat("\n=== PHASE 1 COMPLETE ===\n")

#%% ============================================================================
# PHASE 2D: PRG4 TRIANGULATION
# Three explicit hypotheses
# ============================================================================

cat("\n=== PHASE 2D: PRG4 TRIANGULATION ===\n\n")

#%% Hypothesis 1: Age ↔ PRG4 correlation
cat("--- Hypothesis 1: Age ↔ PRG4 ---\n")
cat("Question: Does PRG4 rescue reverse the age signature?\n\n")

# Get overlapping genes
common_genes_prg4 <- intersect(names(age_extra_vec), prg4_rescue$gene)
cat(sprintf("Overlapping genes (age x PRG4): %d\n", length(common_genes_prg4)))

if(length(common_genes_prg4) > 100) {
  # Create matched vectors
  age_vals <- age_extra_vec[common_genes_prg4]
  prg4_lfc_vec <- setNames(prg4_rescue$rescue_lfc, prg4_rescue$gene)
  prg4_vals <- prg4_lfc_vec[common_genes_prg4]
  
  # Correlation
  cor_test <- cor.test(age_vals, prg4_vals, method = "spearman")
  
  cat(sprintf("\nAge vs PRG4 rescue correlation:\n"))
  cat(sprintf("  Spearman rho = %.3f\n", cor_test$estimate))
  cat(sprintf("  p-value = %.2e\n", cor_test$p.value))
  
  h1_interpretation <- if(cor_test$estimate < 0 & cor_test$p.value < 0.05) {
    "SUPPORTED: PRG4 reverses age signature"
  } else if(cor_test$estimate > 0 & cor_test$p.value < 0.05) {
    "OPPOSITE: PRG4 enhances age signature"
  } else {
    "NOT SIGNIFICANT"
  }
  
  h1_result <- data.table(
    hypothesis = "H1: Age vs PRG4",
    statistic = cor_test$estimate,
    pvalue = cor_test$p.value,
    n = length(common_genes_prg4),
    interpretation = h1_interpretation
  )
  print(h1_result)
  
  # Create scatter plot
  plot_data_h1 <- data.table(
    gene = common_genes_prg4,
    age_t = age_vals,
    prg4_lfc = prg4_vals,
    is_ARE = common_genes_prg4 %in% ARE_genes
  )
  
  p_h1 <- ggplot(plot_data_h1, aes(x = age_t, y = prg4_lfc)) +
    geom_point(aes(color = is_ARE), alpha = 0.5, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "#D35400"),
                       labels = c("Other", "ARE gene")) +
    labs(
      title = "H1: Age vs PRG4 Rescue",
      subtitle = sprintf("Spearman rho = %.3f, p = %.2e", cor_test$estimate, cor_test$p.value),
      x = "Age effect (t-statistic)",
      y = "PRG4 rescue (log2FC)",
      color = ""
    ) +
    theme(
      axis.text = element_text(color = "black"),
      legend.position = "bottom"
    )
  
  ggsave(file.path(results_dir, "H1_age_vs_prg4_scatter.pdf"), p_h1, 
         width = 5, height = 4.5)
  ggsave(file.path(results_dir, "H1_age_vs_prg4_scatter.tiff"), p_h1,
         width = 5, height = 4.5, dpi = 300, compression = "lzw")
  cat("  Saved: H1_age_vs_prg4_scatter.pdf/tiff\n")
  
} else {
  cat("Insufficient overlapping genes for correlation\n")
  h1_result <- data.table(
    hypothesis = "H1: Age vs PRG4",
    statistic = NA_real_, pvalue = NA_real_, n = length(common_genes_prg4),
    interpretation = "INSUFFICIENT DATA"
  )
}

#%% Hypothesis 2: Age ↔ ARE
cat("\n--- Hypothesis 2: Age ↔ ARE ---\n")
cat("Question: Are ARE genes dysregulated with age?\n\n")

# GSEA of ARE genes in age signature
ARE_list <- list(ARE = ARE_genes)

# Check overlap first
are_in_signature <- sum(ARE_genes %in% names(age_extra_vec))
cat(sprintf("ARE genes in signature: %d / %d\n", are_in_signature, length(ARE_genes)))

gsea_age_are <- fgsea(
  pathways = ARE_list,
  stats = age_extra_vec,
  minSize = 5,
  maxSize = 500
)

cat("GSEA: ARE genes in age signature (extramacular):\n")

if(nrow(gsea_age_are) > 0 && !is.na(gsea_age_are$NES[1])) {
  print(gsea_age_are[, .(pathway, pval, padj, NES, size)])
  
  h2_interpretation <- if(gsea_age_are$NES < 0 & gsea_age_are$padj < 0.05) {
    "SUPPORTED: ARE genes DOWN with age"
  } else if(gsea_age_are$NES > 0 & gsea_age_are$padj < 0.05) {
    "OPPOSITE: ARE genes UP with age"
  } else {
    "NOT SIGNIFICANT"
  }
  
  h2_result <- data.table(
    hypothesis = "H2: Age vs ARE",
    statistic = gsea_age_are$NES,
    pvalue = gsea_age_are$padj,
    n = gsea_age_are$size,
    interpretation = h2_interpretation
  )
  print(h2_result)
  
  # GSEA plot
  pdf(file.path(results_dir, "H2_age_ARE_gsea.pdf"), width = 5, height = 4)
  print(
    plotEnrichment(ARE_list[["ARE"]], age_extra_vec) +
      labs(title = "H2: ARE genes in Age Signature",
           subtitle = sprintf("NES = %.2f, padj = %.3f", gsea_age_are$NES, gsea_age_are$padj)) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"))
  )
  dev.off()
  cat("  Saved: H2_age_ARE_gsea.pdf\n")
} else {
  cat("  GSEA returned no results (insufficient overlap or other issue)\n")
  h2_result <- data.table(
    hypothesis = "H2: Age vs ARE",
    statistic = NA_real_,
    pvalue = NA_real_,
    n = are_in_signature,
    interpretation = "INSUFFICIENT OVERLAP"
  )
  print(h2_result)
}

#%% Hypothesis 3: PRG4 ↔ ARE
cat("\n--- Hypothesis 3: PRG4 ↔ ARE ---\n")
cat("Question: Does PRG4 upregulate ARE genes?\n\n")

# Create PRG4 rescue signature vector
prg4_vec <- setNames(prg4_rescue$rescue_lfc, prg4_rescue$gene)
prg4_vec <- prg4_vec[!is.na(prg4_vec)]
prg4_vec <- sort(prg4_vec, decreasing = TRUE)

# Check overlap
are_in_prg4 <- sum(ARE_genes %in% names(prg4_vec))
cat(sprintf("ARE genes in PRG4 signature: %d / %d\n", are_in_prg4, length(ARE_genes)))

# GSEA
gsea_prg4_are <- fgsea(
  pathways = ARE_list,
  stats = prg4_vec,
  minSize = 5,
  maxSize = 500
)

cat("GSEA: ARE genes in PRG4 rescue signature:\n")

if(nrow(gsea_prg4_are) > 0 && !is.na(gsea_prg4_are$NES[1])) {
  print(gsea_prg4_are[, .(pathway, pval, padj, NES, size)])
  
  h3_interpretation <- if(gsea_prg4_are$NES > 0 & gsea_prg4_are$padj < 0.05) {
    "SUPPORTED: PRG4 upregulates ARE genes"
  } else if(gsea_prg4_are$NES < 0 & gsea_prg4_are$padj < 0.05) {
    "OPPOSITE: PRG4 downregulates ARE genes"
  } else {
    "NOT SIGNIFICANT"
  }
  
  h3_result <- data.table(
    hypothesis = "H3: PRG4 vs ARE",
    statistic = gsea_prg4_are$NES,
    pvalue = gsea_prg4_are$padj,
    n = gsea_prg4_are$size,
    interpretation = h3_interpretation
  )
  print(h3_result)
  
  # GSEA plot
  pdf(file.path(results_dir, "H3_prg4_ARE_gsea.pdf"), width = 5, height = 4)
  print(
    plotEnrichment(ARE_list[["ARE"]], prg4_vec) +
      labs(title = "H3: ARE genes in PRG4 Rescue Signature",
           subtitle = sprintf("NES = %.2f, padj = %.3f", gsea_prg4_are$NES, gsea_prg4_are$padj)) +
      theme_classic() +
      theme(axis.text = element_text(color = "black"))
  )
  dev.off()
  cat("  Saved: H3_prg4_ARE_gsea.pdf\n")
} else {
  cat("  GSEA returned no results (insufficient overlap or other issue)\n")
  h3_result <- data.table(
    hypothesis = "H3: PRG4 vs ARE",
    statistic = NA_real_,
    pvalue = NA_real_,
    n = are_in_prg4,
    interpretation = "INSUFFICIENT OVERLAP"
  )
  print(h3_result)
}

#%% Triangulation summary
cat("\n=== TRIANGULATION SUMMARY ===\n\n")

triangulation_results <- rbind(h1_result, h2_result, h3_result)
print(triangulation_results)

fwrite(triangulation_results, file.path(results_dir, "prg4_triangulation_results.csv"))

# Overall interpretation
cat("\n")
supported_count <- sum(grepl("SUPPORTED", triangulation_results$interpretation))
if(supported_count == 3) {
  cat("*** ALL THREE HYPOTHESES SUPPORTED ***\n")
  cat("Interpretation: Aging is associated with decline of antioxidant response.\n")
  cat("PRG4 restores antioxidant capacity, reversing age-associated transcriptomic damage.\n")
} else if(supported_count >= 1) {
  cat(sprintf("*** %d of 3 HYPOTHESES SUPPORTED ***\n", supported_count))
  cat("Mixed results - see individual hypothesis interpretations above.\n")
} else {
  cat("*** NO HYPOTHESES SUPPORTED ***\n")
  cat("The triangulation model is not supported by this data.\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat(sprintf("All outputs saved to: %s\n", results_dir))

# List output files
cat("\nOutput files:\n")
for(f in list.files(results_dir)) {
  cat(sprintf("  - %s\n", f))
}
