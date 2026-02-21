#%%
# Comprehensive GSEA Analysis - GSE135092
# 1. Full C2 curated collection
# 2. GSEA on DE genes (mean-based)
# 3. Directional variance GSEA (increased vs decreased)

library(dplyr)
library(readr)
library(fgsea)
library(msigdbr)

#%%
# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")

cat("Comprehensive GSEA Analysis\n")
cat(rep("=", 70), "\n\n", sep = "")

#%%
# Load results with gene symbols
results_mac <- read_csv(file.path(results_dir, "RPE_Macula_all_results_symbols.csv"), 
                        show_col_types = FALSE)
results_nonmac <- read_csv(file.path(results_dir, "RPE_nonMacula_all_results_symbols.csv"),
                           show_col_types = FALSE)

cat(sprintf("Loaded: %d genes\n\n", nrow(results_mac)))

#%%
# Download MSigDB gene sets
cat("Downloading MSigDB gene sets...\n")

# Hallmark
msig_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets <- split(msig_hallmark$gene_symbol, msig_hallmark$gs_name)
cat(sprintf("  Hallmark: %d gene sets\n", length(hallmark_sets)))

# GO BP
msig_gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gobp_sets <- split(msig_gobp$gene_symbol, msig_gobp$gs_name)
cat(sprintf("  GO BP: %d gene sets\n", length(gobp_sets)))

# C2 - FULL curated collection (not just CP!)
msig_c2_full <- msigdbr(species = "Homo sapiens", collection = "C2")
c2_full_sets <- split(msig_c2_full$gene_symbol, msig_c2_full$gs_name)
cat(sprintf("  C2 Full Curated: %d gene sets\n\n", length(c2_full_sets)))

#%%
# Function to run GSEA
run_gsea <- function(results_df, gene_sets, tissue_name, db_name, gene_type, direction = NULL) {
  label <- ifelse(is.null(direction), gene_type, paste0(gene_type, "_", direction))
  cat(sprintf("Running: %s - %s - %s\n", tissue_name, db_name, label))
  
  # Filter based on gene type and direction
  if (gene_type == "variance") {
    if (!is.null(direction)) {
      if (direction == "increased") {
        filtered <- results_df %>% filter(!is.na(levene_pval), var_ratio > 1)
      } else {
        filtered <- results_df %>% filter(!is.na(levene_pval), var_ratio < 1)
      }
    } else {
      filtered <- results_df %>% filter(!is.na(levene_pval))
    }
    # Rank by variance p-value
    ranked_list <- filtered %>%
      mutate(rank_stat = -log10(levene_pval) * sign(var_ratio - 1)) %>%
      arrange(desc(rank_stat))
  } else { # DE genes
    filtered <- results_df %>% filter(!is.na(P.Value))
    # Rank by DE p-value and direction
    ranked_list <- filtered %>%
      mutate(rank_stat = -log10(P.Value) * sign(logFC)) %>%
      arrange(desc(rank_stat))
  }
  
  ranks <- setNames(ranked_list$rank_stat, ranked_list$gene_symbol)
  
  # Run GSEA
  gsea_res <- fgsea(pathways = gene_sets,
                    stats = ranks,
                    minSize = 15,
                    maxSize = 500)
  
  gsea_res <- gsea_res %>%
    filter(padj < 0.25) %>%
    arrange(padj, desc(abs(NES)))
  
  cat(sprintf("  Found %d pathways (FDR<0.25)\n", nrow(gsea_res)))
  
  return(gsea_res)
}

#%%
# 1. Full C2 Curated on Variance Genes
cat("\n1. C2 Full Curated on Variance Genes\n")
cat(rep("-", 70), "\n", sep = "")

c2_var_mac <- run_gsea(results_mac, c2_full_sets, "Macula", "C2_Full", "variance")
c2_var_nonmac <- run_gsea(results_nonmac, c2_full_sets, "non-Macula", "C2_Full", "variance")

write_csv(c2_var_mac, file.path(results_dir, "gsea/variance_GSEA_C2Full_macula.csv"))
write_csv(c2_var_nonmac, file.path(results_dir, "gsea/variance_GSEA_C2Full_nonmacula.csv"))

cat("\n")

#%%
# 2. GSEA on DE Genes
cat("2. GSEA on DE Genes (Mean-Based)\n")
cat(rep("-", 70), "\n", sep = "")

# Hallmark
de_hall_mac <- run_gsea(results_mac, hallmark_sets, "Macula", "Hallmark", "DE")
de_hall_nonmac <- run_gsea(results_nonmac, hallmark_sets, "non-Macula", "Hallmark", "DE")

write_csv(de_hall_mac, file.path(results_dir, "gsea/DE_GSEA_Hallmark_macula.csv"))
write_csv(de_hall_nonmac, file.path(results_dir, "gsea/DE_GSEA_Hallmark_nonmacula.csv"))

# GO BP
de_gobp_mac <- run_gsea(results_mac, gobp_sets, "Macula", "GO_BP", "DE")
de_gobp_nonmac <- run_gsea(results_nonmac, gobp_sets, "non-Macula", "GO_BP", "DE")

write_csv(de_gobp_mac, file.path(results_dir, "gsea/DE_GSEA_GOBP_macula.csv"))
write_csv(de_gobp_nonmac, file.path(results_dir, "gsea/DE_GSEA_GOBP_nonmacula.csv"))

# C2 Full
de_c2_mac <- run_gsea(results_mac, c2_full_sets, "Macula", "C2_Full", "DE")
de_c2_nonmac <- run_gsea(results_nonmac, c2_full_sets, "non-Macula", "C2_Full", "DE")

write_csv(de_c2_mac, file.path(results_dir, "gsea/DE_GSEA_C2Full_macula.csv"))
write_csv(de_c2_nonmac, file.path(results_dir, "gsea/DE_GSEA_C2Full_nonmacula.csv"))

cat("\n")

#%%
# 3. Directional Variance GSEA
cat("3. Directional Variance GSEA (Increased vs Decreased)\n")
cat(rep("-", 70), "\n", sep = "")

# Hallmark - Increased variance
var_inc_hall_mac <- run_gsea(results_mac, hallmark_sets, "Macula", "Hallmark", "variance", "increased")
var_inc_hall_nonmac <- run_gsea(results_nonmac, hallmark_sets, "non-Macula", "Hallmark", "variance", "increased")

write_csv(var_inc_hall_mac, file.path(results_dir, "gsea/variance_increased_GSEA_Hallmark_macula.csv"))
write_csv(var_inc_hall_nonmac, file.path(results_dir, "gsea/variance_increased_GSEA_Hallmark_nonmacula.csv"))

# Hallmark - Decreased variance
var_dec_hall_mac <- run_gsea(results_mac, hallmark_sets, "Macula", "Hallmark", "variance", "decreased")
var_dec_hall_nonmac <- run_gsea(results_nonmac, hallmark_sets, "non-Macula", "Hallmark", "variance", "decreased")

write_csv(var_dec_hall_mac, file.path(results_dir, "gsea/variance_decreased_GSEA_Hallmark_macula.csv"))
write_csv(var_dec_hall_nonmac, file.path(results_dir, "gsea/variance_decreased_GSEA_Hallmark_nonmacula.csv"))

# GO BP - Increased variance
var_inc_gobp_mac <- run_gsea(results_mac, gobp_sets, "Macula", "GO_BP", "variance", "increased")
var_inc_gobp_nonmac <- run_gsea(results_nonmac, gobp_sets, "non-Macula", "GO_BP", "variance", "increased")

write_csv(var_inc_gobp_mac, file.path(results_dir, "gsea/variance_increased_GSEA_GOBP_macula.csv"))
write_csv(var_inc_gobp_nonmac, file.path(results_dir, "gsea/variance_increased_GSEA_GOBP_nonmacula.csv"))

# GO BP - Decreased variance
var_dec_gobp_mac <- run_gsea(results_mac, gobp_sets, "Macula", "GO_BP", "variance", "decreased")
var_dec_gobp_nonmac <- run_gsea(results_nonmac, gobp_sets, "non-Macula", "GO_BP", "variance", "decreased")

write_csv(var_dec_gobp_mac, file.path(results_dir, "gsea/variance_decreased_GSEA_GOBP_macula.csv"))
write_csv(var_dec_gobp_nonmac, file.path(results_dir, "gsea/variance_decreased_GSEA_GOBP_nonmacula.csv"))

cat("\n")

#%%
# 4. Create summary table
cat("4. Creating Summary Table\n")
cat(rep("-", 70), "\n", sep = "")

summary_df <- data.frame(
  Tissue = rep(c("RPE Macula", "RPE non-Macula"), each = 9),
  Analysis = rep(c(
    "Variance (All)", "Variance (Increased)", "Variance (Decreased)",
    "DE (Mean)", "Variance C2 Full", "DE C2 Full",
    "Variance GO BP", "DE GO BP", "DE Hallmark"
  ), 2),
  Hallmark = c(
    10, nrow(var_inc_hall_mac), nrow(var_dec_hall_mac), nrow(de_hall_mac), NA, NA, NA, NA, nrow(de_hall_mac),
    18, nrow(var_inc_hall_nonmac), nrow(var_dec_hall_nonmac), nrow(de_hall_nonmac), NA, NA, NA, NA, nrow(de_hall_nonmac)
  ),
  GO_BP = c(
    116, nrow(var_inc_gobp_mac), nrow(var_dec_gobp_mac), nrow(de_gobp_mac), NA, NA, nrow(var_inc_gobp_mac) + nrow(var_dec_gobp_mac), nrow(de_gobp_mac), NA,
    514, nrow(var_inc_gobp_nonmac), nrow(var_dec_gobp_nonmac), nrow(de_gobp_nonmac), NA, NA, nrow(var_inc_gobp_nonmac) + nrow(var_dec_gobp_nonmac), nrow(de_gobp_nonmac), NA
  ),
  C2_Full = c(
    nrow(c2_var_mac), NA, NA, nrow(de_c2_mac), nrow(c2_var_mac), nrow(de_c2_mac), NA, NA, NA,
    nrow(c2_var_nonmac), NA, NA, nrow(de_c2_nonmac), nrow(c2_var_nonmac), nrow(de_c2_nonmac), NA, NA, NA
  )
)

write_csv(summary_df, file.path(results_dir, "GSEA_comprehensive_summary.csv"))

cat("\n")
print(summary_df)

cat("\n", rep("=", 70), "\n", sep = "")
cat("COMPREHENSIVE GSEA COMPLETE!\n")
cat(rep("=", 70), "\n", sep = "")

cat("\nFiles saved in:", file.path(results_dir, "gsea/"), "\n")
cat("Summary table:", file.path(results_dir, "GSEA_comprehensive_summary.csv"), "\n")
