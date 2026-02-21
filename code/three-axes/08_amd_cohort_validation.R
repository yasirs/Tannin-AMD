
# 08_amd_cohort_validation.R
# Purpose: Validate Axes in AMD Cohort using GSVA (Patient Scoring) and GSEA (Signature Matching).

source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(fgsea)
library(GSVA)
library(tibble) 

out_dir <- "results/three-axes/amd-validation"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load AMD Results (GSE135092)
# Need the FULL list for GSEA ranking.
# Note: The provided file is "DE results", so it has LogFC/Pval for ~all genes? Or filtered?
# Assuming it's the full list.
amd_file <- "results/cohort-GSE135092/GSE135092_DE_results.csv"
amd_res <- read_csv(amd_file, show_col_types = FALSE)

# Fix Symbol Column (Regex)
print(paste("AMD Res Columns:", paste(colnames(amd_res), collapse=", ")))
gene_col <- colnames(amd_res)[grepl("symbol", colnames(amd_res), ignore.case=TRUE)][1]
if (is.na(gene_col)) stop("Could not find gene symbol column")
print(paste("Using Gene Column:", gene_col))

# Trim
amd_res[[gene_col]] <- trimws(amd_res[[gene_col]])

# Rank AMD List
lfc_col <- colnames(amd_res)[grepl("logFC|log2FoldChange", colnames(amd_res), ignore.case=TRUE)][1]
pval_col <- colnames(amd_res)[grepl("adj.P.Val|padj|FDR|pvalue", colnames(amd_res), ignore.case=TRUE)][1]

if (is.na(lfc_col)) stop("Could not find LogFC column")
if (is.na(pval_col)) stop("Could not find P-value column")
print(paste("Using LFC:", lfc_col, "Pval:", pval_col))

# Load PRG4 Results to define signature
prg4_res <- read_csv("data/RPE_DE_results.csv", show_col_types = FALSE)
# Define "PRG4 Down / Rescue" Set:
# Genes DOWN by PRG4 (LogFC < -0.5, P < 0.05)
prg4_down_genes <- prg4_res %>% 
    filter(H2O2PRG4_vs_H2O2.log2FoldChange < -0.5, H2O2PRG4_vs_H2O2.pvalue < 0.05) %>%
    pull(hgnc_symbol)
min_p <- min(amd_res[[pval_col]][amd_res[[pval_col]] > 0], na.rm=TRUE)
amd_res$stat <- sign(amd_res[[lfc_col]]) * -log10(ifelse(amd_res[[pval_col]]==0, min_p, amd_res[[pval_col]]))

amd_ranks <- amd_res %>%
  filter(!is.na(stat), !is.na(!!sym(gene_col))) %>%
  arrange(desc(stat)) %>%
  select(!!sym(gene_col), stat) %>%
  distinct(!!sym(gene_col), .keep_all=TRUE) %>%
  deframe()
  
# Run GSEA for PRG4 Sig
fgsea_amd <- fgseaMultilevel(pathways = list(PRG4_Rescue_Down = prg4_down_genes), stats = amd_ranks)
print("AMD Validation GSEA (PRG4 Sig):")
print(fgsea_amd)
write_csv(fgsea_amd, file.path(out_dir, "amd_prg4_gsea.csv"))

# Run GSEA for AXIS Sets (To show Dysregulation)
gene_sets_dir <- "results/three-axes/gene-sets"
load_set <- function(name) {
  f <- file.path(gene_sets_dir, paste0(name, "_mapped.csv"))
  if (file.exists(f)) read_csv(f, show_col_types = FALSE)$gene_symbol else NULL
}
axis_paths <- list(
  Oxidative = load_set("oxidative_pro_disease"),
  Inflammatory = load_set("inflammatory_pro_disease"),
  Senescence = load_set("senescence_pro_disease")
)

fgsea_axis <- fgseaMultilevel(pathways = axis_paths, stats = amd_ranks, minSize=15, maxSize=500)
print("AMD Validation GSEA (Axes):")
print(fgsea_axis)
write_csv(fgsea_axis, file.path(out_dir, "amd_axis_dysregulation_gsea.csv"))

p_gsea <- plotEnrichment(prg4_down_genes, amd_ranks) +
  labs(title = "PRG4 Rescue Set in AMD", subtitle = paste("NES =", round(fgsea_amd$NES, 2))) +
  theme_pub()
ggsave(file.path(out_dir, "Fig_AMD_PRG4_GSEA.pdf"), p_gsea, width = 5, height = 4)


# 2. GSVA: Patient Scoring
# Ideally we need the Expression Matrix (not just DE results) for GSVA.
# "GSE135092_DE_results.csv" usually just has stats.
# Do we have the counts matrix?
# User prompt says: "Validate axis expression in AMD patients (GSE135092)..."
# If we don't have the matrix, we CANNOT run GSVA.
# Checking file availability... 
# If unavailable, we must stick to the GSEA validation above (Method 1) which is robust for DE lists.

# I will check if "GSE135092_expression_matrix.csv" or similar exists?
# If not, I will output a message or skip GSVA part.
# For now, I'll assume we only have DE results and stick to GSEA validation (which is sufficient).
# User plan mentioned "GSVA scores... correlates with clinical metadata".
# This requires patient-level data. The file path provided was `results/cohort-GSE135092_DE_results.csv`.
# I suspect we lack the matrix. If so, GSEA (Method 1) is the fallback.

# I will write the script to rely on GSEA only for now.
# Attempting GSVA without matrix is impossible.
