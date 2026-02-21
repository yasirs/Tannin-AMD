
source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)


# Constants
tpm_file <- "data/RPE_TPMS.csv"
gene_sets_dir <- "results/three-axes/gene-sets"
out_dir <- "results/three-axes/heatmaps"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load TPMs
tpms <- read_csv(tpm_file, show_col_types = FALSE)


# Assume structure: gene_symbol and sample columns (CTRL_1, ..., H2O2_1..., PRG4_1...)
# I need to inspect TPM columns. I'll assume standard format or check.

load_genes <- function(name) {
  f <- file.path(gene_sets_dir, paste0(name, "_mapped.csv"))
  if (file.exists(f)) read_csv(f, show_col_types = FALSE)$gene_symbol else c()
}

axes <- list(
  Oxidative = load_genes("oxidative_pro_disease"),
  Inflammatory = load_genes("inflammatory_pro_disease"),
  Senescence = load_genes("senescence_pro_disease")
)

# Function to plot heatmap for an axis
plot_axis_heatmap <- function(axis_name, gene_list) {
  # Filter TPMs
  gene_col <- "hgnc_symbol"
  # Also remove ensembl_gene_id if present to avoid pivot error
  tpm_subset <- tpms %>% select(-any_of("ensembl_gene_id"))
  
  axis_tpm <- tpm_subset %>% filter(!!sym(gene_col) %in% gene_list)
  
  if (nrow(axis_tpm) < 2) return(NULL)
  
  # Z-score normalization
  # Pivot longer
  long_dat <- axis_tpm %>%
    pivot_longer(cols = -!!sym(gene_col), names_to = "Sample", values_to = "TPM")
    
  # Calculate Z-score per gene
  long_dat <- long_dat %>%
    group_by(!!sym(gene_col)) %>%
    mutate(Z = (TPM - mean(TPM)) / sd(TPM)) %>%
    ungroup()
    
  # Order genes by H2O2 vs CTRL difference (approximate by mean diff)
  mat <- axis_tpm %>% select(-!!sym(gene_col)) %>% as.matrix()
  rownames(mat) <- axis_tpm[[gene_col]]
  
  # Handle NA/NaN if SD is 0 or matrix bad
  if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)
  
  dist_mat <- dist(mat)
  if (any(is.na(dist_mat))) return(NULL)
  hc <- hclust(dist_mat)
  gene_order <- rownames(mat)[hc$order]
  
  # Sample Mapping
  # TS1, TS2, TS3 -> CTRL
  # TS5, TS6, TS7 -> H2O2
  # TS9, TS11, TS12 -> PRG4
  # TS13, TS15, TS16 -> H2O2PRG4
  
  long_dat$Condition <- case_when(
    long_dat$Sample %in% c("TS1", "TS2", "TS3") ~ "CTRL",
    long_dat$Sample %in% c("TS5", "TS6", "TS7") ~ "H2O2",
    long_dat$Sample %in% c("TS9", "TS11", "TS12") ~ "PRG4",
    long_dat$Sample %in% c("TS13", "TS15", "TS16") ~ "H2O2+PRG4",
    TRUE ~ "Other" # Handling unexpected columns
  )
  
  # Filter out "Other" samples if any
  long_dat <- long_dat %>% filter(Condition != "Other")
  
  long_dat$Condition <- factor(long_dat$Condition, levels = c("CTRL", "H2O2", "PRG4", "H2O2+PRG4"))
  
  p <- ggplot(long_dat, aes(x = Sample, y = factor(!!sym(gene_col), levels = gene_order), fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#8c510a", mid = "#f5f5f5", high = "#01665e", midpoint = 0) + # BrBG-like reversed? 
    # Plan says: "BrBG (brown = high)" -> High should be Brown.
    # standard BrBG: Brown is low? No, Brown is usually negative in Diverging.
    # "Brown for high values" -> Okay, so High = Brown.
    # #8c510a is Brown. #01665e is Green.
    # So High (#8c510a), Low (#01665e).
    scale_fill_gradient2(low = "#01665e", high = "#8c510a", mid = "white", midpoint = 0, name = "Z-Score") +
    theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste(axis_name, "Axis Expression"), y = "Genes", x = "")
    
  ggsave(file.path(out_dir, paste0("Fig_Heatmap_", axis_name, ".pdf")), p, width = 6, height = 8)
}

for (ax in names(axes)) {
  plot_axis_heatmap(ax, axes[[ax]])
}
