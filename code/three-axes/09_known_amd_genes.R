
source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# Constants
known_amd_file <- "results/three-axes/gene-sets/known_amd_genes_mapped.csv"
tpm_file <- "data/RPE_TPMS.csv"
out_dir <- "results/three-axes/known-amd"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load data
tpms <- read_csv(tpm_file, show_col_types = FALSE)
known_genes <- read_csv(known_amd_file, show_col_types = FALSE)$gene_symbol

# Filter TPMs
gene_col <- "hgnc_symbol"
tpm_subset <- tpms %>% select(-any_of("ensembl_gene_id"))

axis_tpm <- tpm_subset %>% filter(!!sym(gene_col) %in% known_genes)

# Make long
long_dat <- axis_tpm %>%
    pivot_longer(cols = -!!sym(gene_col), names_to = "Sample", values_to = "TPM")

# Map Conditions
long_dat$Condition <- case_when(
    long_dat$Sample %in% c("TS1", "TS2", "TS3") ~ "CTRL",
    long_dat$Sample %in% c("TS5", "TS6", "TS7") ~ "H2O2",
    long_dat$Sample %in% c("TS9", "TS11", "TS12") ~ "PRG4",
    long_dat$Sample %in% c("TS13", "TS15", "TS16") ~ "H2O2+PRG4",
    TRUE ~ "Other"
)
long_dat <- long_dat %>% filter(Condition != "Other")
long_dat$Condition <- factor(long_dat$Condition, levels = c("CTRL", "H2O2", "PRG4", "H2O2+PRG4"))

# Plot Barplots per Gene
# Select top 9 relevant ones for a facet plot? Or all? 
# CFH, C3, ARMS2, HTRA1, RPE65, BEST1...
# Let's plot all in a faceted way.
p_bar <- ggplot(long_dat, aes(x = Condition, y = TPM, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~hgnc_symbol, scales = "free_y") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none") +
  labs(title = "Known AMD Gene Expression", y = "TPM")

ggsave(file.path(out_dir, "Fig_Known_AMD_genes_barplot.pdf"), p_bar, width = 10, height = 8)

# Heatmap
# Calculate Z-scores
z_dat <- long_dat %>%
    group_by(!!sym(gene_col)) %>%
    mutate(Z = (TPM - mean(TPM)) / sd(TPM)) %>%
    ungroup()

p_heat <- ggplot(z_dat, aes(x = Sample, y = hgnc_symbol, fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#01665e", high = "#8c510a", mid = "white", midpoint = 0) +
    facet_grid(~Condition, scales = "free_x", space = "free_x") +
    theme_pub() +
    theme(axis.text.x = element_blank()) +
    labs(title = "Known AMD Gene Expression (Z-score)", x = "")

ggsave(file.path(out_dir, "Fig_Known_AMD_genes_heatmap.pdf"), p_heat, width = 8, height = 6)
