
source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# Constants
out_dir <- "results/three-axes/senescence-deep-dive"

# Markers
# Senescence: CDKN1A (p21), CDKN2A (p16), LMNB1 (Loss), GLB1 (Beta-gal)
# SASP: IL6, CXCL8, MMP1, MMP3, SERPINE1 (PAI-1)
markers <- c("CDKN1A", "CDKN2A", "LMNB1", "GLB1", 
             "IL6", "CXCL8", "MMP1", "MMP3", "SERPINE1")

# Load TPMs
tpm_file <- "data/RPE_TPMS.csv"
tpms <- read_csv(tpm_file, show_col_types = FALSE)

# Filter
gene_col <- "hgnc_symbol"
tpm_subset <- tpms %>% select(-any_of("ensembl_gene_id"))
axis_tpm <- tpm_subset %>% filter(!!sym(gene_col) %in% markers)

# Long format
long_tpms <- axis_tpm %>%
  pivot_longer(cols = -!!sym(gene_col), names_to = "Sample", values_to = "TPM")

long_tpms$Condition <- case_when(
    long_tpms$Sample %in% c("TS1", "TS2", "TS3") ~ "CTRL",
    long_tpms$Sample %in% c("TS5", "TS6", "TS7") ~ "H2O2",
    long_tpms$Sample %in% c("TS9", "TS11", "TS12") ~ "PRG4",
    long_tpms$Sample %in% c("TS13", "TS15", "TS16") ~ "H2O2+PRG4",
    TRUE ~ "Other"
)
long_tpms <- long_tpms %>% filter(Condition != "Other")
long_tpms$Condition <- factor(long_tpms$Condition, levels = c("CTRL", "H2O2", "PRG4", "H2O2+PRG4"))

# Plot Barplots
p <- ggplot(long_tpms, aes(x = Condition, y = TPM, fill = Condition)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  facet_wrap(~hgnc_symbol, scales = "free_y") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position = "none") +
  labs(title = "Senescence & SASP Markers", y = "TPM")

ggsave(file.path(out_dir, "Fig_Senescence_Markers_Barplot.pdf"), p, width = 10, height = 8)
