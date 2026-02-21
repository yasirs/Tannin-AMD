
source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# Constants
out_dir <- "results/three-axes/senescence-deep-dive"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load SenMayo Genes
senmayo_mapped <- read_csv("results/three-axes/gene-sets/senescence_pro_disease_mapped.csv", show_col_types = FALSE)
senmayo_genes <- senmayo_mapped$gene_symbol

# Load DE Results
de_res <- read_csv("data/RPE_DE_results.csv", show_col_types = FALSE)

# Filter
sen_de <- de_res %>% filter(hgnc_symbol %in% senmayo_genes)

# 1. Enrichment Plot (Volcano style or Boxplot)
# Compare SenMayo genes vs Background
# H2O2 Induction
de_res$IsSenMayo <- de_res$hgnc_symbol %in% senmayo_genes

p1 <- ggplot(de_res, aes(x = IsSenMayo, y = H2O2_vs_CTRL.log2FoldChange, fill = IsSenMayo)) +
  geom_boxplot() +
  theme_pub() +
  labs(title = "SenMayo Genes Induction (H2O2)", x = "Is SenMayo?", y = "Log2 Fold Change") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + 
  scale_fill_manual(values = c("grey", "#009E73")) # Green for senescence

ggsave(file.path(out_dir, "Fig_SenMayo_Induction_Boxplot.pdf"), p1, width = 4, height = 5)

# PRG4 Restoration
p2 <- ggplot(de_res, aes(x = IsSenMayo, y = H2O2PRG4_vs_H2O2.log2FoldChange, fill = IsSenMayo)) +
  geom_boxplot() +
  theme_pub() +
  labs(title = "SenMayo Genes Rescue (PRG4)", x = "Is SenMayo?", y = "Log2 Fold Change (PRG4 Effect)") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) + 
  scale_fill_manual(values = c("grey", "#009E73"))

ggsave(file.path(out_dir, "Fig_SenMayo_Rescue_Boxplot.pdf"), p2, width = 4, height = 5)

# 2. Heatmap of Top SenMayo Genes
# Top induced
top_sen <- sen_de %>% arrange(desc(H2O2_vs_CTRL.log2FoldChange)) %>% head(30)
top_genes <- top_sen$hgnc_symbol

# Load TPMs
tpms <- read_csv("data/RPE_TPMS.csv", show_col_types = FALSE)
# Map samples
# Need to pivot
gene_col <- "hgnc_symbol"
tpm_subset <- tpms %>% select(-any_of("ensembl_gene_id"))

long_tpms <- tpm_subset %>%
  filter(!!sym(gene_col) %in% top_genes) %>%
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

# Z-score
long_tpms <- long_tpms %>%
    group_by(!!sym(gene_col)) %>%
    mutate(Z = (TPM - mean(TPM)) / sd(TPM)) %>%
    ungroup()

p3 <- ggplot(long_tpms, aes(x = Sample, y = hgnc_symbol, fill = Z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#01665e", high = "#8c510a", mid = "white", midpoint = 0) +
    facet_grid(~Condition, scales = "free_x", space = "free_x") +
    theme_pub() +
    theme(axis.text.x = element_blank()) +
    labs(title = "Top SenMayo Gene Expression", x = "")

ggsave(file.path(out_dir, "Fig_SenMayo_Heatmap.pdf"), p3, width = 6, height = 8)
