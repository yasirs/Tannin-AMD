
# 16_summary_figure.R
# Purpose: Generate Summary Bubble Plot using GSEA NES from all phases.

source("code/three-axes/utils/plotting_theme.R")
library(dplyr)
library(ggplot2)
library(readr)

out_dir <- "results/three-axes/summary"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Load NES Results from various steps

# 1. H2O2 Induction & PRG4 Rescue (From 06_axis_scores.R output)
# Note: 06 now outputs "axis_nes_summary.csv" which has H2O2_vs_CTRL and H2O2PRG4_vs_H2O2
nes_file <- "results/three-axes/axis-scores/axis_nes_summary.csv"

if (file.exists(nes_file)) {
    nes_dat <- read_csv(nes_file, show_col_types = FALSE)
    
    # Extract Induction
    ind_res <- nes_dat %>% 
        filter(Contrast == "H2O2_vs_CTRL") %>%
        mutate(Type = "H2O2 Induction", value = NES, pval = padj, Axis = pathway) %>%
        select(Axis, Type, value, pval)
        
    # Extract Rescue
    rec_res <- nes_dat %>% 
        filter(Contrast == "H2O2PRG4_vs_H2O2") %>%
        mutate(Type = "PRG4 Rescue", value = NES, pval = padj, Axis = pathway) %>%
        select(Axis, Type, value, pval)
} else {
    ind_res <- data.frame()
    rec_res <- data.frame()
}

# 2. AMD Dysregulation (From 08_amd_cohort_validation.R)
amd_file <- "results/three-axes/amd-validation/amd_axis_dysregulation_gsea.csv"
if (file.exists(amd_file)) {
    amd_dat <- read_csv(amd_file, show_col_types = FALSE)
    # pathway col has axis name
    # Clean axis names if needed (e.g. "oxidative_pro_disease" -> "Oxidative")
    # Actually the script 08 passed lists like "Oxidative" = ... so names should be clean.
    
    val_res <- amd_dat %>%
        mutate(Type = "AMD Dysregulation", value = NES, pval = padj, Axis = pathway) %>%
        select(Axis, Type, value, pval)
} else {
    val_res <- data.frame()
}

# Combine
all_res <- rbind(ind_res, rec_res, val_res)

# Plot Bubble Plot
# Size = -log10(FDR), Color = NES
all_res$LogP <- -log10(all_res$pval)
all_res$LogP[is.infinite(all_res$LogP)] <- 5 # Cap infinite
all_res$LogP[is.na(all_res$LogP)] <- 0

# Ensure Factor Order for X-axis
all_res$Type <- factor(all_res$Type, levels = c("H2O2 Induction", "AMD Dysregulation", "PRG4 Rescue"))

p <- ggplot(all_res, aes(x = Type, y = Axis, size = LogP, color = value)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", name = "NES") +
  theme_pub() +
  labs(title = "Three Axes Summary (GSEA NES)", x = "", y = "", size = "-log10(FDR)") +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_size(range = c(3, 10))

ggsave(file.path(out_dir, "Fig_Summary_Bubble.pdf"), p, width = 6, height = 4)
