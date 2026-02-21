#%%
# Resize Age GSEA Summary Plot
# Loads existing GSEA CSV results and replots with increased width (24 inches)

library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)

# Configuration
base_dir <- "/home/ysuhail/work/Tannin-AMD"
results_dir <- file.path(base_dir, "results/cohort-GSE135092/rpe_covariate_de")
gsea_dir <- file.path(results_dir, "gsea")
figs_dir <- file.path(results_dir, "figures")

# Load Results
cat("Loading GSEA results from CSV...\n")
gsea_mac_h <- read_csv(file.path(gsea_dir, "Age_GSEA_Hallmark_macula.csv"), show_col_types = FALSE)
gsea_mac_c2 <- read_csv(file.path(gsea_dir, "Age_GSEA_C2_macula.csv"), show_col_types = FALSE)
gsea_nonmac_h <- read_csv(file.path(gsea_dir, "Age_GSEA_Hallmark_nonmacula.csv"), show_col_types = FALSE)
gsea_nonmac_c2 <- read_csv(file.path(gsea_dir, "Age_GSEA_C2_nonmacula.csv"), show_col_types = FALSE)

# Plotting Function
plot_gsea <- function(df, title) {
  if(nrow(df) == 0) return(NULL)
  
  # Filter for positive NES (UP with Age)
  df_up <- df %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(10)
  
  if(nrow(df_up) == 0) return(NULL)
  
  p <- ggplot(df_up, aes(x = reorder(pathway, NES), y = NES, fill = padj)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(low = "red", high = "blue", name = "FDR") +
    labs(title = title, x = "", y = "NES") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    )
  return(p)
}

# Create Plots
p1 <- plot_gsea(gsea_mac_h, "Age UP - Macula (Hallmark)")
p2 <- plot_gsea(gsea_mac_c2, "Age UP - Macula (C2)")
p3 <- plot_gsea(gsea_nonmac_h, "Age UP - non-Macula (Hallmark)")
p4 <- plot_gsea(gsea_nonmac_c2, "Age UP - non-Macula (C2)")

# Combine and Save
plist <- list(p1, p2, p3, p4)
plist <- plist[!sapply(plist, is.null)]

cat(sprintf("Generated %d plots.\n", length(plist)))

if(length(plist) > 0) {
  p_combined <- grid.arrange(grobs = plist, ncol = 2)
  
  # Use width = 24 (approximates 50% increase from 16)
  cat("Saving as PDF and TIFF (width=24)...\n")
  ggsave(file.path(figs_dir, "Age_GSEA_summary.pdf"), p_combined, width = 24, height = 12)
  ggsave(file.path(figs_dir, "Age_GSEA_summary.tiff"), p_combined, 
         width = 24, height = 12, dpi = 300, compression = "lzw")
  cat("Done.\n")
} else {
  cat("Nothing to plot.\n")
}
