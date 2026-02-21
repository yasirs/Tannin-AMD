library(ggplot2)
library(dplyr)
library(tidyr)

comparison_data <- read.csv("results/cohort-GSE29801/dry_amd_de/comprehensive_comparison.csv")

comparison_long <- comparison_data %>%
  pivot_longer(cols = c(Features_FDR05, Features_P001),
               names_to = "Threshold",
               values_to = "Count") %>%
  filter(!is.na(Count))

comparison_long$Threshold <- ifelse(comparison_long$Threshold == "Features_FDR05",
                                   "FDR < 0.05", "p < 0.01")

method_order <- c(
  "Original DE (individual stages)",
  "Improved DE: Combined Early",
  "Improved DE: Combined Intermediate", 
  "Improved DE: All Dry AMD",
  "Differential Variability (F-test)",
  "Disease Variance (Age-Adjusted)",
  "Age Variance",
  "Age Main Effect",
  "Age × Disease Interaction"
)

comparison_long$Method <- factor(comparison_long$Method, levels = method_order)

p <- ggplot(comparison_long, aes(x = Method, y = Count, fill = Threshold)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Tissue, scales = "free_y", ncol = 1) +
  labs(
    title = "Comprehensive Signal Detection Comparison",
    x = "Analysis Method",
    y = "Number of Genes Detected",
    fill = "Significance"
  ) +
  theme_classic(base_family = "serif") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  ) +
  scale_fill_manual(values = c("FDR < 0.05" = "#1f78b4", "p < 0.01" = "#e31a1c"))

ggsave("results/cohort-GSE29801/dry_amd_de/comprehensive_comparison_plot.tiff", 
       p, width = 12, height = 10, compression = "lzw", dpi = 300)

message("✓ Saved: comprehensive_comparison_plot.tiff")
