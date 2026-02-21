
library(ggplot2)

# Publication theme
theme_pub <- function() {
  theme_classic() +
    theme(
      text = element_text(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# Color palettes
axis_colors <- c(
  "Oxidative" = "#E69F00",  # Orange
  "Inflammatory" = "#56B4E9", # Sky Blue
  "Senescence" = "#009E73"   # Bluish Green
)
