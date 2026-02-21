
source("code/three-axes/utils/plotting_theme.R")
library(ggplot2)

# Coordinates for boxes
boxes <- data.frame(
  x = c(0, -1, 0, 1, 0, 0),
  y = c(3, 2, 2, 2, 0, 4), # Sen, Inf, Ox on level 2?
  # Let's organize hierarchy
  # Level 3: H2O2 (Top)
  # Level 2: Axes (Ox, Inf, Sen)
  # Level 1: AMD (Bottom)
  # Level 2.5: PRG4 (Side/Inhibition)
  label = c("Oxidative", "Inflammatory", "Senescence", "AMD Pathobiology", "H2O2 Treatment", "PRG4 Rescue")
)

# Refined positions
nodes <- data.frame(
  id = c("H2O2", "Ox", "Inf", "Sen", "AMD", "PRG4"),
  x = c(0, -1, 0, 1, 0, -2),
  y = c(3, 1.5, 1.5, 1.5, 0, 1.5),
  label = c("H2O2", "Oxidative\nAxis", "Inflammatory\nAxis", "Senescence\nAxis", "AMD\nPathology", "PRG4"),
  fill = c("#ffcccb", "#E69F00", "#56B4E9", "#009E73", "#E6E6FA", "#90EE90"),
  type = c("Trigger", "Axis", "Axis", "Axis", "Outcome", "Treatment")
)

edges <- data.frame(
  from_x = c(0, 0, 0, -1, 0, 1, -2, -2, -2),
  from_y = c(2.7, 2.7, 2.7, 1.2, 1.2, 1.2, 1.5, 1.5, 1.5),
  to_x = c(-1, 0, 1, 0, 0, 0, -1.3, -0.3, 0.7), # Approximate for PRG4 inhibition
  to_y = c(1.8, 1.8, 1.8, 0.3, 0.3, 0.3, 1.5, 1.5, 1.5),
  type = c("Induction", "Induction", "Induction", "Promotion", "Promotion", "Promotion", "Inhibition", "Inhibition", "Inhibition")
)

p <- ggplot() +
  # Edges
  geom_segment(data = edges[edges$type!="Inhibition",], aes(x = from_x, y = from_y, xend = to_x, yend = to_y), 
               arrow = arrow(length = unit(0.3, "cm")), size = 1) +
  # Inhibition edges (tee) - simulate with segment + cross bar? Or just red color
  geom_segment(data = edges[edges$type=="Inhibition",], aes(x = from_x, y = from_y, xend = to_x, yend = to_y), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), color = "red", linetype = "dashed", size = 1.2) +
   
  # Nodes
  geom_rect(data = nodes, aes(xmin = x - 0.4, xmax = x + 0.4, ymin = y - 0.3, ymax = y + 0.3, fill = fill), color = "black") +
  geom_text(data = nodes[nodes$id!="PRG4",], aes(x = x, y = y, label = label), size = 4) +
  geom_text(data = nodes[nodes$id=="PRG4",], aes(x = x, y = y, label = label), size = 4, fontface="bold") +
  
  scale_fill_identity() +
  theme_void() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  labs(title = "Mechanism of PRG4 Rescue in AMD") +
  xlim(-2.5, 1.5) + ylim(-0.5, 3.5)

ggsave("results/three-axes/summary/Fig_Mechanism_Diagram.pdf", p, width = 8, height = 6)
