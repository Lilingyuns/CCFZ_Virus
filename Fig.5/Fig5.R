#### START ####
rm(list = ls())
library(ggplot2)
library(dplyr)
library(readr)

data <- read_tsv("Fig6.tsv")
data$COUNT <- as.numeric(as.character(data$COUNT))

####  Fig.5 a ####
p <- ggplot(data, aes(x = data$COUNT, y = data$Envrionment)) +
  geom_bar(stat = "identity", fill = "#4d7e54") +
  theme_minimal() +
  labs(
    x = "The probability that vOTU encodes at least one MRG (%)",
    y = NULL,
    title = NULL
  ) +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face = "plain", colour = "black"),
    axis.title.y = element_text(size = 12, face = "plain"),
    axis.text.x = element_text(size = 12, face = "plain", hjust = 0.5, colour = "black"),
    axis.title.x = element_text(size = 12, face = "plain"),
    legend.text = element_text(size = 10, face = "italic"),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(vjust = 0, hjust = 0.5, size = 14, lineheight = 1.4),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  guides(fill = guide_legend(ncol = 4))
