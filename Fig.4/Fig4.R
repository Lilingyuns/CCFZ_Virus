#### START ####
rm(list = ls())
library(tidyverse)
library(vegan)
library(vegan)
library(ggplot2)
library(ggalt)
library(dplyr)
library(reshape2)
library(ggpubr)

#### Fig.5 b ####
data <- read.table("03AMG_abundance_KO_merged_to_draw.csv", header = TRUE, sep = ",", row.names = 1)
group <- read.delim("metadata.txt", header = TRUE, sep = "\t", row.names = 1)

data <- decostand(data, method='hellinger')
mycol <- c("#5baeff", "#FF69B4")
distance <- vegdist(data, method = 'bray')
nmds <- metaMDS(distance, k = 2)
stress <- nmds$stress
df <- as.data.frame(nmds$points)


df <- cbind(df, group)


anosim = anosim(data, group$Type1, permutations = 999, distance = "bray")
print(anosim)
anosim_text <- paste(paste("R = ", round(anosim$statistic, 2)), "  p = 0.001")
stress_text <- paste("Stress  = ", round(stress, 4),"\n")


Fig5b <- ggplot(df, aes(MDS1, MDS2)) +
  geom_point(aes(color = Type1), size = 4, alpha = 0.7) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  theme_bw()+
  stat_ellipse(data = df, aes(color = Type1), level=0.95, linetype="dashed", linewidth=1.5, show.legend = FALSE) +
ggtitle(paste(stress_text, anosim_text))+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.x = element_text(size = 14, face = "plain", ), 
        axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.y = element_text(size = 14, face = "plain", ), 
        legend.text = element_text(size = 14, face = "plain"),
        legend.title = element_blank(),
        legend.position = 'bottom',
        plot.title = element_text (vjust = -20, hjust = 0.1, size = 14, lineheight = 1.4))

#### Fig.5 c ####
data_A <- read.table("03AMG_abundance_freq_to_draw.csv", header = TRUE, sep = ",")
data_A <- data_A[1:49, ]  


metadata <- read.delim("metadata.txt", header = TRUE, sep = "\t")


data_long <- data_A %>%
  pivot_longer(cols = -Sample, names_to = "Pathway", values_to = "value") %>%
  mutate(Sample = as.character(Sample), Pathway = as.character(Pathway))


data_B <- data_long %>%
  left_join(metadata, by = "Sample") %>%
  select(labels = Pathway, variable = Type1, value)


data_B <- as_tibble(data_B)


dif <- data_B %>%
  mutate(labels = gsub('\\.', ' ', labels))


dif_summary <- dif %>%
  group_by(labels, variable) %>%
  summarize(mean_value = mean(value), .groups = "drop")


Fig5c <- ggplot(dif_summary, aes(x = labels, y = ifelse(variable == "Sediment", mean_value, -mean_value), fill = variable)) +
  geom_bar(stat = "identity") +  
  geom_text(
    aes(label = signif(mean_value, 2), vjust = ifelse(variable == "Sediment", -0.6, 1.3)), 
    size = 4, 
    family = "sans" 
  ) +
  scale_fill_manual(
    name = "Sample Type", 
    values = c("#FF69B4", "#5baeff"), 
    breaks = c("Sediment", "Nodule"), 
    labels = c("Sediment", "Nodule")  
  ) +
  scale_y_continuous(
    limits = c(-2.5, 2.5),  
    labels = abs,           
    expand = expansion(mult = c(0.1, 0.1))  
  ) +
  labs(x = "", y = "Viral frequency (%)", title = "") + 
 
  stat_compare_means(
    data = dif, 
    aes(x = labels, y = value, group = variable), 
    method = "wilcox.test",      
    label = "p.signif",         
    vjust = 20, size = 5,        
    hide.ns = TRUE               
  ) +
  theme_bw() +  
  theme(
    legend.position = "none",  
    panel.grid = element_blank(),  
    axis.text.x = element_text(size = 12, face = "plain", hjust = 1, colour = "black", angle = 60), 
    axis.title.x = element_blank(),  
    axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"), 
    axis.title.y = element_text(size = 14, face = "plain"),  
    text = element_text(size = 14, family = "sans") 
  )

