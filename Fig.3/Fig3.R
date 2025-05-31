#### START ####
rm(list = ls())
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(vegan)
library(reshape2)
library(patchwork)
library(ggsci)
library(tidyr)

#### Fig.3 b ####
data <- read.csv("VHR_to_plot.csv")
colnames(data) <- c("Species", "BarValue1", "BarValue2", "LineValue1", "LineValue2")
data <- data %>%
  mutate(AverageLineValue = rowMeans(select(., LineValue1, LineValue2), na.rm = TRUE))

data <- data %>%
  arrange(AverageLineValue)

species_order <- as.character(data$Species)
data_long <- data %>% 
  pivot_longer(cols = c("BarValue1", "BarValue2"), names_to = "Sample", values_to = "Value") %>%
  na.omit()
scale_factor <- 70 / max(data$LineValue1, data$LineValue2, na.rm = TRUE)


Fig3b <- ggplot() +
  geom_bar(data = data_long, aes(x = Species, y = Value, fill = Sample), 
           stat = "identity", position = "dodge", alpha = 0.4) +
  geom_line(data = data, aes(x = Species, 
                             y = LineValue1 * scale_factor, group = 1), 
            color = "#FF69B4", linewidth = 1, na.rm = TRUE) +
  geom_point(data = data %>% filter(!is.na(LineValue1)), 
             aes(x = Species, 
                 y = LineValue1 * scale_factor), 
             color = "#FF69B4", size = 3) +
  geom_line(data = data, aes(x = Species, 
                             y = LineValue2 * scale_factor, group = 1), 
            color = "#5baeff", linewidth = 1, na.rm = TRUE) +
  geom_point(data = data %>% filter(!is.na(LineValue2)), 
             aes(x = Species, 
                 y = LineValue2 * scale_factor), 
             color = "#5baeff", size = 3) +
  scale_x_discrete(limits = species_order) +
  scale_y_continuous(name = "Host relative abundance (%)", 
                     limits = c(0, 70),
                     sec.axis = sec_axis(~ . / scale_factor, name = "Log10(VHR)"),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("#FF69B4", "#5baeff"), labels = c("Sample 1", "Sample 2")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(hjust = 0.5, size = 14, colour = 'black', angle = 0),
    axis.text.y = element_text(size = 14, colour = 'black'),
    legend.position = 'none'
  ) +
  coord_flip() +
  labs(x = "")


#### Fig.3 c-f ####
abundance <- read.csv("VHR_to_plot_lifestyle_0219.csv")
split_data <- split(abundance, abundance$Lifestyle)


df1 <- split_data[[1]]
df2 <- split_data[[2]]
df3 <- split_data[[3]]



result_Fig3c <- cor.test(df1$Host_abundance_log, df1$Virus_abundance_log, method = "spearman")
result_Fig3d <- cor.test(df1$Host_abundance_log, df1$VHR, method = "spearman")

result_Fig3e <- cor.test(df2$Host_abundance_log, df2$Virus_abundance_log, method = "spearman")
result_Fig3f <- cor.test(df2$Host_abundance_log, df2$VHR, method = "spearman")


mycol <- c("red", "#faad61", "#71abc8")


Fig3c <- ggplot(df1, aes(x = Host_abundance_log, y = Virus_abundance_log, color = Lifestyle, fill = Lifestyle)) +
  geom_point(color = "grey", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", aes(fill = Lifestyle), se = FALSE, size = 1.6) +
  scale_color_manual(values = "red") +
  scale_fill_manual(values = "red") +
  labs(x = "Log10(HostMAG abundance)", y = "Log10(Viral abundance)") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.x = element_text(size = 14, face = "plain"), 
        axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.y = element_text(size = 14, face = "plain"), 
        legend.text = element_text(size = 14, face = "plain"),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(vjust = 0, hjust = 0.5, size = 14, lineheight = 1.4)) +
  theme(legend.position = "none")


Fig3d <- ggplot(df1, aes(x = Host_abundance_log, y = VHR, color = Lifestyle, fill = Lifestyle)) +
  geom_point(color = "grey", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", aes(fill = Lifestyle), se = FALSE, size = 1.6) +
  scale_color_manual(values = "red") +
  scale_fill_manual(values =  "red") +
  labs(x = "Log10(HostMAG abundance)", y = "Log10(VHR)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.x = element_text(size = 14, face = "plain", ), 
        axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.y = element_text(size = 14, face = "plain", ), 
        legend.text = element_text(size = 14, face = "plain"),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text (vjust = 0, hjust = 0.5, size = 14, lineheight = 1.4))+
  theme(legend.position="none")

Fig3e <- ggplot(df2, aes(x = Host_abundance_log, y = Virus_abundance_log, color = Lifestyle, fill = Lifestyle)) +
  geom_point(color = "grey", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", aes(fill = Lifestyle), se = FALSE, size = 1.6) +
  scale_color_manual(values = "#faad61") +
  scale_fill_manual(values = "#faad61") +
  labs(x = "Log10(HostMAG abundance)", y = "Log10(Viral abundance)") +
  scale_y_continuous(breaks = seq(floor(min(df2$Virus_abundance_log)), ceiling(max(df2$Virus_abundance_log)), by = 1)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.x = element_text(size = 14, face = "plain"), 
        axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.y = element_text(size = 14, face = "plain"), 
        legend.text = element_text(size = 14, face = "plain"),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(vjust = 0, hjust = 0.5, size = 14, lineheight = 1.4)) +
  theme(legend.position = "none")


Fig3f <- ggplot(df2, aes(x = Host_abundance_log, y = VHR, color = Lifestyle, fill = Lifestyle)) +
  geom_point(color = "grey", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", aes(fill = Lifestyle), se = FALSE, size = 1.6) +
  scale_color_manual(values = "#faad61") +
  scale_fill_manual(values =  "#faad61") +
  labs(x = "Log10(HostMAG abundance)", y = "Log10(VHR)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.x = element_text(size = 14, face = "plain", ), 
        axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.y = element_text(size = 14, face = "plain", ), 
        legend.text = element_text(size = 14, face = "plain"),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text (vjust = 0, hjust = 0.5, size = 14, lineheight = 1.4))+
  theme(legend.position="none")

