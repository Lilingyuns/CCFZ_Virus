#### START ####
rm(list = ls())
library(ggplot2)
library(dplyr)
library(broom)
library(patchwork)



#### Fig.6 c ####
data <- read.csv("Co_data_b_0mm.csv", sep = ",", header = TRUE)
data$Time <- factor(data$Time, levels = unique(data$Time))
stat_test <- data %>%
  group_by(Time) %>%
  do(tidy(t.test(OD600 ~ Strain, data = .))) %>%
  select(Time, p.value) %>%
  mutate(
    sig = case_when(
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


global_max_od <- max(data$OD600)  
y_pos_unified <- data %>%
  distinct(Time) %>%  
  mutate(y_pos = global_max_od * 0.93) %>%  
  left_join(stat_test, by = "Time") 


c_Co_0mm <- ggplot(data, aes(x = Time, y = OD600, color = Strain)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,  dodge.width = 0.8), 
             size = 1.55, alpha = 0.6 ) +  
  geom_vline(xintercept = seq(1.5, length(unique(data$Time)) - 0.5, 1),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.75, position = position_dodge(0.8)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.8),
               linewidth = 0.8) +
  geom_text(
    data = y_pos_unified,  
    aes(x = Time, y = y_pos, label = sig),
    inherit.aes = FALSE, 
    size = 5, 
    vjust = 0        
  ) +
  scale_color_manual(values = c("#EE9539", "#2BC380")) +
  labs(y = "OD600", color = "Strain") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank()
  )


data <- read.csv("Cd_data_b_0.5mm.csv", sep = ",", header = TRUE)
data$Time <- factor(data$Time, levels = unique(data$Time))
stat_test <- data %>%
  group_by(Time) %>%
  do(tidy(t.test(OD600 ~ Strain, data = .))) %>%
  select(Time, p.value) %>%
  mutate(
    sig = case_when(
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


global_max_od <- max(data$OD600)  
y_pos_unified <- data %>%
  distinct(Time) %>% 
  mutate(y_pos = global_max_od * 0.93) %>% 
  left_join(stat_test, by = "Time")  


c_Cd_0.5mm <- ggplot(data, aes(x = Time, y = OD600, color = Strain)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,  dodge.width = 0.8), 
             size = 1.55, alpha = 0.6 ) +
  geom_vline(xintercept = seq(1.5, length(unique(data$Time)) - 0.5, 1),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.75, position = position_dodge(0.8)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.8),
               linewidth = 0.8) +
  geom_text(
    data = y_pos_unified, 
    aes(x = Time, y = y_pos, label = sig),
    inherit.aes = FALSE, 
    size = 5, 
    vjust = 0        
  ) +
  scale_color_manual(values = c("#EE9539", "#2BC380")) +
  labs(y = "OD600", color = "Strain") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank()
  )



data <- read.csv("Cd_data_b_0.75mm.csv", sep = ",", header = TRUE)
data$Time <- factor(data$Time, levels = unique(data$Time))

stat_test <- data %>%
  group_by(Time) %>%
  do(tidy(t.test(OD600 ~ Strain, data = .))) %>%
  select(Time, p.value) %>%
  mutate(
    sig = case_when(
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

global_max_od <- max(data$OD600)  
y_pos_unified <- data %>%
  distinct(Time) %>%  
  mutate(y_pos = global_max_od * 0.88) %>%  
  left_join(stat_test, by = "Time")  


c_Cd_0.75mm <- ggplot(data, aes(x = Time, y = OD600, color = Strain)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,  dodge.width = 0.8), 
             size = 1.55, alpha = 0.6 ) +
  geom_vline(xintercept = seq(1.5, length(unique(data$Time)) - 0.5, 1),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.75, position = position_dodge(0.8)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.8),
               linewidth = 0.8) +
  geom_text(
    data = y_pos_unified,  
    aes(x = Time, y = y_pos, label = sig),
    inherit.aes = FALSE, 
    size = 5, 
    vjust = 0       
  ) +
  scale_color_manual(values = c("#EE9539", "#2BC380")) +
  
  labs(y = "OD600", color = "Strain") +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank()
  )



data <- read.csv("Co_data_b_1mm.csv", sep = ",", header = TRUE)
data$Time <- factor(data$Time, levels = unique(data$Time))


stat_test <- data %>%
  group_by(Time) %>%
  do(tidy(t.test(OD600 ~ Strain, data = .))) %>%
  select(Time, p.value) %>%
  mutate(
    sig = case_when(
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


global_max_od <- max(data$OD600)  
y_pos_unified <- data %>%
  distinct(Time) %>%  
  mutate(y_pos = global_max_od * 0.88) %>%  
  left_join(stat_test, by = "Time")  


c_Co_1mm <- ggplot(data, aes(x = Time, y = OD600, color = Strain)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,  dodge.width = 0.8), 
             size = 1.55, alpha = 0.6 ) +
  geom_vline(xintercept = seq(1.5, length(unique(data$Time)) - 0.5, 1),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.5, position = position_dodge(0.8)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.8),
               linewidth = 0.8) +
  geom_text(
    data = y_pos_unified,  
    aes(x = Time, y = y_pos, label = sig),
    inherit.aes = FALSE, 
    size = 5, 
    vjust = 0        
  ) +
  scale_color_manual(values = c("#EE9539", "#2BC380")) +
  labs(y = "OD600", color = "Strain") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank()
  )

#### Fig.6 d ####
data <- read.csv("Cu_data_b_0mm.csv", sep = ",", header = TRUE)
data$Time <- factor(data$Time, levels = unique(data$Time))
stat_test <- data %>%
  group_by(Time) %>%
  do(tidy(t.test(OD600 ~ Strain, data = .))) %>%
  select(Time, p.value) %>%
  mutate(
    sig = case_when(
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


global_max_od <- max(data$OD600)  
y_pos_unified <- data %>%
  distinct(Time) %>%  
  mutate(y_pos = global_max_od * 0.95) %>%  
  left_join(stat_test, by = "Time") 


y_Cu_0mm <- ggplot(data, aes(x = Time, y = OD600, color = Strain)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,  dodge.width = 0.8), 
             size = 1.55, alpha = 0.6 ) +
  geom_vline(xintercept = seq(1.5, length(unique(data$Time)) - 0.5, 1),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.75, position = position_dodge(0.8)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.8),
               linewidth = 0.8) +
  geom_text(
    data = y_pos_unified,  
    aes(x = Time, y = y_pos, label = sig),
    inherit.aes = FALSE, 
    size = 5, 
    vjust = 0        
  ) +
  scale_color_manual(values = c("#EE9539", "#127EE1")) +
  labs(y = "OD600", color = "Strain") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank()
  )


data <- read.csv("Cu_data_b_1mm.csv", sep = ",", header = TRUE)
data$Time <- factor(data$Time, levels = unique(data$Time))
stat_test <- data %>%
  group_by(Time) %>%
  do(tidy(t.test(OD600 ~ Strain, data = .))) %>%
  select(Time, p.value) %>%
  mutate(
    sig = case_when(
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )


global_max_od <- max(data$OD600)  
y_pos_unified <- data %>%
  distinct(Time) %>%  
  mutate(y_pos = global_max_od * 0.93) %>%  
  left_join(stat_test, by = "Time") 


y_Cu_1mm <- ggplot(data, aes(x = Time, y = OD600, color = Strain)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,  dodge.width = 0.8), 
             size = 1.55, alpha = 0.6 ) +   
  geom_vline(xintercept = seq(1.5, length(unique(data$Time)) - 0.5, 1),
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  stat_summary(fun = mean, geom = "point", size = 1.75, position = position_dodge(0.8)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(0.8),
               linewidth = 0.8) +
  geom_text(
    data = y_pos_unified,  
    aes(x = Time, y = y_pos, label = sig),
    inherit.aes = FALSE, 
    size = 5, 
    vjust = 0        
  ) +
  scale_color_manual(values = c("#EE9539", "#127EE1")) +
  labs(y = "OD600", color = "Strain") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank()
  )
