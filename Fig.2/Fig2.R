#### START ####
rm(list = ls())
library(vegan)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork) 
library(ggsci) 
library(cowplot)
library(dplyr)  
library(gridExtra) 
library(tidyr)

#### Fig.2 a b####
data <- t(read.table("vOTU_aligned_90_abundance.csv", header = TRUE, sep = ",", row.names = 1))
group <- read.delim("metadata.txt", header = TRUE, sep = "\t")


richness <- rowSums(data > 0)
shannon <- diversity(data, index = "shannon")
pielou <- shannon / log(richness, exp(1))
simpson<- diversity(data, index = 'simpson')


data_richness = data.frame(richness)
data_shannon = data.frame(shannon)
data_pielou = data.frame(pielou)
data_simpson = data.frame(simpson)

spe_alpha <- cbind(
  data_richness,
  data_shannon,
  data_pielou,
  data_simpson,
  group
)


mycol <- c("#5baeff", "#FF69B4")

Fig2a <- ggplot(spe_alpha, aes(x = Type1, y = shannon, color = Type1)) +
  geom_boxplot(aes(color = Type1), lwd = 0.8, alpha = 0.4) + 
  geom_jitter(aes(color = Type1), 
              position = position_jitterdodge(jitter.height = 0, 
                                              jitter.width = 0.2, 
                                              dodge.width = 0.0),
              alpha = 0.4) + 
  labs(y = "Shannon index (Virus)") + 
  scale_color_manual(values = mycol)+
  scale_fill_manual(values = mycol)+
  theme_bw() + 
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     vjust = 1 , hjust = 0.5, size = 5) + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.y = element_text(size = 14, face = "plain", )) + 
  xlab(NULL) +
  theme(legend.position="none")

test_result <- wilcox.test(richness ~ Type1, data = spe_alpha)

Fig2b <- ggplot(spe_alpha, aes(x = Type1, y = richness, color = Type1)) +
  geom_boxplot(aes(color = Type1), lwd = 0.8, alpha = 0.4) + 
  geom_jitter(aes(color = Type1), 
              position = position_jitterdodge(jitter.height = 0, 
                                              jitter.width = 0.2, 
                                              dodge.width = 0.0),
              alpha = 0.4) + 
  labs(y = "Richness index (Virus)") + 
  scale_color_manual(values = mycol)+
  scale_fill_manual(values = mycol)+
  theme_bw() + 
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     vjust = 1 , hjust = 0.5, size = 5) + 
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.y = element_text(size = 14, face = "plain", )) + 
  xlab(NULL) +
  theme(legend.position="none")

test_result <- wilcox.test(shannon ~ Type1, data = spe_alpha)

#### Fig.2 c ####
data1 <- read.table("vOTU_aligned_90_abundance_alpha_alpha.csv", header = TRUE, sep = ",", row.names = 1)
group1 <- read.delim("metadata_alpha_alpha.txt", header = TRUE, sep = "\t")
data2 <- read.table("MAG_aligned_90_abundance.csv", header = TRUE, sep = ",", row.names = 1)
group2 <- read.delim("metadata_MAG.txt", header = TRUE, sep = "\t")


richness1 <- rowSums(data1 > 0)
shannon1 <- diversity(data1, index = "shannon")
pielou1 <- shannon1 / log(richness1, exp(1))
simpson1<- diversity(data1, index = 'simpson')

richness2 <- rowSums(data2 > 0)
shannon2 <- diversity(data2, index = "shannon")
pielou2 <- shannon1 / log(richness2, exp(1))
simpson2<- diversity(data2, index = 'simpson')

data_richness1 = data.frame(richness1)
data_shannon1 = data.frame(shannon1)
data_pielou1 = data.frame(pielou1)
data_simpson1 = data.frame(simpson1)

spe_alpha1 <- cbind(
  data_richness1,
  data_shannon1,
  data_pielou1,
  data_simpson1,
  group1
)

data_richness2 = data.frame(richness2)
data_shannon2 = data.frame(shannon2)
data_pielou2 = data.frame(pielou2)
data_simpson2 = data.frame(simpson2)

spe_alpha2 <- cbind(
  data_richness2,
  data_shannon2,
  data_pielou2,
  data_simpson2,
  group2
)

result_p1 <- cor.test(spe_alpha1$shannon1, spe_alpha2$shannon2, method = "spearman")

df1 <- data.frame(shannon1 = spe_alpha1$shannon1, shannon2 = spe_alpha2$shannon2)

Fig2c <- ggplot(df1, aes(x = shannon2, y = shannon1)) +
  geom_point(color = "grey", size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", color = "black", fill = "black", level = 0.95, se = TRUE, size = 1.6, alpha = 0.2) +
  labs(x = "Shannon index (Prokaryotes)", y = "Shannon index (Virus)") +
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

#### Fig.2 d e ####
data <- read.table("MAG_aligned_90_abundance_DEPTH.csv", header = TRUE, sep = ",", row.names = 1)
group <- read.delim("metadata_MAG_Depth.txt", header = TRUE, sep = "\t", row.names = 1)
vdata <- read.table("vOTU_abundance_t_Depth.csv", header = TRUE, sep = ",", row.names = 1)
vgroup <- read.delim("metadata_vOTU_Depth.txt", header = TRUE, sep = "\t", row.names = 1)


richness <- rowSums(data > 0)
shannon <- diversity(data, index = "shannon")
pielou <- shannon / log(richness, exp(1))
simpson<- diversity(data, index = 'simpson')

vrichness <- rowSums(vdata > 0)
vshannon <- diversity(vdata, index = "shannon")
vpielou <- vshannon / log(vrichness, exp(1))
vsimpson<- diversity(vdata, index = 'simpson')


data_richness = data.frame(richness)
data_shannon = data.frame(shannon)
data_pielou = data.frame(pielou)
data_simpson = data.frame(simpson)

vdata_richness = data.frame(vrichness)
vdata_shannon = data.frame(vshannon)
vdata_pielou = data.frame(vpielou)
vdata_simpson = data.frame(vsimpson)

spe_alpha <- cbind(
  data_richness,
  data_shannon,
  data_pielou,
  data_simpson,
  group
)


vspe_alpha <- cbind(
  vdata_richness,
  vdata_shannon,
  vdata_pielou,
  vdata_simpson,
  vgroup
)

colnames(vspe_alpha) <- colnames(spe_alpha)
spe_alpha$Group <- "Microbiome(MAG)"
vspe_alpha$Group <- "Virus(vOTU)"

result_Fig2_d1 <- cor.test(spe_alpha$shannon, spe_alpha$Depth, method = "spearman")
result_Fig2_d1
result_Fig2_d2 <- cor.test(vspe_alpha$shannon, vspe_alpha$Depth, method = "spearman")
result_Fig2_d2

result_Fig2_e1 <- cor.test(spe_alpha$richness, spe_alpha$Depth, method = "spearman")
result_Fig2_e1
result_Fig2_e2 <- cor.test(vspe_alpha$richness, vspe_alpha$Depth, method = "spearman")
result_Fig2_e2


Spe_alpha <- merge(spe_alpha, vspe_alpha, by = intersect(names(spe_alpha), names(vspe_alpha)), all = TRUE)

mycol <- c("#7a70b5", "#81b3a9")

Fig2d <- ggplot(Spe_alpha, aes(x = Depth, y = shannon, color = Group, fill = Group)) +
  geom_point(aes(color = Group), size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", aes(fill = Group), level = 0.95, se = TRUE, size = 1.6) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  labs(x = "Depth (cmbsf)", y = "Shannon index")+
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

Fig2e <- ggplot(Spe_alpha, aes(x = Depth, y = richness, color = Group, fill = Group)) +
  geom_point(aes(color = Group), size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", aes(fill = Group), level = 0.95, se = TRUE, size = 1.6) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  labs(x = "Depth (cmbsf)", y = "Richness index")+
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


#### Fig.2 f ####
data <- t(read.table("vOTU_aligned_90_abundance.csv", header = TRUE, sep = ",", row.names = 1))
group <- read.delim("metadata.txt", header = TRUE, sep = "\t")
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

Fig2f <- ggplot(df, aes(MDS1, MDS2)) +
  geom_point(aes(color = Type1), size = 2, alpha = 0.7) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  theme_bw()+
  stat_ellipse(data = df, aes(color = Type1), level=0.95, linetype="dashed", linewidth=1.5, show.legend = FALSE)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.x = element_text(size = 14, face = "plain", ), 
        axis.text.y = element_text(size = 14, face = "plain", hjust = 0.5, colour = "black"),
        axis.title.y = element_text(size = 14, face = "plain", ), 
        legend.text = element_text(size = 14, face = "plain"),
        legend.title = element_blank())+
  theme(legend.position="none")

#### Fig.2 g h ####
species_taxo <- read.csv("vOTU_taxonomy_a.csv")
species_taxo_long <- species_taxo %>%
  pivot_longer(cols = -Name, names_to = "Viral_Family", values_to = "Abundance")

sample_order <- c("All", "Sediment", "Nodule")
species_taxo_long$Name <- factor(species_taxo_long$Name, levels = sample_order)

viral_family_order <- c("Unclassified","Classified")

species_taxo_long$Viral_Family <- factor(species_taxo_long$Viral_Family, levels = viral_family_order)

num_families <- length(unique(species_taxo_long$Viral_Family))
mycol <- c("#DCDCDC","#1B9E77")


Fig2g <- ggplot(data = species_taxo_long, aes(x = Name, y = Abundance, fill = Viral_Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = mycol) +
  theme_minimal() +
  labs(x = "Sample",
       y = "Percentage of classified vOTUs (%)",
       fill = "Viral Family") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain", colour = "black"),  
    axis.title.x = element_text(size = 12, face = "plain"), 
    axis.text.y = element_text(size = 10, face = "plain", hjust = 0.5, colour = "black"),
    axis.title.y = element_text(size = 12, face = "plain"), 
    legend.text = element_text(size = 10, face = "plain"),
    legend.title = element_blank(),
    plot.title = element_text(vjust = 0, hjust = 0.5, size = 14, lineheight = 1.4),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid = element_blank(),  
    axis.ticks.y = element_line(color="black", size=0.8, lineend = 0.8),  
    axis.ticks.x = element_line(color="black", size=0.8, lineend = 0.8)   
  ) +
  theme(legend.position = "none" ) +
  xlab(NULL) +
  coord_cartesian(ylim = c(-0.2, 100.2)) + 
  scale_y_continuous(expand = c(0, 0))  


species_abundance <- read.csv("abundance_to_draw_merged.csv")


species_abundance_long <- species_abundance %>%
  pivot_longer(cols = -X, names_to = "Viral_Family", values_to = "Abundance")

species_abundance_long$Viral_Family <- gsub("\\.\\d+$", "", species_abundance_long$Viral_Family)

species_abundance_long$Viral_Family <- gsub("\\.", " ", species_abundance_long$Viral_Family)
species_abundance_long$Abundance[is.na(species_abundance_long$Abundance)] <- 0


viral_family_order2 <- c("unclassified",	
                         "Ackermannviridae",	"Autographiviridae", "Autolykiviridae",	
                         "Demerecviridae",	"Drexlerviridae",	"Fuselloviridae",	"Herelleviridae",	
                         "Inoviridae",	"Kyanoviridae",	"Lavidaviridae",	"Mimiviridae",	"Phycodnaviridae",	
                         "Rountreeviridae",	"Straboviridae",	"Zobellviridae","Unc  Caudovirales at family level")
mycol2 <- c("#DCDCDC", "#47A98C", "#D57D3A", "#796FAA", "#C33F96", "#A7CB5D",
                     "#CBCDFA", "#86AECE", "#E08A86", "blue", "#d62728", "#F0C1CA",
                     "#FFFFB3", "#bcbd22", "#c49c94", "#BC80BD", "#56B4E9")
                     

species_abundance_long$Viral_Family <- factor(species_abundance_long$Viral_Family, levels = viral_family_order2)
species_abundance_long$X <- factor(species_abundance_long$X, levels = c("All", "Sediment", "Nodule"))


Fig2h <- ggplot(data = species_abundance_long, aes(x = X, y = Abundance, fill = Viral_Family)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_y_continuous(limits = c(0, 101), expand = c(0, 0)) + 
  scale_fill_manual(values = mycol2) + 
  theme_minimal() +
  labs(x = "Sample",
       y = "Relative abundance (%)",
       fill = "Viral Family") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "plain", colour = "black"),  
    axis.title.x = element_text(size = 12, face = "plain"), 
    axis.text.y = element_text(size = 10, face = "plain", hjust = 0.5, colour = "black"),
    axis.title.y = element_text(size = 12, face = "plain"), 
    legend.text = element_text(size = 10, face = "plain"),
    legend.title = element_blank(),
    plot.title = element_text(vjust = 0, hjust = 0.5, size = 14, lineheight = 1.4),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid = element_blank(),  
    axis.ticks.y = element_line(color="black", size=0.8, lineend = 0.8), 
    axis.ticks.x = element_line(color="black", size=0.8, lineend = 0.8)) +
  theme(legend.position = "none" ) +
  xlab(NULL) +
  coord_cartesian(ylim = c(-0.2, 100.2)) +  
  scale_y_continuous(expand = c(0, 0))   

