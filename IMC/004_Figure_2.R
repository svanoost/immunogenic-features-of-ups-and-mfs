### Script for visualizing the relative abundance heatmap & boxplots
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# load necessary packages
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

### Load data ####
load("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/IMC_heatmap_boxplots.RData")

#### Figure 2 ####
# Plot the relative abundance heatmap
draw(ht_list, row_title = "Phenotypes", row_title_gp = gpar(font = 1, fontsize = 10), 
     column_title = "Samples", column_title_gp = gpar(font = 1, fontsize = 10),
     heatmap_legend_side = "bottom",
     adjust_annotation_extension = FALSE)

# Plot the accompanying relative abundance boxplots, separated per cell type
ggplot(data = rel_counts, aes(x = pheno_factor, y = rel_counts, fill = Diagnosis))+
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(breaks = c("MFS", "USTS"), values = c("darkorchid", "darkgoldenrod3"))+
  scale_x_discrete(labels = NULL)+
  xlab(NULL)+
  ylab("Relative abundance (%)")+
  facet_grid(~cluster_factor, scales = "free", space = "free")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid = element_blank())
