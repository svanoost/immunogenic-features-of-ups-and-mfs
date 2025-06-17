### Script for visualization of summed immune cell densities from imaging mass cytometry data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggpubr)

# Set working directory
master.location <- setwd(master.location)

##### Load data ####
load("../analysis_files/IMC_summed_cell_populations.RData")

#### Supplementary figure 3A ####
# Boxplots of summed cell populations in low- and high-grade myxofibrosarcoma and UPS
ggboxplot(summed_counts, x = "grade", y = "sum_counts", fill = "Diagnosis")+ 
  facet_wrap(~cell_type, scales = "free_y", nrow = 2)+
  scale_fill_manual(values = c("darkorchid", "darkgoldenrod3"), 
                    breaks = c("MFS", "UPS"))+
  scale_x_discrete(labels = c("Low grade MFS", "High grade MFS", "UPS"))+
  xlab(NULL)+
  ylab("Mean cell density / mm2")+
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", label = "p.signif")+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
