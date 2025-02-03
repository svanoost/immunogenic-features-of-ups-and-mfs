## Script for Supplementary figure 1, associations with pathologic response
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required packages
library(ggplot2)
library(ggpubr)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
clinical_data <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/STS_clinical_data.tsv", header = TRUE)

# Assign colors for the figures
anno_colors <- list(diagnosis = c("USTS" = "darkgoldenrod3", "MFS" = "darkorchid"),
                     response = c("Good" = "#1A9641", "Poor" = "#D7191C"))

#### Supplementary figure 1A ####
# Association between fibrosis after radiotherapy, tumor depth and subtype
g1 <- ggplot(data = clinical_data, 
       aes(x = response, y = Fibrosis, fill = response))+
  geom_boxplot()+
  scale_fill_manual(name = "Pathologic response", 
                    breaks = c("Good", "Poor"), 
                    values = anno_colors$response)+
  facet_grid(~diagnosis+depth, scales = "free", space = "free")+
  theme_bw()+
  ylab("Fibrosis (%)")+
  xlab(NULL)+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")

# Association between necrosis after radiotherapy, tumor depth and subtype
g2 <- ggplot(data = clinical_data, 
       aes(x = response, y = Necrosis, fill = response))+
  geom_boxplot()+
  scale_fill_manual(name = "Pathologic response", 
                    breaks = c("Good", "Poor"), 
                    values = anno_colors$response)+
  facet_grid(~diagnosis+depth, scales = "free", space = "free")+
  theme_bw()+
  ylab("Necrosis (%)")+
  xlab(NULL)+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom")

ggarrange(g1, g2, nrow = 2)

#### Supplementary figure 1B ####
# Association between max tumor size and tumor depth in myxofibrosarcoma specifically
ggplot(data = clinical_data[clinical_data$diagnosis == "MFS",], 
       aes(x = response, y = size, fill = response))+
  geom_boxplot()+
  scale_fill_manual(name = "Pathologic response", 
                    breaks = c("Good", "Poor"), 
                    values = anno_colors$response)+
  facet_grid(~depth, scales = "free", space = "free")+
  theme_bw()+
  ylab("Max tumor size (cm)")+
  xlab(NULL)+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom")
