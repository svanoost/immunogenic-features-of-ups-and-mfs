### Script for comparing pre- vs post-treatment cell densities from the imaging mass cytometry data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggpubr)

#### Load data ####
load("../analysis_files/Pre_vs_Post_boxplots.RData")

#### Figure 4 ####
# Create paired boxplots for myeloid cells in UPS
g1 <- ggpaired(df[df$Diagnosis == "UPS" & 
                     df$phenotype %in% myeloid,], 
         x = "status", y = "count", id = "Pat_ID", color = "status",
         line.color = "gray", line.size = 0.4, linetype = "dashed")+ 
  facet_wrap(~phenotype, scales = "free_y", ncol = 3)+
  scale_color_manual(name = "Sample", 
                     values = c("#00BBA7", "#008273"), 
                     breaks = c("untreated", "treated"))+
  xlab("Phenotype")+
  ylab("Mean cell density / mm2")+
  stat_compare_means(paired = TRUE, method = "t.test", label.y = 0)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")

# Create paired boxplots for T cells in UPS
g2 <- ggpaired(df[df$Diagnosis == "UPS" & 
                           df$phenotype %in% Tcells,], 
               x = "status", y = "count", id = "Pat_ID", color = "status",
               line.color = "gray", line.size = 0.4, linetype = "dashed")+ 
  facet_wrap(~phenotype, scales = "free_y", ncol = 2)+
  scale_color_manual(name = "Sample", 
                     values = c("#00BBA7", "#008273"),
                     breaks = c("untreated", "treated"))+
  xlab("Phenotype")+
  ylab("Mean cell density / mm2")+
  stat_compare_means(paired = TRUE, method = "t.test", label.y = 0)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")

# Combine both boxplots and align. Title is only used in the script, not in the final figures
annotate_figure(ggarrange(g1, g2, common.legend = T, widths = c(1.465, 1)), 
                top = text_grob("Myeloid cells & T cells in UPS", size = 14, face = "bold"))

#### Supplementary figure 5A ####
# Filter for FDR significant phenotypes other than myeloid cells and T cells, which were only significant in UPS. 
# Myxofibrosarcomas are taken along as a comparison.
signif <- signif[6:9]

# Create paired boxplots for other significant phenotypes in UPS
g5 <- ggpaired(df[df$Diagnosis == "UPS" & 
                    df$phenotype %in% signif,], 
               x = "status", y = "count", id = "Pat_ID", color = "status",
               line.color = "gray", line.size = 0.4, linetype = "dashed")+ 
  facet_wrap(~phenotype, scales = "free_y", ncol = 2)+
  scale_color_manual(name = "Sample", 
                     values = c("#00BBA7", "#008273"), 
                     breaks = c("untreated", "treated"))+
  xlab("Phenotype")+
  ylab("Mean cell density / mm2")+
  stat_compare_means(paired = TRUE, method = "t.test", label.y = 0)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")

# Create paired boxplots for other significant phenotypes in myxofibrosarcoma
g6 <- ggpaired(df[df$Diagnosis == "MFS" & 
                    df$phenotype %in% signif,], 
               x = "status", y = "count", id = "Pat_ID", color = "status",
               line.color = "gray", line.size = 0.4, linetype = "dashed")+ 
  facet_wrap(~phenotype, scales = "free_y", ncol = 2)+
  scale_color_manual(name = "Sample", 
                     values = c("#00BBA7", "#008273"), 
                     breaks = c("untreated", "treated"))+
  xlab("Phenotype")+
  ylab("Mean cell density / mm2")+
  stat_compare_means(paired = TRUE, method = "t.test", label.y = 0)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")

# Combine both boxplots and align. Title is only used in the script, not in the final figures
annotate_figure(ggarrange(g5, g6, common.legend = T),
                top = text_grob("Other significant phenotypes in UPS and myxofibrosarcoma", size = 14, face = "bold"))

#### Supplementary figure 5B ####
# Create paired boxplots for myeloid cells in myxofibrosarcoma
g3 <- ggpaired(df[df$Diagnosis == "MFS" & 
                           df$phenotype %in% myeloid,], 
               x = "status", y = "count", id = "Pat_ID", color = "status",
               line.color = "gray", line.size = 0.4, linetype = "dashed")+ 
  facet_wrap(~phenotype, scales = "free_y", ncol = 3)+
  scale_color_manual(name = "Sample", 
                     values = c("#00BBA7", "#008273"), 
                     breaks = c("untreated", "treated"))+
  xlab("Phenotype")+
  ylab("Mean cell density / mm2")+
  stat_compare_means(paired = TRUE, method = "t.test", label.y = 0)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")

# Create paired boxplots for T cells in UPS
g4 <- ggpaired(df[df$Diagnosis == "MFS" & 
                           df$phenotype %in% Tcells,], 
               x = "status", y = "count", id = "Pat_ID", color = "status",
               line.color = "gray", line.size = 0.4, linetype = "dashed")+ 
  facet_wrap(~phenotype, scales = "free_y", ncol = 2)+
  scale_color_manual(name = "Sample", 
                     values = c("#00BBA7", "#008273"), 
                     breaks = c("untreated", "treated"))+
  xlab("Phenotype")+
  ylab("Mean cell density / mm2")+
  stat_compare_means(paired = TRUE, method = "t.test", label.y = 0)+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank(),
        legend.position = "bottom")   

# Combine both boxplots and align. Title is only used in the script, not in the final figures
annotate_figure(ggarrange(g3, g4, common.legend = T, widths = c(1.465, 1)), 
                top = text_grob("Myeloid cells & T cells in myxofibrosarcoma", size = 14, face = "bold"))
