## Script for Figure 1 and Supplementary figure 2, SIC and ICR heatmaps & stacked bar plots
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required packages
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(viridis)

#### Load data ####
load("../analysis_files/ICR_and_SIC_clusters.RData")

# Assign colors for the heatmap annotations and the bar plots
anno_colors <- list(Diagnosis = c("DDLPS" = viridis(9)[1], 
                                   "MFS" = "darkorchid",
                                   "STLMS" = viridis(9)[7],
                                   "ULMS" = viridis(9)[9],
                                   "UPS" = "darkgoldenrod3",
                                   "MPNST" = viridis(9)[3],
                                   "SS" = viridis(9)[5]),
                     Dataset = c("LCCO" = "black",
                                 "TCGA" = "grey"),
                     ICR = c("High" = "#E41A1C", "Medium" = "#4DAF4A", "Low" = "#377EB8"),
                     grade = c("high" = "#134980", "low" = "#FFFF99"),
                     depth = c("deep" = "#80cdc1", "superficial" = "#a6611a", "unknown" = "white"),
                     SIC = c("A" = "#2b7ab3", "B" = "#acd7e6", "C" = "#359e46", "D" = "#f8ae61", "E" = "#d61f26"))

#### Supplementary figure 2A ####
# Heatmap of the SIC classiciation of all samples
# Create top annotation for the heatmap
ha1 <- HeatmapAnnotation(df = design_all[, c("Diagnosis", "Dataset", "depth", "grade", "SIC")], 
                         col = anno_colors)

# Plot the heatmap with the SIC classification
Heatmap(Density.MCP,
        name = "Expression (Z-Score)",
        column_title = "Samples",
        row_title = "Cell populations",
        show_column_names = FALSE,
        top_annotation = ha1,
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

#### Supplementary figure 2B ####
# Stacked bar plots of the SIC classifications
ggplot(data = per_subtype_SIC)+
  geom_bar(aes(x = Diagnosis_factor, y = perc, fill = SIC), stat = "identity")+
  scale_fill_manual(name = "SIC", 
                    breaks = c("E", "D", "C", "B", "A"), 
                    values = c("#d61f26", "#f8ae61", "#359e46", "#acd7e6", "#2b7ab3"))+
  facet_grid(~Dataset, scales = "free", space = "free")+
  geom_text(aes(x = Diagnosis_factor, y = perc, label = n), position = position_stack(vjust = 0.5))+
  xlab("Diagnosis")+
  ylab("Percentage (%)")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.text = element_text(color = "black"))

#### Figure 1A ####
# Heatmap of the SIC classiciation of all samples
# Create top annotation for the heatmap
ha2 <- HeatmapAnnotation(df = design_all[, c("Diagnosis", "Dataset", "depth", "grade", "SIC", "ICR")], 
                         col = anno_colors)

# Plot the heatmap with the ICR clusters
Heatmap(Density.ICR,
        name = "Expression (Z-Score)",
        column_title = "Samples",
        row_title = "ICR genes",
        show_column_names = FALSE,
        top_annotation = ha2,
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

#### Figure 1B ####
# Stacked bar plots of the ICR clusters
ggplot(data = per_subtype_ICR)+
  geom_bar(aes(x = Diagnosis_factor, y = perc, fill = ICR), stat = "identity")+
  scale_fill_manual(name = "ICR cluster", 
                    breaks = c("High", "Medium", "Low"), 
                    values = c("#E41A1C", "#4DAF4A", "#377EB8"))+
  facet_grid(~Dataset, scales = "free", space = "free")+
  geom_text(aes(x = Diagnosis_factor, y = perc, label = n), position = position_stack(vjust = 0.5))+
  xlab("Diagnosis")+
  ylab("Percentage (%)")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45),
        axis.text = element_text(color = "black"))
