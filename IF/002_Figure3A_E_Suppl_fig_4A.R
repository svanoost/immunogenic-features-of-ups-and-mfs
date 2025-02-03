### Script for visualizing the T cell (IF) and macrophage (IHC) data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
# Load cell counts from the IF and IHC stainings
load("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/IF_cell_counts_boxplots.RData")
load("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/IF_cell_counts_heatmap.RData")

# Load T cell counts from IMC vs IF for correlation
IMC_IF <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/IMC_vs_IF_Tcell_counts.tsv", header = TRUE)

#### Figure 3A ####
# Plot the boxplots with the T cells counts, % of PD-1 positivity and the macrophage counts
# statistics were added manually
ggplot(data = df, aes(x = phenotype, y = count, fill = diagnosis))+
  geom_boxplot()+
  xlab("Phenotype")+
  ylab(NULL)+
  facet_wrap(~pheno_factor, scales = "free",
             strip.position = "left",
             labeller = as_labeller(c(total_Tcells = "Mean cell density / mm2", 
                                      perc_PD1_CD4 = "PD-1 positivity (%)",
                                      perc_PD1_CD8 = "PD-1 positivity (%)",
                                      `CD68+` = "Mean cell density / mm2",
                                      `CD68+CD163+` = "Mean cell density / mm2")))+
  scale_fill_manual(name = "Diagnosis", 
                    breaks = c("MFS", "USTS"), 
                    values = c("darkorchid", "darkgoldenrod3"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.placement = "outside")

#### Figure 3E ####
# Plot the heatmap showing the association between macrophage and T cell infiltration
# Create the row annotation for the samples
ha_1 <- rowAnnotation(df = ann_row, 
                      col = list(Diagnosis = c("MFS" = "darkorchid", "USTS" = "darkgoldenrod3")),
                      na_col = "white")


# Create the heatmap for double positive macrophages and total T cells, separated per subtype
Heatmap(t(Density.Z),
        name = "Density (Z-Score)",
        border = FALSE,
        row_split = ann_row$Diagnosis,
        show_column_names = TRUE,
        show_row_names = FALSE,
        left_annotation = ha_1,
        row_gap = unit(2, "mm"),
        row_title = "Phenotypes",
        column_title = NULL,
        clustering_method_rows = "ward.D",
        column_names_rot = -45,
        cluster_columns = FALSE,
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

#### Supplementary figure 4A #### 
# Plot the correlation between IMC and IF T cell counts
# Calculate the Pearson correlation and statistical significance
cor <- cor.test(IMC_IF[, "total_Tcells_IF"], IMC_IF[, "total_Tcells_IMC"])

# Create the correlation and P value annotation
labels <- c(paste0("~italic(R)^2 ==", round(cor$estimate, 2)),
             paste0("~italic(P) ==", round(cor$p.value, 5)))

# Create the correlation plot
ggplot(data = IMC_IF, aes(x = total_Tcells_IF, y = total_Tcells_IMC))+
  geom_smooth(method = "lm", color  ="black", alpha = 0.3)+
  geom_point()+
  xlab("total T cell density (IF)")+
  ylab("total T cell density (IMC)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")+
  annotate(geom = "text", x = 0, y = c(900, 0.95*900), 
           label = labels, parse = TRUE,
           hjust = 0)
