### Script for calculating the relative abundance of immune cell populations and prepare heatmap & boxplots
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# load necessary packages
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
# Previously generated data frame with the cell counts
load("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/IMC_cell_counts_with_annotation.RData")
df[df$phenotype == "Vessels", "clusters"] <- "Rest"

#  design file with the sample annotations
design_file <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/IMC_sample_annotation.tsv", 
                          header = TRUE)

# data frame with the identified phenotypes and cell types
phenotypes <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/phenotype_annotations.tsv", 
                         header = TRUE)

#### Create heatmap annotation data frames ####
# create phenotype annotation to separate the heatmaps
ann_row <- data.frame(df[1:32, "clusters",FALSE])
colnames(ann_row) <- "Cell.Type"
row.names(ann_row) <- df[1:32, "phenotype"]
ann_row[ann_row$Cell.Type == "", "Cell.Type"] <- "Tumor"

# Remove low-grade myxofibrosarcomas
ann_col <- design_file[design_file$grade_numeric != "1",]
row.names(ann_col) <- ann_col$L_ID

#### Calculate the relative abundance of all phenotypes ####
# Calculate the relative phenotype abundance for all phenotypes
rel_counts <- df[df$L_ID %in% ann_col$L_ID, c("L_ID", "Diagnosis", "grade_numeric", "count", "phenotype", "clusters")] %>% 
  group_by(L_ID) %>%
  mutate(total_counts = sum(count)) %>%
  group_by(L_ID, phenotype) %>%
  mutate(rel_counts = count/total_counts*100) %>%
  unique()

# Reshape the data frame for visualization
rel_pivot <- pivot_wider(rel_counts[, c("L_ID", "phenotype", "rel_counts")], names_from = L_ID, values_from = rel_counts)

#### Calculate the Z-score of the relative abundance ####
# Create a matrix to calculate the Z-score per phenotype
mat <- as.matrix(rel_pivot[, -1])
row.names(mat) <- rel_pivot$phenotype

# Check whether sample IDs correspond between the phenotypes and annotations
all(colnames(mat) == ann_col$L_ID)

# Calculate the Z-score
Density.Z = mat
for(j in 1: nrow(Density.Z))  {
  Density.Z[j,] = (mat[j,]-mean(mat[j,]))/sd(mat[j,])
}

#### Specify the annotation colors and the row order of the heatmap ####
# Row order from the heatmap after hierarchical clustering. This will be used to align the heatmap with boxplots
pheno_order <- c("KI67_Tcells", "CD8_KI67_Tcells", "Tcells", "CD8_Tcells", "CD57_Tcells",   # T cells
                 "Tregs", "gd_Tcells", "CD57_ILCs", "ILCs", "KI67_ILCs", "CD56_ILCs",       # Innate lymphoid cells
                 "Bcells", "PlasmaBcells","Vessels", "Granulocytes",                        # Other cell types
                 "Macrophages", "CD163_HLADR_Macrophages", "CD163_Macrophages",             # Macrophages
                 "CD163_CD204_Macrophages", "CD163_CD204_HLADR_Macrophages",                # Macrophages
                 "CD163_Monocytes", "Monocytes", "CD163_CD204_Monocytes",                   # Monocytes
                 "CD163_HLADR_Monocytes", "CD163_CD204_HLADR_Monocytes",                    # Monocytes
                 "HLADR_DCs", "CD11c_HLADR_DCs")                                            # Dendritic cells

# Create the annotation list and colors for the heatmap
anno_colors <- list(Diagnosis = c("USTS" = "darkgoldenrod3", "MFS" = "darkorchid"),
                    response = c("Good" = "#1a9641", "Poor" = "#d7191c"),
                    grade_numeric = c("1" = "#ffff99", "2" = "#99ccff", "3" = "#3333ff"),
                    ICR = c("High" = "#E41A1C", "Medium" = "#4DAF4A", "Low" = "#377EB8"))

#### Create the relative abundance heatmap with boxplots ####
# Cell types will be separated for improved visualization due to the large variety in relative abundance
# Make heatmap 1, with the T cells
# Filter the specified cell type from the row annotation data frame
hmap_1 <- ann_row %>% filter(Cell.Type == "Tcells")

# Retrieve those phenotypes from the Z-score matrix
hmap_1 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_1),] 

# Generate a top annotation for the combined heatmap
hb_1 <- HeatmapAnnotation(df = ann_col[, c("Diagnosis", "response", "grade_numeric", "ICR")],
                      col = anno_colors)

# Create the heatmap for the specified cell type
# Columns are clustered but rows are predefined after an initial hierarchical clustering
h1 <- Heatmap(hmap_1[pheno_order[1:7],],
              name = "Z-Score", 
              border = TRUE,
              column_names_gp = gpar(fontsize = 10),
              column_gap = unit(2, "mm"),
              column_split = ann_col$Diagnosis,
              column_title = " ",
              row_names_side = "left",
              top_annotation = hb_1,
              clustering_method_columns = "ward.D",
              column_names_rot = -45,
              show_column_dend = T,
              cluster_rows = FALSE,
              show_row_names = T,
              show_column_names = T,
              show_row_dend = F,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")) 
)

# Make heatmap 2, with the innate lymphoid cells
hmap_2 <- ann_row %>% filter(Cell.Type == "ILCs")
hmap_2 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_2),] 

h2 <- Heatmap(hmap_2[pheno_order[8:11],],
              name = "Z-Score",
              border = TRUE,
              column_names_gp = gpar(fontsize = 10),
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              row_names_side = "left",
              clustering_method_columns = "ward.D",
              column_names_rot = -45,
              show_column_dend = T,
              cluster_rows = FALSE,
              show_row_names = T,
              show_column_names = T,
              show_row_dend = F,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

# Make heatmap 3, with the other cell types
hmap_3 <- ann_row %>% filter(Cell.Type == "Rest")
hmap_3 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_3),] 

h3 <- Heatmap(hmap_3[pheno_order[12:15],],
              name = "Z-Score",
              border = TRUE,
              column_names_gp = gpar(fontsize = 10),
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              row_names_side = "left",
              clustering_method_columns = "ward.D",
              column_names_rot = -45,
              show_column_dend = T,
              cluster_rows = FALSE,
              show_row_names = T,
              show_column_names = T,
              show_row_dend = F,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

# Make heatmap 4, with the macrophages
hmap_4 <- ann_row %>% filter(Cell.Type == "Macrophages")
hmap_4 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_4),] 

h4 <- Heatmap(hmap_4[pheno_order[16:20],],
              name = "Z-Score",
              border = TRUE,
              column_names_gp = gpar(fontsize = 10),
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              row_names_side = "left",
              clustering_method_columns = "ward.D",
              column_names_rot = -45,
              show_column_dend = T,
              cluster_rows = FALSE,
              show_row_names = T,
              show_column_names = T,
              show_row_dend = F,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

# Make heatmap 5, with the monocytes
hmap_5 <- ann_row %>% filter(Cell.Type == "Monocytes")
hmap_5 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_5),] 

h5 <- Heatmap(hmap_5[pheno_order[21:25],],
              name = "Z-Score",
              border = TRUE,
              column_names_gp = gpar(fontsize = 10),
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              row_names_side = "left",
              clustering_method_columns = "ward.D",
              column_names_rot = -45,
              show_column_dend = T,
              cluster_rows = FALSE,
              show_row_names = T,
              show_column_names = T,
              show_row_dend = F,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

# Make heatmap 6, with the dendritic cells
hmap_6 <- ann_row %>% filter(Cell.Type == "DCs")
hmap_6 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_6),]

h6 <- Heatmap(hmap_6[pheno_order[26:27],],
              name = "Z-Score",
              border = TRUE,
              column_names_gp = gpar(fontsize = 10),
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              row_names_side = "left",
              clustering_method_columns = "ward.D",
              column_names_rot = -45,
              show_column_dend = T,
              cluster_rows = FALSE,
              show_row_names = T,
              show_column_names = T,
              show_row_dend = F,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
)

#### Combine the heatmaps ###X
ht_list <- h1 %v% h2 %v% h3 %v% h4 %v% h5 %v% h6

#### Order data frame for boxplot visualization ####
# remove all tumor/stroma cell populations
rel_counts <- rel_counts[!grepl(rel_counts$phenotype, pattern = "Tumor"),]

# Make a factor of all the phenotypes in the same order as the heatmap
rel_counts$pheno_factor <- factor(rel_counts$phenotype, levels = pheno_order)

# Make a factor of the cell type clusters to match the heatmap
rel_counts$cluster_factor <- factor(rel_counts$clusters, levels = c("Tcells", "ILCs", 
                                                                    "Rest", "Macrophages",
                                                                    "Monocytes", "DCs"))

#### Save data ####
# Heatmap & data for boxplots
save(ht_list, rel_counts, file = "I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/IMC_heatmap_boxplots.RData")
