## Script for batch correction of the RNA-seq data with limma
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required packages
library(dplyr)
library(MCPcounter)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
# Load preprocessed and log2 quantile-normalized gene expression data.
load("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/Log2_Quantile_Normalized_Expression.RData")

# Load metadata describing samples (e.g., diagnosis, dataset) and set sample IDs as row names.
design_all <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/design_file_LCCO_TCGA.tsv", header = TRUE)
row.names(design_all) <- design_all$L_ID

#### Batch correction ####
# Correct batch effects using the `limma` package to adjust expression values for dataset-specific biase
mat <- limma::removeBatchEffect(RNASeq.QN.LOG2[, design_all$L_ID], design_all$Dataset)
df <- as.data.frame(mat)

#### MCP counter for Sarcoma Immune Classes (SIC) ####
# Perform MCPcounter analysis on sarcoma-specific genes to estimate immune cell abundance.
MCP.df <- as.matrix(MCPcounter.estimate(df[, design_all$L_ID], 
                                        featuresType = "HUGO_symbols", 
                                        probesets = "SARC"))

# Standardize (z-score) each immune cell population's abundance across samples.
Density.MCP = MCP.df[1:9,] # exclude fibroblasts (column 10)
for(j in 1: nrow(Density.MCP))  {
  Density.MCP[j,] = (MCP.df[j,]-mean(MCP.df[j,]))/sd(MCP.df[j,])
}

# Perform hierarchical clustering to divide samples into 5 SICs
dend <- as.dendrogram(hclust(dist(t(Density.MCP)), method = "ward.D"))
dend_order <- dendextend::cutree(dend, k = 5)

# Create a data frame with sample IDs and assigned SIC clusters
anno_col <- data.frame(L_ID = colnames(Density.MCP), SIC = dend_order)
anno_col$SIC <- as.character(anno_col$SIC)

# Merge the SIC annotations with the design file.
design_all <- left_join(design_all, anno_col)

# Assign descriptive names to SIC clusters
design_all[design_all$SIC == "3", "SIC"] <- "E" # Immune-high
design_all[design_all$SIC == "4", "SIC"] <- "D" # Immune-medium high
design_all[design_all$SIC == "1", "SIC"] <- "C" # Highly vascularized
design_all[design_all$SIC == "5", "SIC"] <- "B" # Immune-low
design_all[design_all$SIC == "2", "SIC"] <- "A" # Immune desert

#### Immunologic constant of rejection (ICR) signature ####
# Define a list of genes associated with the ICR signature.
ICR_genes = c("IFNG", "TBX21", "CD8A", "CD8B", 
              "STAT1", "IRF1","CXCL9", "CXCL10",
              "CCL5","GNLY", "PRF1", "GZMA",
              "GZMB", "GZMH","CD274", "CTLA4",
              "FOXP3", "IDO1") #  IL12B and PDCD1 are not in the data set

# Extract the ICR gene expression data
mat2 <- as.matrix(df[ICR_genes, design_all$L_ID])

# Standardize (z-score) ICR gene expression values across samples.
Density.ICR = mat2
for(j in 1: nrow(Density.ICR))  {
  Density.ICR[j,] = (mat2[j,]-mean(mat2[j,]))/sd(mat2[j,])
}

# Perform hierarchical clustering to divide samples into 3 clusters
dend <- as.dendrogram(hclust(dist(t(Density.ICR)), method = "ward.D"))
dend_order <- dendextend::cutree(dend, k = 3)

# Create a data frame with sample IDs and assigned ICR clusters
anno_col <- data.frame(L_ID = colnames(Density.ICR), ICR = dend_order)
anno_col$ICR <- as.character(anno_col$ICR)

# Merge ICR annotations with the design file
design_all <- left_join(design_all, anno_col)

# Assign descriptive names to ICR clusters
design_all[design_all$ICR == "1", "ICR"] <- "High"
design_all[design_all$ICR == "2", "ICR"] <- "Low"
design_all[design_all$ICR == "3", "ICR"] <- "Medium"

#### Calculate the distribution of all SIC clusters across both datasets ####
# Filter the design data to exclude samples with "low" grade
per_subtype_SIC <- design_all[design_all$grade != "low",] %>% 
  group_by(Diagnosis, Dataset) %>%                      # Group by Diagnosis and Dataset
  dplyr::count(SIC) %>%                                 # Count samples per SIC within each group
  group_by(Diagnosis, Dataset) %>%                      # Re-group by Diagnosis and Dataset
  mutate(total = sum(n)) %>%                            # Compute total samples in each group
  group_by(Diagnosis, Dataset, SIC) %>%                 # Group by Diagnosis, Dataset, and SIC
  mutate(perc = round(n / total * 100, 2))              # Calculate percentage for each SIC

# These percentages are used for Supplementary figure 2B
# Assign a specific order to the SIC categories for consistent visualization
per_subtype_SIC$SIC <- factor(per_subtype_SIC$SIC, levels = c("E", "D", "C", "B", "A"))

# Assign a specific order to the Diagnosis categories to match a desired sequence
per_subtype_SIC$Diagnosis_factor <- factor(per_subtype_SIC$Diagnosis, levels = c("USTS", "MFS", 
                                                                                 "DDLPS",
                                                                                 "STLMS", "ULMS",
                                                                                 "MPNST",
                                                                                 "SS"))

# Same as above but aggregates across all datasets for each Diagnosis
# These percentages are used in the text
per_subtype_SIC_total <- design_all[design_all$grade != "low",] %>% 
  group_by(Diagnosis) %>%                               # Group by Diagnosis only
  dplyr::count(SIC) %>%                                 # Count samples per SIC for each Diagnosis
  group_by(Diagnosis) %>%                               # Re-group by Diagnosis
  mutate(total = sum(n)) %>%                            # Compute total samples for each Diagnosis
  group_by(Diagnosis, SIC) %>%                          # Group by Diagnosis and SIC
  mutate(perc = round(n / total * 100, 2))              # Calculate percentage for each SIC

#### Save aggregated percentages of SIC distribution ####
write.table(per_subtype_SIC_total, "I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/output_files/SIC_percentages.txt", 
            row.names = FALSE, sep = "\t")

#### Calculate the distribution of all ICR clusters across both datasets ####
# Filter the design data to exclude samples with "low" grade
per_subtype_ICR <- design_all[design_all$grade != "low",] %>% 
  group_by(Diagnosis, Dataset) %>%                      # Group by Diagnosis and Dataset
  dplyr::count(ICR) %>%                                 # Count samples per ICR cluster within each group
  group_by(Diagnosis, Dataset) %>%                      # Re-group by Diagnosis and Dataset
  mutate(total = sum(n)) %>%                            # Compute total samples in each group
  group_by(Diagnosis, Dataset, ICR) %>%                 # Group by Diagnosis, Dataset, and ICR cluster
  mutate(perc = round(n/total*100, 2))                  # Calculate percentage for each ICR cluster

# These percentages are used for Figure 1B
# Assign a specific order to the ICR clusters for consistent visualization
per_subtype_ICR$ICR <- factor(per_subtype_ICR$ICR, levels = c("High", "Medium", "Low"))

# Assign a specific order to the Diagnosis categories to match a desired sequence
per_subtype_ICR$Diagnosis_factor <- factor(per_subtype_ICR$Diagnosis, levels = c("USTS", "MFS", 
                                                                         "DDLPS",
                                                                         "STLMS", "ULMS",
                                                                         "MPNST",
                                                                         "SS"))

# Reorder the data frame by Diagnosis and then by ICR cluster
per_subtype_ICR <- per_subtype_ICR %>%
  arrange(desc(ICR))                                    # order by ICR (in the specified order)

# Same as above but aggregates across all datasets for each Diagnosis
# These percentages are used in the text
per_subtype_ICR_total <- design_all[design_all$grade != "low",] %>% 
  group_by(Diagnosis) %>%                               # Group by Diagnosis only
  dplyr::count(ICR) %>%                                 # Count samples per ICR cluster for each Diagnosis
  group_by(Diagnosis) %>%                               # Re-group by Diagnosis
  mutate(total = sum(n)) %>%                            # Compute total samples for each Diagnosis
  group_by(Diagnosis, ICR) %>%                          # Group by Diagnosis and ICR cluster
  mutate(perc = round(n/total*100, 2))                  # Calculate percentage for each ICR cluster

#### Save aggregated percentages of ICR clusters ####
write.table(per_subtype_ICR_total, "I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/output_files/ICR_percentages.txt", 
            row.names = FALSE, sep = "\t")

#### Save data for Figure 1 and Supplementary Figure 2
save(design_all, Density.ICR, Density.MCP, per_subtype_SIC, per_subtype_ICR, 
     file = "I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/ICR_and_SIC_clusters.RData")
