### Script for calculations to sum immune cell densities from imaging mass cytometry data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(reshape2)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
# data frame with the cell counts for all samples
cell_counts <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/IMC_phenotype_counts.tsv", 
                          header = TRUE)

#  design file with the sample annotations
design_file <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/IMC_sample_annotation.tsv", 
                          header = TRUE)

# data frame with the identified phenotypes and cell types
phenotypes <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/phenotype_annotations.tsv", 
                         header = TRUE)

#### Merge data frames and prepare for visualization ####
# reshape the cell count data frame and make values numeric
df <- melt(cell_counts, id.vars = "phenotype", variable.name = "L_ID", value.name = "count")
df$count <- as.numeric(df$count)

# Join the sample annotation with the cell counts
df <- left_join(design_file, df, "L_ID")

# Join the samples with the phenotype annotation
df <- left_join(df, phenotypes, "phenotype")

# Make a factor of the grading for visualization
df$grade <- factor(df$grade, levels = c("Low_MFS", "High_MFS", "High_USTS"))

# Group all phenotypes from major cell populations and sum the cell counts
summed_counts <- df[, c("L_ID", "Diagnosis", "grade", "count", "cell_type")] %>% 
  group_by(L_ID, cell_type) %>%
  mutate(sum_counts = sum(count)) %>%
  select (-count) %>%
  unique()

# Make a factor of the cell types for visualization
summed_counts$cell_type <- factor(summed_counts$cell_type, levels = c("Bcells", "PlasmaBcells", "Tcells", "ILCs",
                                                                      "Granulocytes", "DCs", "Macrophages", "Monocytes",
                                                                      "Vessels", "Tumor"))

# Specify the comparisons to be made for statistical testing
my_comparisons <- list( c("Low_MFS", "High_MFS"), c("High_MFS", "High_USTS") )

#### Save data ####
# save data for Supplementary figure 3A
save(summed_counts, my_comparisons, 
     file = "I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/IMC_summed_cell_populations.RData")

# save data for figure 2
save(df, file = "I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/IMC_cell_counts_with_annotation.RData")
