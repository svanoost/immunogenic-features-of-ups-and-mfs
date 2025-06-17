### Script for comparing pre- vs post-treatment cell densities from the imaging mass cytometry data
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(reshape2)
library(dplyr)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
# data frame with the cell counts for all samples
cell_counts <- read.delim("../input_files/IMC_phenotype_counts.tsv", 
                          header = TRUE)

#  design file with the sample annotations of pre vs post 
design_file <- read.delim("../input_files/IMC_pre_post_sample_annotation.tsv", 
                          header = TRUE)

# data frame with the identified phenotypes and cell types
phenotypes <- read.delim("../input_files/phenotype_annotations.tsv", 
                         header = TRUE)

#### Merge data frames and prepare for visualization ####
# reshape the cell count data frame and make values numeric
df <- melt(cell_counts, id.vars = "phenotype", variable.name = "L_ID", value.name = "count")
df$count <- as.numeric(df$count)

# Join the sample annotation with the cell counts
df <- left_join(design_file, df, "L_ID")

# Join the samples with the phenotype annotatio
df <- left_join(df, phenotypes, "phenotype")

# Make a factor of the treatment status for visualization
df$status <- factor(df$status, levels = c("untreated", "treated"))

# Specify the myeloid cells and T cells separately
myeloid <- phenotypes[phenotypes$cell_type == "Macrophages" | phenotypes$cell_type == "Monocytes", "phenotype"]
Tcells <- phenotypes[phenotypes$cell_type == "Tcells", "phenotype"]

#### Perform statistical tests ####
# Perform paired t-test with a Benjamini-Hochberg correction for myxofibrosarcoma ####
# Create data frame for testing
test <- data.frame(diagnosis = NA, phenotype = NA, mean.untreated = NA, mean.treated = NA, t.test = NA, p.val = NA)

# Create a vector of the phenotypes
phenotype_list <- unique(df$phenotype)

# Loop through the different phenotypes
for(i in 1:length(phenotype_list)){
  the_test <- t.test(df[df$Diagnosis == "MFS" & df$phenotype == phenotype_list[i] &
                          df$status == "untreated", "count"], 
                     df[df$Diagnosis == "MFS" & df$phenotype == phenotype_list[i] &
                          df$status == "treated", "count"], paired = TRUE)
  test[i, "diagnosis"] = "MFS"
  test[i, "phenotype"] = phenotype_list[i]
  test[i, "mean.untreated"] = mean(df[df$Diagnosis == "MFS" & df$phenotype == phenotype_list[i] &
                                        df$status == "untreated", "count"])
  test[i, "mean.treated"] = mean(df[df$Diagnosis == "MFS" & df$phenotype == phenotype_list[i] &
                                      df$status == "treated", "count"])
  test[i, "t.test"] = the_test$statistic
  test[i, "p.val"] = the_test$p.value
}

# Perform correction for multiple comparisons
test$FDR <- p.adjust(test$p.val, method = "fdr", n = nrow(test))

# Create separate vector with the significant phenotypes in myxofibrosarcoma
signif <- test[test$FDR < 0.05, "phenotype"]
signif                                         # no significant changes were identified in myxofibrosarcoma

# Save the statistics of the pre- vs post-treatment cell densities
write.table(test, "../output_files/MFS_Statistics_Pre_vs_Post_Phenotypes.txt",
            sep = "\t", row.names = FALSE)

# Perform paired t-test with a Benjamini-Hochberg correction for UPS ####
# Create data frame for testing
test <- data.frame(diagnosis = NA, phenotype = NA, mean.untreated = NA, mean.treated = NA, t.test = NA, p.val = NA)

# Loop through the different phenotypes
for(i in 1:length(phenotype_list)){
  the_test <- t.test(df[df$Diagnosis == "UPS" & df$phenotype == phenotype_list[i] &
                                 df$status == "untreated", "count"], 
                     df[df$Diagnosis == "UPS" & df$phenotype == phenotype_list[i] &
                                 df$status == "treated", "count"], paired = TRUE)
  test[i, "diagnosis"] = "UPS"
  test[i, "phenotype"] = phenotype_list[i]
  test[i, "mean.untreated"] = mean(df[df$Diagnosis == "UPS" & df$phenotype == phenotype_list[i] &
                                df$status == "untreated", "count"])
  test[i, "mean.treated"] = mean(df[df$Diagnosis == "UPS" & df$phenotype == phenotype_list[i] &
                                df$status == "treated", "count"])
  test[i, "t.test"] = the_test$statistic
  test[i, "p.val"] = the_test$p.value
}

# Perform correction for multiple comparisons
test$FDR <- p.adjust(test$p.val, method = "fdr", n = nrow(test))

# Create separate vector with the significant phenotypes in UPS
signif <- test[test$FDR < 0.05, "phenotype"]
signif

# Save the statistics of the pre- vs post-treatment cell densities
write.table(test, "../output_files/UPS_Statistics_Pre_vs_Post_Phenotypes.txt",
            sep = "\t", row.names = FALSE)

#### Save data ####
# Save the data for paired boxplots 
save(df, signif, myeloid, Tcells, 
     file = "../analysis_files/Pre_vs_Post_boxplots.RData")
