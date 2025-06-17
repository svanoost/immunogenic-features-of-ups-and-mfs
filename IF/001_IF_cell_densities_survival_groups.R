### Script to create groups based on the T cell (IF) and macrophage (IHC) stainings
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(reshape2)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
# Load T cell and macrophage counts from the IF and double IHC stainings
cell_counts <- read.delim("../input_files/IF_phenotype_counts.tsv", header = TRUE)

# Load clinical data of the whole cohort
clinical_data <- read.delim("../input_files/STS_clinical_data.tsv", header = TRUE)

# Reshape the cell count data from wide to long format for analysis
df <- melt(cell_counts, id.vars = "phenotype", variable.name = "L_ID", value.name = "count")

# Define relevant immune cell subtypes to retain for analysis
keep <- c("CD68+", "CD68+CD163+", "total_Tcells", "perc_PD1_CD4", "perc_PD1_CD8")

# Merge with clinical data to associate samples with subtype
df <- left_join(df, clinical_data[, c("L_ID", "diagnosis")])
df <- df[df$phenotype %in% keep,]

# Convert phenotype variable to a factor with specified order
df$pheno_factor <- factor(df$phenotype, levels = c("total_Tcells", "perc_PD1_CD4", 
                                                                 "perc_PD1_CD8", "CD68+", "CD68+CD163+"))
# Perform unpaired t-tests ####
# Create data frame for testing myxofibrosarcoma vs UPS
test <- data.frame(phenotype = NA, mean.MFS = NA, mean.UPS = NA, t.test = NA, p.val = NA)

# Loop through the different phenotypes
for(i in 1:length(keep)){
  the_test <- t.test(df[df$diagnosis == "MFS" & df$phenotype == keep[i], "count"], 
                     df[df$diagnosis == "UPS" & df$phenotype == keep[i], "count"])
  test[i, "phenotype"] = keep[i]
  test[i, "mean.MFS"] = mean(df[df$diagnosis == "MFS" & df$phenotype == keep[i], "count"], na.rm = TRUE)
  test[i, "mean.UPS"] = mean(df[df$diagnosis == "UPS" & df$phenotype == keep[i], "count"], na.rm = TRUE)
  test[i, "t.test"] = the_test$statistic
  test[i, "p.val"] = the_test$p.value
}

# Save the statistical test results
write.table(test, "../output_files/Statistics_IF_phenotypes.txt",
            sep = "\t", row.names = FALSE)

# Perform unpaired t-test to compare PD-1 positivity between CD8- and CD8+ T cells
the_test <- t.test(df[df$phenotype == "perc_PD1_CD4", "count"], 
                   df[df$phenotype == "perc_PD1_CD8", "count"])
test_PD1 <- data.frame()
test_PD1[1, "phenotype"] = "% of PD-1+CD3+CD8- vs PD-1+CD3+CD8+"
test_PD1[1, "mean.%.CD8-"] = mean(df[df$phenotype == "perc_PD1_CD4", "count"], na.rm = TRUE)
test_PD1[1, "mean.%.CD8+"] = mean(df[df$phenotype == "perc_PD1_CD8", "count"], na.rm = TRUE)
test_PD1[1, "t.test"] = the_test$statistic
test_PD1[1, "p.val"] = the_test$p.value

# Save the statistical test results
write.table(test_PD1, "../output_files/Statistics_IF_PD1.txt",
            sep = "\t", row.names = FALSE)

#### Define groups based on median densities for survival analysis ####
# Transpose and clean data for survival grouping
groups <- as.data.frame(t(cell_counts[, -1]))
colnames(groups) <- cell_counts$phenotype
groups$L_ID <- row.names(groups)

# Merge with clinical data to associate samples wiht subtype
groups <- left_join(groups[, c("L_ID", "CD68+", "CD68+CD163+", "total_Tcells")], 
                       clinical_data[, c("L_ID", "diagnosis")])

# Compute median T cell density for both subtypes
UPS_median <- median(groups[groups$diagnosis == "UPS", "total_Tcells"], na.rm = TRUE)
MFS_median <- median(groups[groups$diagnosis == "MFS", "total_Tcells"], na.rm = TRUE)

# Set all samples at "low"
groups$Tcell_groups <- "Low"

# Assign "high" to samples above the median
groups[groups$diagnosis == "UPS" & groups$total_Tcells > UPS_median, "Tcell_groups"] <- "High"
groups[groups$diagnosis == "MFS" & groups$total_Tcells > MFS_median, "Tcell_groups"] <- "High"

# Merge groups back to clinical data
# Do this here since 1 sample has no macrophage information
clinical_data <- left_join(clinical_data, groups[, c("L_ID", "Tcell_groups")])

# Remove sample without macrophage data
groups <- groups[!is.na(groups$`CD68+`),]

# Compure median single positive macrophage counts for each subtype
UPS_median_single <- median(groups[groups$diagnosis == "UPS", "CD68+"], na.rm = TRUE)
MFS_median_single <- median(groups[groups$diagnosis == "MFS", "CD68+"], na.rm = TRUE)

# Set all samples at "low"
groups$CD68_single_groups <- "Low"

# Assign "high" to samples above the median
groups[groups$diagnosis == "UPS" & groups$`CD68+` > UPS_median_single, 
       "CD68_single_groups"] <- "High"
groups[groups$diagnosis == "MFS" & groups$`CD68+` > MFS_median_single, 
       "CD68_single_groups"] <- "High"

# Compute median double positive macrophage counts for each subtype
UPS_median_double <- median(groups[groups$diagnosis == "UPS", "CD68+CD163+"], na.rm = TRUE)
MFS_median_double <- median(groups[groups$diagnosis == "MFS", "CD68+CD163+"], na.rm = TRUE)

# Set all samples at "low"
groups$CD68_double_groups <- "Low"

# Assign "high" to samples above the median
groups[groups$diagnosis == "UPS" & groups$`CD68+CD163+` > UPS_median_double, 
       "CD68_double_groups"] <- "High"
groups[groups$diagnosis == "MFS" & groups$`CD68+CD163+` > MFS_median_double, 
       "CD68_double_groups"] <- "High"

# Merge macrophage groups with clinical data
clinical_data <- left_join(clinical_data, 
                           groups[, c("L_ID", "CD68_single_groups", "CD68_double_groups")])

# Save the assigned survival groups separately
write.table(clinical_data[, c("L_ID", "diagnosis", "Tcell_groups", "CD68_single_groups", "CD68_double_groups")], 
            "../output_files/IF_survival_groups.txt",
            sep = "\t", row.names = FALSE)

# Separate the clinical data for further analysis
surv_Tcells <- clinical_data[!is.na(clinical_data$Tcell_groups),]
surv_macrophages <- clinical_data[!is.na(clinical_data$CD68_double_groups),]

# Create groups for combined T cell and macrophage groups
# Set all samples at "High and low"
# which are samples with either group high and the other low
surv_macrophages$combined_groups <- "High and low"

# Assign "both high"
surv_macrophages[surv_macrophages$Tcell_groups == "High" & 
                   surv_macrophages$CD68_double_groups == "High", "combined_groups"] <- "Both high"

# Assign "both low"
surv_macrophages[surv_macrophages$Tcell_groups == "Low" & 
                   surv_macrophages$CD68_double_groups == "Low", "combined_groups"] <- "Both low"

#### Create data frame an matrix for heatmap visualization ####
# Create the row annotation for the heatmap
ann_row <- as.data.frame(groups[!is.na(groups$`CD68+`), "diagnosis"])
colnames(ann_row) <- "Diagnosis"
row.names(ann_row) <- groups[!is.na(groups$`CD68+`), "L_ID"]

# Create the matrix for the heatmap visualization
mat <- as.matrix(groups[!is.na(groups$`CD68+`), c("total_Tcells", "CD68+CD163+")])
row.names(mat) <- groups[!is.na(groups$`CD68+`), "L_ID"]

# Transpose the matrix for the calculation of the Z-score
mat_t <- t(mat)

# Calculate the Z-score for the phenotypes
Density.Z = mat_t
for(j in 1: nrow(Density.Z))  {
  Density.Z[j,] = (mat_t[j,]-mean(mat_t[j,]))/sd(mat_t[j,])
}

#### Save data ####
# Save the data for the cell count boxplots
save(df, 
     file = "../analysis_files/IF_cell_counts_boxplots.RData")

# Save the data for the cell count heatmap
save(ann_row, Density.Z, 
     file = "../analysis_files/IF_cell_counts_heatmap.RData")

# Save data for the survival analysis with the assigned groups
save(surv_Tcells, surv_macrophages, 
     file = "../analysis_files/IF_survival_groups.RData")
