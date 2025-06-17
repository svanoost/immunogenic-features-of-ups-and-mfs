## Script for normalizing the RNA-seq data with EDASeq, including the TCGA samples 
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required packages
library(dplyr)
library(EDASeq)
library(base64enc)
library(preprocessCore)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
# Load TPM count table
load("../analysis_files/Combined_Filtered_Counttables_revisedTCGA_LCCO_TPM.RData")

# Load gene information (metadata)
load("../input_files/GeneInfo.Rdata")
geneInfo <- as.data.frame(geneInfo)
geneInfo <- geneInfo[!is.na(geneInfo[, 1]), ]                        # Remove genes without associated information

#### Prepare data set ####
# Rename the loaded RNA expression dataset
gene_expr <- RNA_comb

# Extract gene names (rows of the gene expression matrix)
genes <- as.data.frame(row.names(gene_expr))
colnames(genes) <- "hgnc_symbol"

# Add GC content and filter out duplicates and genes without GC content
genes <- left_join(genes, geneInfo[, c("gc", "hgnc_symbol")])
genes$dups <- duplicated(genes$hgnc_symbol)                          # Identify duplicate gene symbols
genes <- genes[genes$dups != "TRUE", ]                               # Remove duplicate genes
genes_filtered <- genes[!is.na(genes$gc), ]                          # Retain only genes with GC content

# Subset the expression data to include only the filtered genes
expression_filtered <- gene_expr[genes_filtered$hgnc_symbol, ]

#### Filter genes with low expression ####
# Filter genes expressed in at least one LCCO sample (columns 1 to 23 of the dataset)
length(which(apply(expression_filtered[, 1:23], 1, function(x) { any(x > 3) })))  # Count genes with expression > 3
idx <- which(apply(expression_filtered[, 1:23], 1, function(x) { any(x > 3) }))   # Get indices of such genes
gene_expr_filtered <- expression_filtered[idx, ]                                  # Filter the expression data for these genes

#### Match filtered genes with metadata ####
# Retain only genes present in both the expression data and gene information
available.genes <- unique(row.names(gene_expr_filtered))               # Get unique filtered genes
geneInfo <- geneInfo[geneInfo$hgnc_symbol %in% available.genes, ]      # Subset gene info for available genes
geneInfo$dups <- duplicated(geneInfo$hgnc_symbol)                      # Identify duplicates in gene info
geneInfo <- geneInfo %>% filter(dups == "FALSE")                       # Remove duplicates
row.names(geneInfo) <- geneInfo$hgnc_symbol                            # Set row names to gene symbols
geneInfo <- geneInfo[, !colnames(geneInfo) == "dups"]                  # Remove the 'dups' column

# Convert the filtered expression data to a numeric matrix
expression.filtered <- as.matrix(gene_expr_filtered[row.names(geneInfo), ])
mode(expression.filtered) <- "numeric"                                 # Ensure numeric mode for calculations
dim(expression.filtered)                                               # Check dimensions (14369 genes x 229 samples)

#### Save filtered data ####
# Save the filtered expression matrix for future use
save(expression_filtered, file = "../analysis_files/Filtered_Gene_Expression_Matrix.RData")

#### Normalization with EDAseq ####
# Create a SeqExpressionSet object for normalization
RNASeq.expr.set <- newSeqExpressionSet(expression.filtered, featureData = geneInfo)

# Ensure the GC content is numeric
fData(RNASeq.expr.set)[, "gcContent"] <- as.numeric(geneInfo[, "gc"])

# Perform within-lane normalization to adjust for gene-specific effects (e.g., GC content)
RNASeq.expr.set <- withinLaneNormalization(RNASeq.expr.set, "gcContent", which = "upper", offset = TRUE)

# Perform between-lane normalization to adjust for differences in sequencing depth
RNASeq.expr.set <- betweenLaneNormalization(RNASeq.expr.set, which = "upper", offset = TRUE)

# Apply the offset from EDAseq to the raw expression data and log-transform
RNASeq.NORM <- log(expression.filtered + .1) + offst(RNASeq.expr.set)
RNASeq.NORM <- floor(exp(RNASeq.NORM) - .1)                            # Convert back to non-decimal integer values

#### Quantile Normalization ####
# Apply quantile normalization to ensure uniform distribution across samples
RNASeq.NORM.quantiles <- normalize.quantiles(RNASeq.NORM)
RNASeq.NORM.quantiles <- floor(RNASeq.NORM.quantiles)                 # Convert to integer values

# Preserve row and column names after normalization
row.names(RNASeq.NORM.quantiles) <- row.names(RNASeq.NORM)
colnames(RNASeq.NORM.quantiles) <- colnames(RNASeq.NORM)

#### Log transformation ####
# Perform log2 transformation to stabilize variance
RNASeq.QN.LOG2 <- log(RNASeq.NORM.quantiles + 1, 2)

#### Save normalized data ####
# Save the normalized expression matrix for future use
save(RNASeq.QN.LOG2, file = "../analysis_files/Log2_Quantile_Normalized_Expression.RData")
