## Script for processing of RNA-seq data including the TCGA samples 
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(biomaRt)
library(dplyr)
library(IOBR)

# Set working directory
master.location <- setwd(master.location)

#### Load LCCO data ####
# Read raw gene expression count data (LCCO)
gene_expr <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/raw_counts_STS.tsv", 
                        sep = "\t", header = TRUE)

# Add a column with Ensembl gene IDs as row names
gene_expr$ensembl_gene_id <- row.names(gene_expr)

# Load the design file containing sample metadata for LCCO
design_STS <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/design_file_STS.tsv", 
                         sep = "\t", header = TRUE)

# Verify that sample IDs in the count table match those in the design file
all(colnames(gene_expr)[1:23] == design_STS$L_ID)  # Returns TRUE if IDs match

#### Annotate genes with BioMart ####
# Initialize a connection to Ensembl BioMart for human genes
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get available filters and attributes (optional; used for reference)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

# Retrieve Ensembl gene IDs, HGNC symbols, and Entrez gene IDs from BioMart
gene_names_biomart <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), 
                            filters = "ensembl_gene_id", 
                            values = gene_expr$ensembl_gene_id, 
                            mart = ensembl)

# Merge count table with annotated gene symbols and Entrez IDs
gene_expr <- left_join(gene_expr, 
                       gene_names_biomart[, c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol")])

#### Load and prepare TCGA data ####
# Load preprocessed TCGA SARC gene expression data
load("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/TCGA_SARC_Processed.RData")

# Convert TCGA gene metadata to a data frame and standardize column names
Des <- as.data.frame(Des)
colnames(Des) <- c("GeneSymbol", "entrezgene_id")
class(Des$entrezgene_id) <- "integer"  # Ensure Entrez gene IDs are integers

# Convert TCGA expression data to a data frame and add Entrez gene IDs
Data <- as.data.frame(Data)
Data$entrezgene_id <- Des$entrezgene_id

#### Filtering and merging the data sets ####
# Remove duplicate genes in the LCCO data (using Ensembl and Entrez IDs)
gene_expr$dups <- duplicated(gene_expr$ensembl_gene_id)
gene_expr <- gene_expr %>% filter(dups == FALSE)
gene_expr$dups <- duplicated(gene_expr$entrezgene_id)
gene_expr_filtered <- gene_expr %>% filter(dups == FALSE)  # Results in 26,660 genes

# Filter LCCO data to retain genes present in the TCGA dataset
gene_expr_filtered <- gene_expr_filtered %>% 
  filter(entrezgene_id %in% Des$entrezgene_id)  # Results in 19,486 genes

# Merge TCGA data with annotated gene symbols from LCCO
Data_filtered <- left_join(Data, gene_expr_filtered[, c("entrezgene_id", "hgnc_symbol")])

#### TPM conversion and additional filtering ####
# Extract LCCO expression data for TPM conversion using sample IDs
gene_expr_convert <- gene_expr_filtered[, design_STS$L_ID]
row.names(gene_expr_convert) <- gene_expr_filtered$entrezgene_id  # Set row names to Entrez IDs

# Convert raw counts to TPM using the IOBR package
gene_expr_TPM <- IOBR::count2tpm(gene_expr_convert, idType = "Entrez")  # Results in 19,223 genes
gene_expr_TPM$hgnc_symbol <- row.names(gene_expr_TPM)  # Add HGNC symbols as a column

# Filter TCGA data to retain genes present in LCCO data after TPM conversion
Data_filtered <- Data_filtered %>% filter(hgnc_symbol %in% row.names(gene_expr_TPM))

# Filter LCCO TPM data to match TCGA data
gene_expr_TPM_filtered <- gene_expr_TPM %>% 
  filter(hgnc_symbol %in% Data_filtered$hgnc_symbol)  # Results in 17,707 genes

#### Merge filtered data sets ####
# Merge LCCO and TCGA datasets by HGNC symbols
RNA_comb <- left_join(gene_expr_TPM_filtered, Data_filtered, by = "hgnc_symbol")

# Remove non-primary tumor samples from TCGA
RNA_comb <- RNA_comb[, !grepl(colnames(RNA_comb), pattern = "\\.1")]

# Set HGNC symbols as row names and remove gene ID columns
row.names(RNA_comb) <- RNA_comb$hgnc_symbol
RNA_comb <- RNA_comb[, !colnames(RNA_comb) == "hgnc_symbol" & 
                       !colnames(RNA_comb) == "entrezgene_id"]

# Standardize TCGA sample names to the first 15 characters
tcga_ids <- colnames(RNA_comb)[24:288]
tcga_ids <- substr(tcga_ids, 1, 15)
colnames(RNA_comb)[24:288] <- tcga_ids

# Filter RNA_comb to retain 206 revised TCGA samples
revised_samples <- read.delim("I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/input_files/revised_TCGA_samplelist.tsv", 
                              header = TRUE)

all_sample_ids <- c(design_STS$L_ID, revised_samples$L_ID)
RNA_comb <- RNA_comb[, colnames(RNA_comb) %in% all_sample_ids]

#### Save combined data ####
# Save the merged and filtered count tables
save(RNA_comb, file = "I:/BEEN/Siddh/LCCO/Git_repositories/STS_final/analysis_files/Combined_Filtered_Counttables_revisedTCGA_LCCO_TPM.RData")
