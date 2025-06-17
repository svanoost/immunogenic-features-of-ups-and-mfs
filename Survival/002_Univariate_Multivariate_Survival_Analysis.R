### Script to perform univariate and multivariate survival analyses
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(survminer)
library(survival)

# Set working directory
master.location <- setwd(master.location)

# Load clinical data of the whole cohort
clinical_data <- read.delim("../input_files/STS_clinical_final.txt", header = TRUE)

#### Univariate analysis ####
# No clinical characteristics were found significant for UPS
# Analysis was performed for: 
# Tumor size (numeric), 
# Gender at birth (F/M), 
# Age at diagnosis (numeric), 
# Tumor depth (Deep/Superficial),
# Pathologic response to radiotherapy (Good/Poor), 
# Resection margin (complete (R0)/ Incomplete (R1 & R2)), 
# Anatomical location of tumor (Upper extremities/Lower extremities/Other (Head & Neck & Trunk))

# Pathologic response to radiotherapy, tumor size and tumor depth were found significant for myxofibrosarcoma
res.cox <- coxph(Surv(MFS, status_MFS) ~ response, data = clinical_data[clinical_data$subtype == "MFS",])
summary(res.cox)

## Log rank P value = 0.004
## Poor response Hazard Ratio = 0.17, P = 0.011*

res.cox <- coxph(Surv(MFS, status_MFS) ~ size, data = clinical_data[clinical_data$subtype == "MFS",])
summary(res.cox)

## Log rank P value = 0.004
## Tumor size Hazard Ratio = 1.15, P = 0.0057**

res.cox <- coxph(Surv(MFS, status_MFS) ~ depth, data = clinical_data[clinical_data$subtype == "MFS",])
summary(res.cox)

## Log rank P value = 0.03
## Superficial depth Hazard Ratio = 0.14, P = 0.059

#### Multivariate analysis ####
# No clinical characteristic was found to be independent & significant for myxofibrosarcoma
res.cox <- coxph(Surv(MFS, status_MFS) ~ response+size+depth, data = clinical_data[clinical_data$subtype == "MFS",])
summary(res.cox)

## Log rank P value = 0.005
## Poor response Hazard Ratio = 0.46, P = 0.36
## Tumor size Hazard Ratio = 1.1, P = 0.15
## Superficial depth Hazard Ratio = 0.22, P = 0.17
