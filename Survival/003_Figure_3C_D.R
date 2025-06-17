### Script for visualizing the metastasis-free survival based on T cell and double positive macrophage groups
# R version 4.4.0

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(survival)
library(survminer)

# Set working directory
master.location <- setwd(master.location)

#### Load data ####
load("../analysis_files/IF_survival_groups.RData")

#### Figure 3C ####
# Perform the metastasis-free survival analysis for the T cell groups in myxofibrosarcomas
fit <- survfit(Surv(MFS, status_MFS) ~ Tcell_groups, 
               data = surv_Tcells[surv_Tcells$diagnosis == "MFS",])

# Check the summary of the survival analysis
fit

# Plot the Kaplan-Meier curve
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("T cell high", "T cell low"),
           font.legend = c(12),
           risk.table.font = c(5),
           title = "Metastasis-free survival in MFS",
           legend.title = "",
           censor.size = 6,
           xlim = c(0, 90), # use 90 months as max
           break.x.by = 15, # use 15 months as breaks
           palette = c("#E41A1C", "#377EB8"),
)+
  xlab(label = "Time (months)")

# Perform the metastasis-free survival analysis for the T cell groups in UPS
fit <- survfit(Surv(MFS, status_MFS) ~ Tcell_groups, 
               data = surv_Tcells[surv_Tcells$diagnosis == "UPS",])

# Check the summary of the survival analysis
fit

# Plot the Kaplan-Meier curve
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("T cell high", "T cell low"),
           font.legend = c(12),
           risk.table.font = c(5),
           title = "Metastasis-free survival in UPS",
           legend.title = "",
           censor.size = 6,
           xlim = c(0, 90), # use 90 months as max
           break.x.by = 15, # use 15 months as breaks
           palette = c("#E41A1C", "#377EB8"),
)+
  xlab(label = "Time (months)")

#### Figure 3D ####
# Perform the metastasis-free survival analysis for the double positive macrophage groups in myxofibrosarcomas
fit <- survfit(Surv(MFS, status_MFS) ~ CD68_double_groups, 
               data = surv_macrophages[surv_macrophages$diagnosis == "MFS",])

# Check the summary of the survival analysis
fit

# Plot the Kaplan-Meier curve
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("CD68+CD163+ high", "CD68+CD163+ low"),
           font.legend = c(12),
           risk.table.font = c(5),
           title = "Metastasis-free survival in MFS",
           legend.title = "",
           censor.size = 6,
           xlim = c(0, 90), # use 90 months as max
           break.x.by = 15, # use 15 months as breaks
           palette = c("#E41A1C", "#377EB8"),
)+
  xlab(label = "Time (months)")

# Perform the metastasis-free survival analysis for the double positive macrophage groups in myxofibrosarcomas
fit <- survfit(Surv(MFS, status_MFS) ~ CD68_double_groups, 
               data = surv_macrophages[surv_macrophages$diagnosis == "UPS",])

# Check the summary of the survival analysis
fit

# Plot the Kaplan-Meier curve
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("CD68+CD163+ high", "CD68+CD163+ low"),
           font.legend = c(12),
           risk.table.font = c(5),
           title = "Metastasis-free survival in UPS",
           legend.title = "",
           censor.size = 6,
           xlim = c(0, 90), # use 90 months as max
           break.x.by = 15, # use 15 months as breaks
           palette = c("#E41A1C", "#377EB8"),
)+
  xlab(label = "Time (months)")
