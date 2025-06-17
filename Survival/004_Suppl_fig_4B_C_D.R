### Script for visualizing the other survival curves (disease-specific in UPS & metastasis-free in UPS and myxofibrosarcoma)
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

#### Supplementary figure 4B ####
# Perform the disease-specific survival analysis for the T cell groups in UPS
fit <- survfit(Surv(OS, status_DSS) ~ Tcell_groups, 
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
           title = "Disease-specific survival in UPS",
           legend.title = "",
           censor.size = 6,
           xlim = c(0, 90), # use 90 months as max
           break.x.by = 15, # use 15 months as breaks
           palette = c("#E41A1C", "#377EB8"),
)+
  xlab(label = "Time (months)")

# Perform the disease-specific survival analysis for the double positive macrophage groups in UPS
fit <- survfit(Surv(OS, status_DSS) ~ CD68_double_groups, 
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
           title = "Disease-specific survival in UPS",
           legend.title = "",
           censor.size = 6,
           xlim = c(0, 90), # use 90 months as max
           break.x.by = 15, # use 15 months as breaks
           palette = c("#E41A1C", "#377EB8"),
)+
  xlab(label = "Time (months)")

#### Supplementary figure 4C ####
# Create a factor to order the groups from high to low
surv_macrophages$combined_groups <- factor(surv_macrophages$combined_groups,
                                           levels = c("Both high", 
                                                      "High and low",
                                                      "Both low"))

# Perform the metastasis-free survival analysis for the combined T cell and macrophage groups in UPS
fit <- survfit(Surv(MFS, status_MFS) ~ combined_groups, 
               data = surv_macrophages[surv_macrophages$diagnosis == "UPS",])

# Check the summary of the survival analysis
fit

# Plot the Kaplan-Meier curve
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("Both high", "High and low", "Both low"),
           font.legend = c(12),
           risk.table.font = c(5),
           title = "Metastasis-free survival in UPS",
           legend.title = "",
           censor.size = 6,
           xlim = c(0, 90), # use 90 months as max
           break.x.by = 15, # use 15 months as breaks
           palette = c("#E41A1C", "#4DAF4A", "#377EB8"),
)+
  xlab(label = "Time (months)")

#### Supplementary figure 4D ####
# Perform the metastasis-free survival analysis for the single positive macrophage groups in myxofibrosarcomas
fit <- survfit(Surv(MFS, status_MFS) ~ CD68_single_groups, 
               data = surv_macrophages[surv_macrophages$diagnosis == "MFS",])

# Check the summary of the survival analysis
fit

# Plot the Kaplan-Meier curve
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("CD68+ high", "CD68+ low"),
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

# Perform the metastasis-free survival analysis for the single positive macrophage groups in UPS
fit <- survfit(Surv(MFS, status_MFS) ~ CD68_single_groups, 
               data = surv_macrophages[surv_macrophages$diagnosis == "UPS",])

# Check the summary of the survival analysis
fit

# Plot the Kaplan-Meier curve
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("CD68+ high", "CD68+ low"),
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
