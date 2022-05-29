
#title: "Analysis script for survival data in Drosophila"
#author: "Anderson Butler"
#date: "3/10/2021"
#most_recent_update: "4/02/2021"
 
## -------------------------------------------------------------------------------------------------
library(tcltk)
survival_excel_path = tk_choose.files(caption = "Select survival spreadsheet")
survival_key_path = tk_choose.files(caption = "Select key spreadsheet")
set_export_directory = tk_choose.dir(caption = "Select export directory")

setwd(set_export_directory)
run_day = as.character(Sys.Date())


## -------------------------------------------------------------------------------------------------
library(readr)
library(readxl)

# Read excel spreadsheets
survival_data_raw <- read_excel(survival_excel_path, range=cellranger::cell_rows(c(2,NA)))
survival_key_raw <- read_excel(survival_key_path)


## -------------------------------------------------------------------------------------------------
library(tidyverse)

df.atrisk = data.frame(survival_data_raw$Vial,survival_data_raw$initial_n_in_vial)
colnames(df.atrisk) = c("Vial","initial_n_in_vial")

df.survivors = data.frame(survival_data_raw$Vial,survival_data_raw$survivors_in_vial,survival_data_raw$age_at_expt_end)
colnames(df.survivors) = c("Vial","survivors_in_vial","age_at_expt_end")

df.ndead = survival_data_raw
df.ndead = df.ndead[, -grep("age", colnames(df.ndead))]
df.ndead = df.ndead[, -grep("nexcluded", colnames(df.ndead))]
df.ndead = df.ndead[, -grep("survivors", colnames(df.ndead))]
df.ndead = df.ndead[, -grep("initial", colnames(df.ndead))]

df.nexcluded = survival_data_raw
df.nexcluded = df.nexcluded[, -grep("age", colnames(df.nexcluded))]
df.nexcluded = df.nexcluded[, -grep("ndead", colnames(df.nexcluded))]
df.nexcluded = df.nexcluded[, -grep("survivors", colnames(df.nexcluded))]
df.nexcluded = df.nexcluded[, -grep("initial", colnames(df.nexcluded))]

df.age = survival_data_raw
df.age = df.age[, -grep("ndead", colnames(df.age))]
df.age = df.age[, -grep("nexcluded", colnames(df.age))]
df.age = df.age[, -grep("survivors", colnames(df.age))]
df.age = df.age[, -grep("initial", colnames(df.age))]



## -------------------------------------------------------------------------------------------------
df.final = data.frame(Vial = numeric(0), Days = numeric(0),Dead1Excluded0 = numeric(0))


## -------------------------------------------------------------------------------------------------
# Initialize a temporary dataframe to store the rows for each vial
df.temp = data.frame(Vial = numeric(0), Days = numeric(0),Dead1Excluded0 = numeric(0))

#Loop to generate individual entries for each dead fly
for (x in 1:length(survival_data_raw$Vial)) {
  for (y in 4:ncol(df.ndead)) {
    if (df.ndead[x,y] > 0) {
      myvector1 = c(Vial = df.ndead$Vial[x], Genotype = survival_key_raw$Genotype[x], Treatment = survival_key_raw$Treatment[x], Days = (df.age[x,y+1]), Dead1Excluded0 = 1) #vector of dead fly info
      df.temp = as.data.frame(myvector1) # vector as a df
      colnames(df.temp)[4] = "Days"
      ndead = as.numeric(df.ndead[x,y])
      df.temp = df.temp[rep(seq_len(nrow(df.temp)), each = ndead), ] # Repeat the row until you have 1 row for each dead fly (in that vial, at that timepoint)
      df.final = rbind(df.final,df.temp)
    }
  }
}

#Loop to generate individual entries for each left-censored fly

for (x in 1:length(survival_data_raw$Vial)) {
  for (y in 4:ncol(df.ndead)) {
    if (df.nexcluded[x,y] > 0){
      myvector2 = c(Vial = df.nexcluded$Vial[x], Genotype = survival_key_raw$Genotype[x], Treatment = survival_key_raw$Treatment[x], Days = (df.age[x,y+1]), Dead1Excluded0 = 0) #vector of left-censored fly info
      df.temp = as.data.frame(myvector2) # vector as a df
      colnames(df.temp)[4] = "Days"
      nexcluded = as.numeric(df.nexcluded[x,y])
      df.temp = df.temp[rep(seq_len(nrow(df.temp)), each = nexcluded), ] # Repeat the row until you have 1 row for each dead fly (in that vial, at that timepoint)
      df.final = rbind(df.final,df.temp)
    }
  }
}

#Loop to generate individual entries for each right-censored fly
for (x in 1:length(survival_data_raw$Vial)) {
  if (colnames(df.survivors[2]) == "survivors_in_vial") {
  if (is.na(df.survivors$survivors_in_vial[x])=="FALSE") {
    if (df.survivors[x,2] > 0) {
      myvector3 = as.list(c(Vial = df.survivors$Vial[x], Genotype = survival_key_raw$Genotype[x], Treatment = survival_key_raw$Treatment[x], Days = as.double(df.survivors$age_at_expt_end[x]), Dead1Excluded0 = 0)) #right-censored
      df.temp = as.data.frame(myvector3) # vector as a df
      df.temp$Days = as.numeric(df.temp$Days)
      df.temp$Dead1Excluded0 = as.numeric(df.temp$Dead1Excluded0)
      n_right_censored = as.numeric(df.survivors[x,2])
      df.temp = df.temp[rep(seq_len(nrow(df.temp)), each = n_right_censored), ]
      df.final = rbind(df.final,df.temp)
      }
  } else {
    vial_name = df.atrisk$Vial[x]
    dead_or_censored_so_far = nrow(df.final[df.final$Vial == vial_name,])
    survivors_so_far = df.atrisk$initial_n_in_vial[df.atrisk$Vial==vial_name] - dead_or_censored_so_far # select from df.at
    if (survivors_so_far > 0) {
      most_recent_age_col = ncol(df.age) # for calculating age below
      age_so_far = as.numeric(df.age[x,most_recent_age_col])
      myvector4 = as.list(c(Vial = df.survivors$Vial[x], Genotype = survival_key_raw$Genotype[x], Treatment = survival_key_raw$Treatment[x], Days = age_so_far, Dead1Excluded0 = 0)) #right-censored
      df.temp = as.data.frame(myvector4)
      df.temp$Days = as.numeric(df.temp$Days)
      df.temp$Dead1Excluded0 = as.numeric(df.temp$Dead1Excluded0)
      n_right_censored = as.numeric(survivors_so_far)
      df.temp = df.temp[rep(seq_len(nrow(df.temp)), each = n_right_censored), ]
      df.final = rbind(df.final,df.temp)
      rm(dead_or_censored_so_far,survivors_so_far,most_recent_age_col,age_so_far)
    }
  }
  }
}

my_file_name = paste("SurvivalAnalysisTable_",run_day, sep="")
my_file_name = paste(my_file_name,".csv", sep="")
write.csv(df.final,my_file_name)

# First Survival Curve
## -------------------------------------------------------------------------------------------------
library(dplyr)
library(survminer)
library(survival)
library(ggokabeito)
library(wesanderson)

fit = surv_fit(Surv(Days,Dead1Excluded0) ~ Treatment, data = df.final)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
#scales::show_col(safe_colorblind_palette)

okabe_ito_palette = ggokabeito::palette_okabe_ito()
#scales::show_col(okabe_ito_palette)

#my_palette = wes_palette("Darjeeling1")

library(ggplot2)
ggsurvplot(fit, data = df.final,
           title = "Survival Curves: Treatment",
           legend = "bottom",
           conf.int = TRUE,
           pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "TREATMENT",               # Change legend titles
           #legend.labs = legend_labels,  # Change legend labels
           palette = okabe_ito_palette,                    # Use JCO journal color palette
           risk.table = T,                  # Add No at risk table
           cumevents = T,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
) + guides(colour = guide_legend(nrow = 4))



## -------------------------------------------------------------------------------------------------
library(dplyr)
library(survminer)
library(survival)

fit = surv_fit(Surv(Days,Dead1Excluded0) ~ Genotype, data = df.final)

legend_labels = as.factor(unique(df.final$Genotype))

ggsurvplot(fit, data = df.final,
           title = "Survival Curves: Genotype",
           conf.int = TRUE,
           pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "GENOTYPE",               # Change legend titles
           #legend.labs = legend_labels,  # Change legend labels
           palette = okabe_ito_palette,                    # Use JCO journal color palette
           risk.table = F,                  # Add No at risk table
           cumevents = F,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
)+ guides(colour = guide_legend(nrow = 4))

## -------------------------------------------------------------------------------------------------
df.final$combined_treatment_group = paste(df.final$Genotype,df.final$Treatment)

fit = surv_fit(Surv(Days,Dead1Excluded0) ~ combined_treatment_group, data = df.final)

legend_labels = as.factor(unique(df.final$combined_treatment_group))

ggsurvplot(fit, data = df.final,
           title = "Survival Curves: Genotype",
           conf.int = TRUE,
           pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "GENOTYPE",               # Change legend titles
           #legend.labs = legend_labels,  # Change legend labels
           palette = okabe_ito_palette,                    # Use JCO journal color palette
           risk.table = F,                  # Add No at risk table
           cumevents = F,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
)+ guides(colour = guide_legend(nrow = 4))

summary(fit)




## -------------------------------------------------------------------------------------------------
# Print the pairwise comparisons directly too:
res <- pairwise_survdiff(Surv(Days, Dead1Excluded0) ~ Treatment,
                         data = df.final, p.adjust.method="holm")
res

summ = capture.output(res)
plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
for (i in seq_along(summ)) {
  text(0, 100 - i*4, pos=4, summ[i], cex = 0.5, family='mono')
}

res <- pairwise_survdiff(Surv(Days, Dead1Excluded0) ~ Genotype,
                         data = df.final,  p.adjust.method="holm")
res

summ = capture.output(res)
plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
for (i in seq_along(summ)) {
  text(0, 100 - i*4, pos=4, summ[i], cex = 0.5, family='mono')
}


## -------------------------------------------------------------------------------------------------
#Make a 2nd table with the same headings for HEALTHSPAN TABLE *Ron's term - try to discourage use of term healthspan for survival, 
#in which we use only the data to the half-life (of the longest living group, drug or no-drug), to calculate the statistics.

# Fit a survival curve based on the entire data set, similar to above
df.final$combined_treatment_group = paste(df.final$Genotype,df.final$Treatment)
fit = surv_fit(Surv(Days,Dead1Excluded0) ~ combined_treatment_group, data = df.final)
legend_labels = as.factor(unique(df.final$combined_treatment_group))

# Find the longest median lifespan from the two treatment groups
surv_medians = surv_median(fit,combine=F)
max_median_life_span = max(surv_medians$median)

# Select only the data from up to the median life span (df.early_to_mid_life_survival) as is, and censor the remainder of the data
# Note that censoring increases statistical power by including information about number of flies under observation. This should not be omitted!!
df.early_to_mid_life_survival = df.final[df.final$Days<=max_median_life_span,]
df.censored_by_cutoff = df.final[df.final$Days>max_median_life_span,]
df.censored_by_cutoff$Dead1Excluded0 = 0 # all censored
df.censored_by_cutoff$Days = 31 #all censored on day 31
df.final2 = rbind(df.early_to_mid_life_survival,df.censored_by_cutoff)

# Fit a new survival object to the early life data, without censored flies from later in life (unused)
df.early_to_mid_life_survival$combined_treatment_group = paste(df.early_to_mid_life_survival$Genotype,df.early_to_mid_life_survival$Treatment)
fit = surv_fit(Surv(Days,Dead1Excluded0) ~ combined_treatment_group, data = df.early_to_mid_life_survival)
legend_labels = as.factor(unique(df.early_to_mid_life_survival$combined_treatment_group))

# Fit a new survival object to the early life data, with censored flies from later in life
df.final2$combined_treatment_group = paste(df.final2$Genotype,df.final2$Treatment)
fit = surv_fit(Surv(Days,Dead1Excluded0) ~ combined_treatment_group, data = df.final2)
legend_labels = as.factor(unique(df.final2$combined_treatment_group))

ggsurvplot(fit, data = df.final2,
           title = "Survival Curves: Genotype",
           subtitle = "Early to middle lifespan survival",
           conf.int = TRUE,
           pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "GENOTYPE",               # Change legend titles
           #legend.labs = legend_labels,  # Change legend labels
           palette = okabe_ito_palette,                    # Use JCO journal color palette
           risk.table = F,                  # Add No at risk table
           cumevents = F,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
) + guides(colour = guide_legend(nrow = 4))

# Pairwise comparisons by genotype
res <- pairwise_survdiff(Surv(Days, Dead1Excluded0) ~ Genotype,
                         data = df.final2, p.adjust.method="holm")
res

summ = capture.output(res)
plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
for (i in seq_along(summ)) {
  text(0, 100 - i*4, pos=4, summ[i], cex = 0.5, family='mono')
}

# Pairwise comparisons by drug
res <- pairwise_survdiff(Surv(Days, Dead1Excluded0) ~ Treatment,
                         data = df.final2, p.adjust.method="holm")
res

summ = capture.output(res)
plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
for (i in seq_along(summ)) {
  text(0, 100 - i*4, pos=4, summ[i], cex = 0.5, family='mono')
}


## -------------------------------------------------------------------------------------------------
# Data preparation and computing cox model

# It looks like the package picks whatever the first condition is for each column as the reference by default. This can be changed using factors, 
# should be coded to reflect the experimental design as treatments/genotypes change

# Begin updates from 3/2/2022
## I wrapped code to build a Cox model in a for loop, to keep it from breaking with experimental designs which do not include >1 Genotype and treatment
library(survival)
if (length(unique(df.final2$Genotype))>1 & length(unique(df.final2$Treatment))>1){
  res.cox <- coxph(Surv(Days, Dead1Excluded0) ~ Genotype + Treatment, data =  df.final2)
  # Plot the baseline survival function
  # with showing all individual predicted surv. curves
  #ggadjustedcurves(res.cox, data = df.final2,
  #                    individual.curves = TRUE)
  ggforest(res.cox)
}

## Calculating summary statistics for the combined treatment groups:
result_summary = surv_summary(fit, data = df.final2)
result_table = attr(result_summary, "table")

library(formattable)
formattable(result_table)



