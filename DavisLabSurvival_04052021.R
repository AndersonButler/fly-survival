
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
      myvector1 = c(Vial = df.ndead$Vial[x], Genotype = survival_key_raw$Genotype[x], Treatment = survival_key_raw$Treatment[x], Days = (df.age[x,y]), Dead1Excluded0 = 1) #vector of dead fly info
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
      myvector2 = c(Vial = df.nexcluded$Vial[x], Genotype = survival_key_raw$Genotype[x], Treatment = survival_key_raw$Treatment[x], Days = (df.age[x,y]), Dead1Excluded0 = 0) #vector of left-censored fly info
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
    dead_or_censored_so_far = nrow(df.final[df.final$Vial == x,])
    survivors_so_far = df.atrisk$initial_n_in_vial[x] - dead_or_censored_so_far
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

fit = surv_fit(Surv(Days,Dead1Excluded0) ~ Treatment, data = df.final)

legend_labels = as.factor(unique(df.final$Treatment))

ggsurvplot(fit, data = df.final,
           title = "Survival Curves: Treatment",
           conf.int = TRUE,
           pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
           surv.median.line = "hv",            # Add median survival lines
           legend.title = "TREATMENT",               # Change legend titles
           legend.labs = legend_labels,  # Change legend labels
           palette = "Pastel1",                    # Use JCO journal color palette
           risk.table = F,                  # Add No at risk table
           cumevents = F,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
)


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
           legend.labs = legend_labels,  # Change legend labels
           palette = "Pastel1",                    # Use JCO journal color palette
           risk.table = F,                  # Add No at risk table
           cumevents = F,                   # Add cumulative No of events table
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = FALSE               # Hide tables y axis text
)


## -------------------------------------------------------------------------------------------------
# And let's try to print the pairwise comparisons directly too:
res <- pairwise_survdiff(Surv(Days, Dead1Excluded0) ~ Treatment,
                         data = df.final)
res

summ = capture.output(res)
plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
for (i in seq_along(summ)) {
  text(0, 100 - i*4, pos=4, summ[i], cex = 0.5, family='mono')
}

res <- pairwise_survdiff(Surv(Days, Dead1Excluded0) ~ Genotype,
                         data = df.final)
res

summ = capture.output(res)
plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
for (i in seq_along(summ)) {
  text(0, 100 - i*4, pos=4, summ[i], cex = 0.5, family='mono')
}


## -------------------------------------------------------------------------------------------------
# Data preparation and computing cox model

## I think it picks whatever the first condition is for each column as the reference by default. Let me see if I can change this by using factors...
df.final2 = df.final
df.final2$Treatment = factor(df.final2$Treatment, levels = c("DMSO","DRUG"))

library(survival)
res.cox <- coxph(Surv(Days, Dead1Excluded0) ~ Genotype + Treatment, data =  df.final2)
# Plot the baseline survival function
# with showing all individual predicted surv. curves
#ggadjustedcurves(res.cox, data = df.final2,
#                    individual.curves = TRUE)

ggforest(res.cox)

