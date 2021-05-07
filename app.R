# Anderson's Survival Analysis Script, but Shiny
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)
library(dplyr)
library(survminer)
library(survival)
library(readr)
library(readxl)
library(tidyverse)
library(cellranger)
library(ggplot2)
library(listviewer)
library(reactR)


#choices <- "choice 1, some other choice, and another one"

#shinyApp(
#    fluidPage(selectInput("id", "select", strsplit(choices, ", ")[[1]])),
#    function(...){}
#)

# Define UI for application
ui <- fluidPage(
    # Application title
    titlePanel("Fly Survival Graphs - Customizable"),
    # Sidebar with a slider input for number of n per group
    sidebarLayout(
        sidebarPanel(
            fileInput('survival_data_input', 'Select data file to upload',
                      accept = c(
                          '.xlsx'
                      )
            ),
            fileInput('survival_key_input', 'Select key file to upload',
                      accept = c(
                          '.xlsx'
                      )
            ),
            actionButton("goButton", "Go!", class = "btn-success"),
            tags$hr(),
            p('If you want a sample file to upload,',
              'you can first download the sample data and key files', 
              'and then try uploading them',
              a(href = 'https://github.com/AndersonButler/fly-survival/blob/9df0fdfe61b8e16df8c8c1243f40c931ea05492f/rawdata_template_survival.xlsx', 'data_template.xlsx'), 'or',
              a(href = 'https://github.com/AndersonButler/fly-survival/blob/9df0fdfe61b8e16df8c8c1243f40c931ea05492f/key_template_survival.xlsx', 'key_template_survival.xlsx'),
              'files, and then try uploading them.'),
        ),
        # Show a plot of the generated distribution
        mainPanel(
            tableOutput("raw_table"),
            plotOutput("SurvivalPlotTreatment"),
            plotOutput("SurvivalPlotGenotype"),
            reactjsonOutput( "rjed" ),
            #plotOutput("PairwiseComparisonsTreatment"),
            #plotOutput("PairwiseComparisonsGenotype")
        )
        )
)

# Define server logic
server <- function(input, output) {
    
    # Load the raw data into memory when file selected
    ## -------------------------------------------------------------------------------------------------
    
    # After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath' will contain the local filenames where the data can
    # be found
    
    observeEvent(input$goButton, {
        
        # Load data
        survival_data_path = input$survival_data_input$datapath
        survival_key_path = input$survival_key_input$datapath
        
        survival_data_raw <- read_excel(input$survival_data_input$datapath, range=cellranger::cell_rows(c(2,NA)))
        survival_key_raw <- read_excel(input$survival_key_input$datapath)
        
        # Print table of raw data
        output$raw_table = renderTable(survival_key_raw)
        
        # Math (see DavisLabSurvival_04052021)
        df.atrisk = data.frame(survival_data_raw$Vial,survival_data_raw$initial_n_in_vial)
        colnames(df.atrisk) = c("Vial","initial_n_in_vial")
        
        df.survivors = data.frame(survival_data_raw$Vial,survival_data_raw$survivors_in_vial,survival_data_raw$age_at_expt_end)
        colnames(df.survivors) = c("Vial","survivors_in_vial","age_at_expt_end")
        
        df.ndead = survival_data_raw # Break table into 3 dfs based on data type, with corresponding coordinates
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
        
        df.final = data.frame(Vial = numeric(0), Days = numeric(0),Dead1Excluded0 = numeric(0)) # Initialize temp dataframes
        df.temp = data.frame(Vial = numeric(0), Days = numeric(0),Dead1Excluded0 = numeric(0))
        
        #Loops through and generate individual entries for:
            # Each dead fly
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
        
            # Each left-censored fly
        
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
        
            # Each right-censored fly
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
        
        # Survival Math: Treatment
        fit = surv_fit(Surv(Days,Dead1Excluded0) ~ Treatment, data = df.final)
        
        legend_labels = as.factor(unique(df.final$Treatment))
        
        output$SurvivalPlotTreatment = renderPlot(ggsurvplot(fit, data = df.final,
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
            ))
        
        # Survival Math: Genotype
        fit2 = surv_fit(Surv(Days,Dead1Excluded0) ~ Genotype, data = df.final)
        
        legend_labels2 = as.factor(unique(df.final$Genotype))
        
        output$SurvivalPlotGenotype = renderPlot(ggsurvplot(fit2, data = df.final,
                   title = "Survival Curves: Genotype",
                   conf.int = TRUE,
                   pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
                   surv.median.line = "hv",            # Add median survival lines
                   legend.title = "GENOTYPE",               # Change legend titles
                   legend.labs = legend_labels2,  # Change legend labels
                   palette = "Pastel2",                    # Use JCO journal color palette
                   risk.table = F,                  # Add No at risk table
                   cumevents = F,                   # Add cumulative No of events table
                   tables.height = 0.15,               # Specify tables height
                   tables.theme = theme_cleantable(),  # Clean theme for tables
                   tables.y.text = FALSE               # Hide tables y axis text
        ))
        
        # Survival Math: Pairwise Comparisons for More Than Two Treatments
        res <- pairwise_survdiff(Surv(Days, Dead1Excluded0) ~ Treatment,
                                 data = df.final)
        
        summary = capture.output(res)
        pairwise_summary_plot_treatment = plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
        pairwise_summary_plot_treatment = pairwise_summary_plot_treatment + for (i in seq_along(summary)) {
                text(0, 100 - i*4, pos=4, summ[i], cex = 0.5, family='mono')
            }
        
        output$PairwiseComparisonsTreatment = renderPlot(pairwise_summary_plot_treatment)
 
        

        res2 <- pairwise_survdiff(Surv(Days, Dead1Excluded0) ~ Genotype,
                                 data = df.final)
        summary2= capture.output(res2)
        
        pairwise_summary_plot_genotype = plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=c(0, 100), ylim=c(0, 100), xaxs = 'i', yaxs = 'i')
        pairwise_summary_plot_genotype2 = pairwise_summary_plot_genotype + for (i in seq_along(summary2)) {
            text(0, 100 - i*4, pos=4, summ[i], cex = 0.5, family='mono')
        }
        
        output$PairwiseComparisonsGenotype = renderPlot(pairwise_summary_plot_genotype2)
        

        output$rjed <- renderReactjson({
                reactjson(as.list(summary2))
           })

        
            
        
        
        
        
    })
    


    
}

# Run the application 
shinyApp(ui = ui, server = server)
