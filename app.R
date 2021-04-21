#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)

# Define UI for input of data:
#ui <- fluidPage(
#    textInput("n_group1", "sd_group1", "n_group2","sd_group2"),
#    verbatimTextOutput("value")
#)


# Define UI for application that draws a histogram after simulation
ui <- fluidPage(
    # Application title
    titlePanel("Fly Survival Graphs - Customizable"),
    
    # Sidebar with a slider input for number of n per group
    sidebarLayout(
        sidebarPanel(
            # Input: Slider for the number of bins ----
            sliderInput(inputId = "bins",
                        label = "Number of bins:",
                        min = 5,
                        max = 50,
                        value = 30)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'orange')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
