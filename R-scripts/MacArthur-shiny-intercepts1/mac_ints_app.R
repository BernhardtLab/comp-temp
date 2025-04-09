# a new go at creating a shiny app that slides parameter values in the MacArthur C-R model to view their effect on Chesson parameters
# this script is loosely based on JRB's macarthur shiny app, but I'm starting again from scratch here

# script DOB: 9 April 2025
# script author: Kaleigh Davis, postdoc at Guelph

library(shiny)

# load macarthur source function
source("R-scripts/temp-indep-macarthur-KD.R")

ui <- 
  fluidPage(
  titlePanel("MacArthur to Chesson Expression Evaluator"),
  
  sidebarLayout(
    sidebarPanel(
      # Input widgets for parameters
      #Can I add a title here to declare that these are model intercepts at the reference temperature?
      sliderInput("param1", "C1N_b", min = 0, max = 10, value = 5),
      sliderInput("param2", "Parameter 2:", min = 0, max = 10, value = 5),
      sliderInput("param3", "Parameter 3:", min = 0, max = 10, value = 5),
      sliderInput("param4", "Parameter 4:", min = 0, max = 10, value = 5),
      sliderInput("param5", "Parameter 5:", min = 0, max = 10, value = 5),
      sliderInput("param6", "Parameter 6:", min = 0, max = 10, value = 5),
      sliderInput("param7", "Parameter 7:", min = 0, max = 10, value = 5),
      sliderInput("param8", "Parameter 8:", min = 0, max = 10, value = 5),
      # numericInput("param2", "Parameter 2:", value = 1),
      # Add more inputs as needed
    ),
    
    mainPanel(
      # Output display for results
      verbatimTextOutput("result")
    )
  )
)

server <- function(input, output) {
  # Reactive expression to perform calculations
  computed_result <- reactive({
    # Get parameters from input
    param1 <- input$param1
    param2 <- input$param2
    
    # Perform your complex mathematical function here
    result <- your_complex_function(param1, param2)
    
    # Return the result
    result
  })
  
  # Output the result for display
  output$result <- renderPrint({
    computed_result()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
