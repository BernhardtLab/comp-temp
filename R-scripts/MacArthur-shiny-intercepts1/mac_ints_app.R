# a new go at creating a shiny app that slides parameter values in the MacArthur C-R model to view their effect on Chesson parameters
# this script is loosely based on JRB's macarthur shiny app, but I'm starting again from scratch here

# script DOB: 9 April 2025
# script author: Kaleigh Davis, postdoc at Guelph

library(shiny)
library(ggplot2)

# load macarthur source function
source("R-scripts/temp-indep-macarthur-KD.R")

ui <- 
  fluidPage(
  titlePanel("MacArthur to Chesson Expression Evaluator"), #title
  fluidRow(plotOutput("mctplot"), width = 12), #plot

  # # Add custom CSS for slider height
  # tags$style(HTML("
  #   .slider-height {
  #     height: 10px;  /* Adjust height here */
  #   }
  # ")),
  
  fluidRow(column(h5("Consumption rate intercepts"), #sliders for model inputs (parameters)
                  # sidebarLayout(
                  #   sidebarPanel(
                  # Input widgets for parameters
                  #Can I add a title here to declare that these are model intercepts at the reference temperature?
                  sliderInput("c1N_b", "c1N_b", min = 0, max = 5, value = 0.2, step = 0.05), #0.2;  max = 2 for all, 
                  # class = "slider-height" to implement slider height - not working, and not prioritizing this fix right now
                  sliderInput("c1P_b", "c1P_b", min = 0, max = 5, value = 0.4, step = 0.05), #0.4
                  sliderInput("c2N_b", "c2N_b", min = 0, max = 5, value = 0.2, step = 0.05), #0.4
                  sliderInput("c2P_b", "c2P_b", min = 0, max = 5, value = 0.4, step = 0.05), #0.2
                  width = 3), #between 1 and 12
  
           column(h5("Resource r and K intercepts"),
                  sliderInput("r_N_b", "r_N_b", min = 0, max = 5, value = 0.05, step = 0.05), #0.1, max = 5 for all
                  sliderInput("r_P_b", "r_P_b", min = 0, max = 5, value = 0.05, step = 0.05), #0.05
                  sliderInput("K_N_b", "K_N_b", min = 0, max = 10000, value = 2000, step = 500), #2000
                  sliderInput("K_P_b", "K_P_b", min = 0, max = 10000, value = 2000, step = 500), #2000
                  width = 3),
           
           column(h5("Conversion efficiency intercepts"),
                  sliderInput("v1N_b", "v1N_b", min = 0, max = 5, value = 0.2, step = 0.01), #0.2, max = 1 for all
                  sliderInput("v1P_b", "v1P_b", min = 0, max = 5, value = 0.4, step = 0.01), #0.4
                  sliderInput("v2N_b", "v2N_b", min = 0, max = 5, value = 0.2, step = 0.01), #0.4
                  sliderInput("v2P_b", "v2P_b", min = 0, max = 5, value = 0.4, step = 0.01), #0.2
                  width = 3),
           
           column(h5("Consumer mortality rate intercepts"),
                  sliderInput("m1_b", "m1_b", min = 0, max = 5, value = 0.01, step = 0.01), #0.01, max = 5 for all
                  sliderInput("m2_b", "m2_b", min = 0, max = 5, value = 0.01, step = 0.01), #0.01
                  width = 3)
  ),
    
    mainPanel(
    # Output display for results
    verbatimTextOutput("result")
    ),
)

server <- function(input, output) {
  # Reactive expression to perform calculations
  computed_result <- reactive({
    
    # Get parameters from input
    c1N_b <- input$c1N_b
    c1P_b <- input$c1P_b
    c2N_b <- input$c2N_b
    c2P_b <- input$c2P_b
    r_N_b <- input$r_N_b
    r_P_b <- input$r_P_b
    K_N_b <- input$K_N_b
    K_P_b <- input$K_P_b
    v1N_b <- input$v1N_b
    v1P_b <- input$v1P_b
    v2N_b <- input$v2N_b
    v2P_b <- input$v2P_b
    m1_b <- input$m1_b
    m2_b <- input$m2_b
    
    # Perform your complex mathematical function here
    result <- temp_indep_mac(c1N_b = c1N_b, c1P_b = c1P_b, #consumption rate of N and P at ref temp for species 1
                             c2N_b = c2N_b, c2P_b = c2P_b, #consumption rate of N and P at ref temp for species 2
                             r_N_b = r_N_b, r_P_b = r_P_b, #growth rate for each resource at ref temp
                             K_N_b = K_N_b, K_P_b = K_P_b, #carrying capacity for each resource at ref temp
                             v1N_b = v1N_b, v1P_b = v1P_b, #conversion efficiency for each resource at ref temp for species 1
                             v2N_b = v2N_b, v2P_b = v2P_b, #conversion efficiency for each resource at ref temp for species 2
                             m1_b = m1_b, m2_b = m2_b)
    
    # Return the result
    result
  })
  
  # Output the result for display
  output$result <- renderPrint({
    computed_result()
  })
  
  
# Render the plot of new_stabil_potential vs new_fit_ratio
  output$mctplot <- renderPlot({
    result <- computed_result()
  
    # Create a data frame for ggplot
    plot_data <- data.frame(
      new_stabil_potential = result$new_stabil_potential,
      new_fit_ratio = result$new_fit_ratio
    )
    
    # Generate the plot with ggplot2
    ggplot(plot_data, aes(x = new_stabil_potential, y =  new_fit_ratio)) +
      geom_ribbon(
        # data = data.frame(x = seq(min(plot_data$new_stabil_potential)*0.99, max(plot_data$new_stabil_potential)*1.01, 0.001)), #can't get this to work
        data = data.frame(x = seq(0, 2, 0.001)),
                                aes(x = x,
                                    y = NULL,
                                    ymin = exp(-x),
                                    ymax = 1/(exp(-x))),
                                fill = "grey", color = "black", alpha = 0.2) +
      geom_hline(yintercept = 1, linetype = 5) +
      coord_cartesian(xlim = c(0,2), ylim = c(-2, 5)) +
      geom_point(color = "blue", size = 6) +  # Blue points for the scatter plot
      labs(x = "-log(rho)", y = "log(k2/k1)",
           title = "Plot of New Stabilizing Potential vs New Fit Ratio") +
      theme_minimal() 
  })
}


# Run the application
shinyApp(ui = ui, server = server)
