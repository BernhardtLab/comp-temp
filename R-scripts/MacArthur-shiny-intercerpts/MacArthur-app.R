#
#  You can run the application by clicking
# the 'Run App' button above (green arrow button).

### JB needs to fix the zero growth rate issue!

library(shiny)
library(cowplot)
library(tidyverse)
library(patchwork)
theme_set(theme_cowplot())

source("arrhenius_function.R")
source("temp_dependences_MacArthur.R")
source("plot_MacArthur.R")
source("plot_MacArthur_alpha.R")


params <- c("rN", "rP", "KN", "KP", "c1N", "c1P", "c2N", "c2P", "m1", "m2")
# Define UI for application that draws the abundance graph
ui <- fluidPage(
   
   # Application title
   titlePanel("MacArthur consumer-resource model"),
   fluidRow(column(("Code written by KD (building on JRB & PJK's), any mistakes are made by Kaleigh!"), width = 4)),
  #fluidRow(column(img(src='joeys-macarthur-equations.png'), width = 6)),
   # fluidRow(column(("Two consumer species (C1, C2) compete for two resources, (N, P). 
   # 				 We model the temperature dependence of parameters r, K, c, m and v using an Arrhenius function."), width = 6)),
   		fluidRow(plotOutput("coolplot"), width = 12),
   		fluidRow(
      	column(h4("For panel a, which parameters to dispay?"), offset = 0.1, 
      		selectInput('parameter1', 'Parameter 1', params),
      		selectInput('parameter2', 'Parameter 2', params, selected = params[[2]]),
      		selectInput('parameter3', 'Parameter 3', params, selected = params[[3]]),
      		selectInput('parameter4', 'Parameter 4', params, selected = params[[4]]), width = 2),
      	
      	
      	column(h4("Temp dependence of resource r & K"), offset = .2, 
         sliderInput("r_EaN",
                     "Ea of N's growth rate (rN)", min = 0, max = 1, value = 0.5, step= 0.05),
         sliderInput("r_EaP",
         			"Ea of P's growth rate (rP)", min = 0, max = 1, value = 0.5, step= 0.05),
         sliderInput("K_EaP",
         			"Ea of P's carrying capacity (KP)", min = -1, max = 1, value = -0.3, step= 0.05),
         sliderInput("K_EaN",
         			"Ea of N's carrying capacity (KN)", min = -1, max = 1, value = -0.3, step= 0.05), 
         #KD added this in just to see if it would work -- it does
         # sliderInput("v_EaN",
         #             "Ea of N's conversion efficiency (v)", min = -0.5, max = 0.5, value = 0, step= 0.05), 
               width = 2),
         column(h4("Baseline consumption rates"), offset = 0.2,
         	   sliderInput("c1N_b",
         	   			"C1's consumption rate of N (c1N)", min = 0, max = 1, value = 0.2, step= 0.05),
         	   sliderInput("c1P_b",
         	   			"C1's consumption rate of P (c1P)", min = 0, max = 1, value = 0.4, step= 0.05), 
         	   sliderInput("c2N_b",
         	   			"C2's consumption rate of N (c2N)", min = 0, max = 1, value = 0.4, step= 0.05),
         	   sliderInput("c2P_b",
         	   			"C2's consumption rate of P (c2P)", min = 0, max = 1, value = 0.2, step= 0.05),
         	   width = 2),
         
         column(h4("Temp dependence of consumption rates"), offset = 0.5,
         sliderInput("c_Ea1N",
         			"Ea of C1's consumption of N (c1N)", min = -1, max = 1, value = 0, step= 0.05),
         sliderInput("c_Ea1P",
         			"Ea of C1's consumption of P (c1P)", min = -1, max = 1, value = 0, step= 0.05),
         sliderInput("c_Ea2N",
         			"Ea of C2's consumption of N (c2N)", min = -1, max = 1, value = 0, step= 0.05),
         sliderInput("c_Ea2P",
         			"Ea of C1's consumption of P (c2P)", min = -1, max = 1, value = 0, step= 0.05), 
         width = 2),
         column(h4("Temp dependence of consumer mortality"), offset = 0.5,
         	   sliderInput("m_Ea1",
         	   			"Ea of C1's mortality rate (m1)", min = 0, max = 1, value = 0, step= 0.05),
         	   sliderInput("m_Ea2",
         	   			"Ea of C2's mortality rate (m2)", min = 0, max = 1, value = 0, step= 0.05), 
         	   width = 2),
         column(h4("Conversion efficiencies"), offset = 0.5,
         	   sliderInput("v1N_b",
         	   			"Conversion of N into C1", min = 0.1, max = 1, value = 0.2, step= 0.05),
         	   sliderInput("v2N_b",
         	   			"Conversion of N into C2", min = 0.1, max = 1, value = 0.4, step= 0.05),
         	   sliderInput("v1P_b",
         	   			"Conversion of P into C1", min = 0.1, max = 1, value = 0.4, step= 0.05),
         	   sliderInput("v2P_b",
         	   			"Conversion of P into C2", min = 0.1, max = 1, value = 0.2, step= 0.05),
         	   width = 2)
         
      # ),

      
      # Show a plot of the time series
      # mainPanel(
      #    plotOutput("coolplot"),
      #    width = 0.8
      # )
   )
   )

server <- function(input, output) {
	
	
   output$coolplot <- renderPlot({
   	
   	
 	Data.temperature.r = temp_dependences_MacArthur(T = seq(25, 25.1, by = 0.1), 
   													r_EaN = input$r_EaN, r_EaP = input$r_EaP, 
   													K_EaN = input$K_EaN, K_EaP = input$K_EaP, 
   													c_Ea1N = input$c_Ea1N, c_Ea1P = input$c_Ea1P, 
   													c_Ea2N = input$c_Ea2N, c_Ea2P = input$c_Ea2P, 
   													v_EaN = 0.0, v_EaP = 0.0, 
 													v1N_b = input$v1N_b, v2N_b = input$v2N_b,
 													v1P_b = input$v1P_b, v2P_b = input$v2P_b,
   													m_Ea1 = input$m_Ea1, m_Ea2 = input$m_Ea2,
 													c1N_b = input$c1N_b, c2P_b = input$c2P_b,
 													c1P_b = input$c1P_b, c2N_b = input$c2N_b)
 	
   	Data.parameter.r = Data.temperature.r[, c("T", input$parameter1, input$parameter2, input$parameter3, input$parameter4)] %>% 
   	  gather(value=value, key=parameter, -T)
   	
   	Plot.r = plot_MacArthur(Data.parameter.r, Data.temperature.r)
   	
Plot.r
   	
     
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

