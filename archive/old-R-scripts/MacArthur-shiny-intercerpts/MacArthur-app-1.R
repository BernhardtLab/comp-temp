# KD April 2025
# adapted from JB's shiny app, saved in this repo in the folder MacArthur-shiny
# this script follows MacArthur-app in my folder MacArthur-shiny-intercerpts. There, I made small tweaks to start to understand the structure of the app. Here, I will begin to diverge in larger ways from the original app.

#  You can run the application by running the shiny command at the bottom of the script

library(shiny)
library(cowplot)
library(tidyverse)
library(patchwork)
theme_set(theme_cowplot())

source("R-scripts/MacArthur-shiny-intercerpts/arrhenius_function-KD.R") # this one is fine
# source("R-scripts/MacArthur-shiny/arrhenius_function.R")

source("R-scripts/MacArthur-shiny-intercerpts/temp_dependences_MacArthur-KD.R")
# source("R-scripts/temp-dep-macarthur-KD.R") #swapping Joey's function out for mine isn't working now, but this is ultimately what needs to be done in order to make all of the intercepts slide-able. Right now half are hard-coded in at the function stage and half are written in here. Many (all?) of the hard-coded ones are overwritten at time of deploying app.
# source("R-scripts/MacArthur-shiny/temp_dependences_MacArthur.R")

# source("R-scripts/MacArthur-shiny-intercerpts/temp-indep-macarthur-KD.R") #the plots will not run when I have this not commented out, so leaving it commented out for now

source("R-scripts/MacArthur-shiny-intercerpts/plot_MacArthur-KD.R") #this one is fine
# source("R-scripts/MacArthur-shiny/plot_MacArthur.R")

# source("R-scripts/MacArthur-shiny-intercerpts/plot_MacArthur_alpha.R") #not used for now

params <- c("rN", "rP", "KN", "KP", "c1N", "c1P", "c2N", "c2P", "m1", "m2")

# Define UI for application that draws the abundance graph ------------------------------------
#What I would like to do is remove all the bits of this that are temperature dependent, and feed in a new model that does not simulate over temperature, but rather just executes the equation with the selected parameters. 

#When I just comment out the EA sub-columns, that breaks the code. It produces the error "Error in data.frame: arguments imply differing number of rows: 2, 0".
#I tried trouble shooting this by manually setting all EA parameters in the temp_dependences_macarthur function (in server section) equal to 0, but this broke it in a different way.  "Error in UseMethod: no applicable method for 'gather' applied to an object of class "c('double', 'numeric')" "

ui <- fluidPage(
   
   # Application title
   titlePanel("MacArthur consumer-resource model"),
   fluidRow(column(("Code written by KD (building on JRB & PJK's), any mistakes are made by Kaleigh!"), width = 4)),
  #fluidRow(column(img(src='joeys-macarthur-equations.png'), width = 6)),
   fluidRow(column(("Two consumer species (C1, C2) compete for two resources, (N, P)"), width = 6)),
   		fluidRow(plotOutput("coolplot"), width = 4), #this displays the plots! I couldn't get this to dynamically update in the R app viewer. Had to close out of app and re-launch from the last line of code.
   		fluidRow(

      	column(h4("For panel a, which parameters to display?"), offset = 0.1, #h4 gives header. Not sure what offset is about.
      	       #this next bit gives the drop down options for each parameter
      	       #the inputId argument is referred to in the server object below, inside the Data.parameter.r object
      		selectInput(inputId = 'parameter1', label = 'Parameter 1', choices = params),
      		selectInput('parameter2', 'Parameter 2', params, selected = params[[2]]),
      		selectInput('parameter3', 'Parameter 3', params, selected = params[[3]]),
      		selectInput('parameter4', 'Parameter 4', params, selected = params[[4]]), width = 2),


      # 	column(h4("Temp dependence of resource r & K"), offset = .2,
      #    sliderInput("r_EaN",
      #                "Ea of N's growth rate (rN)", min = 0, max = 0, value = 0.5, step= 0.05),
      #    sliderInput("r_EaP",
      #    			"Ea of P's growth rate (rP)", min = 0, max = 0, value = 0, step= 0.05),
      #    sliderInput("K_EaP",
      #    			"Ea of P's carrying capacity (KP)", min = 0, max = 0, value = -0.3, step= 0.05),
      #    sliderInput("K_EaN",
      #    			"Ea of N's carrying capacity (KN)", min = 0, max = 0, value = -0.3, step= 0.05),
      #          width = 2),
      	
         column(h4("Baseline consumption rates (c)"), offset = 0.2,
         	   sliderInput("c1N_b",
         	   			"C1's consumption rate of N (c1N)", min = 0, max = 1, value = 0.2, step= 0.05),
         	   sliderInput("c1P_b",
         	   			"C1's consumption rate of P (c1P)", min = 0, max = 1, value = 0.4, step= 0.05), 
         	   sliderInput("c2N_b",
         	   			"C2's consumption rate of N (c2N)", min = 0, max = 1, value = 0.4, step= 0.05),
         	   sliderInput("c2P_b",
         	   			"C2's consumption rate of P (c2P)", min = 0, max = 1, value = 0.2, step= 0.05),
         	   width = 2),
         
         # column(h4("Temp dependence of consumption rates"), offset = 0.5,
         # sliderInput("c_Ea1N",
         # 			"Ea of C1's consumption of N (c1N)", min = -1, max = 1, value = 0, step= 0.05),
         # sliderInput("c_Ea1P",
         # 			"Ea of C1's consumption of P (c1P)", min = -1, max = 1, value = 0, step= 0.05),
         # sliderInput("c_Ea2N",
         # 			"Ea of C2's consumption of N (c2N)", min = -1, max = 1, value = 0, step= 0.05),
         # sliderInput("c_Ea2P",
         # 			"Ea of C1's consumption of P (c2P)", min = -1, max = 1, value = 0, step= 0.05),
         # width = 2),
      	
         # column(h4("Temp dependence of consumer mortality"), offset = 0.5,
         # 	   sliderInput("m_Ea1",
         # 	   			"Ea of C1's mortality rate (m1)", min = 0, max = 1, value = 0, step= 0.05),
         # 	   sliderInput("m_Ea2",
         # 	   			"Ea of C2's mortality rate (m2)", min = 0, max = 1, value = 0, step= 0.05),
         # 	   width = 2),
         
         column(h4("Baseline conversion efficiencies (v)"), offset = 0.5,
         	   sliderInput("v1N_b",
         	   			"Conversion of N into C1", min = 0.1, max = 1, value = 0.2, step= 0.05),
         	   sliderInput("v2N_b",
         	   			"Conversion of N into C2", min = 0.1, max = 1, value = 0.4, step= 0.05),
         	   sliderInput("v1P_b",
         	   			"Conversion of P into C1", min = 0.1, max = 1, value = 0.4, step= 0.05),
         	   sliderInput("v2P_b",
         	   			"Conversion of P into C2", min = 0.1, max = 1, value = 0.2, step= 0.05),
         	   width = 2)
         
   )
   )

server <- function(input, output) {
	
	
   output$coolplot <- renderPlot({
   	
   	
 	Data.temperature.r = temp_dependences_MacArthur(T = seq(25, 25.1, by = 0.1), 
   													r_EaN = 0, r_EaP = 0, #these ask for input from the sliders in the UI
   													K_EaN = 0, K_EaP = 0, 
   													c_Ea1N = 0, c_Ea1P = 0, 
   													c_Ea2N = 0, c_Ea2P = 0, 
   													v_EaN = 0.0, v_EaP = 0.0, 
 													v1N_b = input$v1N_b, v2N_b = input$v2N_b,
 													v1P_b = input$v1P_b, v2P_b = input$v2P_b,
   													m_Ea1 = 0, m_Ea2 = 0,
 													c1N_b = input$c1N_b, c2P_b = input$c2P_b,
 													c1P_b = input$c1P_b, c2N_b = input$c2N_b) #I need to figure out how to refer to the sliders here if I haven't pre-defined the input names in the first column with the dropdowns (since I no longer need that dropdown menu at all)
 	
   	Data.parameter.r = Data.temperature.r[, c("T", input$parameter1, input$parameter2, input$parameter3, input$parameter4)] %>% 
   	  gather(value=value, key=parameter, -T)
   
   	
   	Plot.r = plot_MacArthur(Data.parameter.r, Data.temperature.r) #successfully edited the plot_MacArthur function so that it only plots the coexistence space and not all the other stuff
   	
Plot.r
   	
     
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

