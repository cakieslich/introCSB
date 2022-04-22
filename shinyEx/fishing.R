library(deSolve)
library(shiny)

fishMod <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dNpf <- r*Npf*(K-Npf)/K - pfr*Npf
    dNcf <- r*Ncf*(K-Ncf)/K - cfr
    
    if((Npf+dNpf)<0) {dNpf <- -Npf} # Added to make sure populations are non-negative
    if((Ncf+dNcf)<0) {dNcf <- -Ncf} # Added to make sure populations are non-negative
    
    list(c(dNpf, dNcf))
  })
}

runModel <- function(input){
  ## Chaos in the atmosphere
  
  # All paramters in this example are controlled by sliders in the interface
  # As the slider is moved the reactive variable input is updated with new values
  parameters <- c(r = input$r, 
                  K = input$K,
                  pfr = input$pfr, 
                  cfr = input$cfr
                  )
  
  state <- c(Npf = 10, Ncf = 10)
  times <- seq(0, 200, by = 0.1)
  
  # Table containing data for events (perturbations) applied to system
  # In this example the same perturbattion is applied to both variables/populations
  # Multiple events can be defined in the same table
  eventdat <- data.frame(var = c("Npf","Ncf"), # Which variables affected
                         time = input$pTime, # When the even occurs (Controlled by slider)
                         value = input$perturb, # Value of change/perturbation (Controlled by slider)
                         method = "add") # Type of change, here the perturbation is added to the variable value
                                         # at the specified time.
  
  out <- ode(y = state, 
             times = times, 
             func = fishMod, 
             parms = parameters,
             events=list(data = eventdat) # The events argument is used to introduce changes to the system
                                          # at specified times
             )
  
  return(out)
}

# A shiny app has two main components: a user interface (ui) and a server
# The ui defines the layout of the application
# The server conects to interface to the underlying model/code

# ui is ultimately generated in HTML so some familiarity with HTML can help but is not necessary
# Many of the functions used to build ui are named after an HTML equivalent. 
ui <- fluidPage(
    h3("Population Models: Fishing example"), # Header level 3 (1 is largest and font decreases with increasing number)
    br(), # Line break, leaves an empty line
    
   h4("Effect of perturbations on system"), # Header level 4
   tags$ul( # Beginning of a bulleted list
     # li are each items in the bulleted list
     tags$li("In this example, we can explore the affect of different modes of fishing (proportional or 
             constant fishing rates) on a fish population."
     ),
     tags$li("The underlying model assumes the fish have the same growth rate and carrying capacity."
     ),
     tags$li("Fishing rates should be selected so that the two populations reach similar steady states."
     ),
     tags$li("The interface also allows for perturbations to be introduced, by specifying the time and size
              of the perturbation. Do the two modes of fishing affect the fish populations' resilliance to 
             perturbation? What is the maximum perturbation from which the fish population can recover 
             for each of the two modes of fishing?"
     )
   ),
   br(),
   # Begins a sidebar layout where there is a small sidebar on the left and a main (larger) panel on the right
   sidebarLayout(
     sidebarPanel( # Begins the sidebar panel
       tabsetPanel( # Begins a tabbed layoutt
         tabPanel("Parameters", # Starts tab panel for parameter definition
                  # Each sliderInput defines the a slider that controls an input to the simulation
                  br(),
                  sliderInput("r", # Name of variable in input object 
                              "Growth rate (r)", # Label shown on app
                              min = 0, # Min value on slider
                              max = 1, # Max value on slider
                              step = 0.1, # Step size on slider
                              value = 0.1 # Initial set value on slider
                              ),
                  br(),
                  sliderInput("K", "Carrying capacity (K)", min = 1, max = 10,value = 6),
                  br(),
                  sliderInput("pfr", "Proportional fishing rate (pfr)", min = 0, max = 1,value = 0.03),
                  br(),
                  sliderInput("cfr", "Constant fishing rate (cfr)", min = 0, max = 1,value = 0.13)
         ),
         tabPanel("Perturbation", # Starts tab panel for perturbation definition
                  br(),
                  sliderInput("perturb", "Perturbation", min = -10, max = 10,value = 0,step = 0.1),
                  br(),
                  sliderInput("pTime", "Time when perturbation occurs", min = 0, max = 100,value = 50,step = 0.1)
         )
       )
     ),
     mainPanel( # Starts main panel
       plotOutput("plotSim", # Name of plot in output object
                  height = 350), # Defines height of plot
       withMathJax(), # Needed to be able to include Latex equations
       br(),
       h5("Proportional fishing model"),
       # Latex equations are entered as strings (in quotes) starting and ending with two dollar signs ($$)
       "$$dNpf = r \\cdot Npf \\frac{K-Npf}{K} - pfr \\cdot Npf $$", # Latex equation for proporotional fishing
       br(),
       h5("Constant fishing model"),
       "$$dNcf = r \\cdot Ncf \\frac{K-Ncf}{K} - cfr $$" # Latex equation for constant fishing
     )
   )
)
# The server is written as a function with three inputs
# input is the object where all inputs from the interface (sliders etc) will be stored
# output is the object where all plots, text outputs etc. will be stored/accessed
# session can be used to define other variables/options needed for the app 
# (session is not used in this example, but is needed in the function definition)
server <- function(input, output, session) {

  # Reactive code is defined using the reactive command, which generates a reactive function
  # Here we are simply using reactive to generate a reactive instance of our runModel function
  # outModel is a function that calls/runs our runModel function using the latest values in the input object
  outModel <- reactive({
    runModel(input)
  })
  
  # Plots are generated using renderPlot and are saved in the output object with unique names
  # This plot is named plotSim and this name is used in the ui to place the plot on the page
  output$plotSim <- renderPlot({
    out <- outModel() # Call reactive function to get updated simulation 
    
    # Plots results of simulation and is ran everytime the input is updated
    par(mar=c(4,4,1,1))
    plot(out[,1],out[,2],ylim=c(0,10),xlab="Time",ylab="Population",col=4,lwd=2,type="l",cex.lab=1.3,cex.axis=1.3)
    lines(out[,1],out[,3],xlab="Time",ylab="",col=2,lwd=2,type="l")
    legend("topright",c("Npf","Ncf"),lwd = 2, col=c(4,2))
    
  })

}
shinyApp(ui = ui, server = server)
