library(deSolve)
library(shiny)
library(DiagrammeR)

MAPK <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    if(t > startS & t < endS){
      S <- A*sin(2*pi*freq*t) + A/2  
    } else{
      S <- 0
    }
    
    dX <- -kon*X*E*S + koff*XE + kdp*XP
    dE <- -kon*X*E*S + koff*XE + kcat*XE
    dXE <- kon*X*E*S - koff*XE - kcat*XE
    dXP <- kcat*XE - kon*Y*XP + koff*YXP -kdp*XP + kcat*YXP
    dY <- -kon*Y*XP + koff*YXP + kdp*YP
    dYXP <- kon*Y*XP - koff*YXP - kcat*YXP
    dYP <- kcat*YXP - kon*Z*YP + koff*ZYP - kdp*YP + kcat*ZYP
    dZ <- -kon*Z*YP + koff*ZYP + kdp*ZP
    dZYP <- kon*Z*YP - koff*ZYP - kcat*ZYP
    dZP <- kcat*ZYP - kdp*ZP
    
    list(c(dX, dE, dXE, dXP, dY, dYXP, dYP, dZ, dZYP, dZP))
  })
}

runModel <- function(input){
  ## Chaos in the atmosphere
  
  # All paramters in this example are controlled by sliders in the interface
  # As the slider is moved the reactive variable input is updated with new values
  parameters <- c(kon = input$kon, 
                  koff = input$koff,
                  kcat = input$kcat, 
                  kdp = input$kdp,
                  A = input$A,
                  freq = input$freq,
                  startS = input$startS,
                  endS = input$endS
                  )
  
  state <- c(X = 50, 
             E = 5, 
             XE = 0, 
             XP = 0, 
             Y = 50, 
             YXP = 0, 
             YP = 0, 
             Z = 50, 
             ZYP = 0, 
             ZP = 0)
  
  times <- seq(0, 100, by = 0.1)
  
  out <- ode(y = state, 
             times = times, 
             func = MAPK, 
             parms = parameters
             )
  
  return(out)
}

# A shiny app has two main components: a user interface (ui) and a server
# The ui defines the layout of the application
# The server conects to interface to the underlying model/code

# ui is ultimately generated in HTML so some familiarity with HTML can help but is not necessary
# Many of the functions used to build ui are named after an HTML equivalent. 
ui <- fluidPage(
    h3("Mitogen-Activated Protein Kinase Cascades"), # Header level 3 (1 is largest and font decreases with increasing number)
    br(), # Line break, leaves an empty line
    
   h4("Simplified MAPK cascade model"), # Header level 4
   tags$ul( # Beginning of a bulleted list
     # li are each items in the bulleted list
     tags$li("In this simplified MAPK model, we have assumed that each layer only consists of one 
             phosphorylation reaction (in reality the second and third layers require two phosphorylations)."
     ),
     tags$li("This simplified model also assumes the rates for phosphorylation/dephosphorylation are the same 
             for all three levels of the MAPK cascade."
     ),
     tags$li("MAPK cascades aid in intrcellular signaling through amplification and noise reduction of incoming signals."
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
                  sliderInput("kon", # Name of variable in input object 
                              "Association rate (kon)", # Label shown on app
                              min = 0, # Min value on slider
                              max = 1, # Max value on slider
                              step = 0.1, # Step size on slider
                              value = 0 # Initial set value on slider
                              ),
                  br(),
                  sliderInput("koff", "Dissociation rate (koff)", min = 0, max = 1,value = 0,step = 0.1),
                  br(),
                  sliderInput("kcat", "Catalytic rate (kcat)", min = 0, max = 3,value = 0,step = 0.1),
                  br(),
                  sliderInput("kdp", "Dephosphorylation rate (kdp)", min = 0, max = 10,value = 0,step = 0.1)
         ),
         tabPanel("Signal", # Starts tab panel for perturbation definition
                  br(),
                  sliderInput("A", "Signal amplitude (A)", min = 0, max = 1,value = 0.1,step = 0.1),
                  br(),
                  sliderInput("freq", "Sginal frequence (freq)", min = 0, max = 1,value = 0.5,step = 0.1),
                  br(),
                  sliderInput("startS", "Signal start time (startS)", min = 0, max = 100, value = 25,step = 0.1),
                  br(),
                  sliderInput("endS", "Signal end time (endS)", min = 0, max = 100, value = 75,step = 0.1)
         )
       )
     ),
     mainPanel( # Starts main panel
       tabsetPanel( # Begins a tabbed layoutt
         tabPanel("Plot",
                  plotOutput("plotSim", # Name of plot in output object
                              height = 350), # Defines height of plot
                  br(),
                  column(4,downloadButton("report", "Download Report"))
         ),
         tabPanel("Diagram",
                  grVizOutput("diagram", width = "100%", height = "400px") 
         )
       )
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
    plot(out[,1],out[,5],ylim=c(0,ceiling(max(out[,c(5,8,11)]))),xlab="Time",ylab="Population",col=4,lwd=2,type="l",cex.lab=1.3,cex.axis=1.3)
    lines(out[,1],out[,8],xlab="Time",ylab="",col=2,lwd=2,type="l")
    lines(out[,1],out[,11],xlab="Time",ylab="",col="forestgreen",lwd=2,type="l")
    legend("topright",c("XP","YP","ZP"),lwd = 2, col=c("blue","red","forestgreen"))
    
  })
  
  output$diagram <- renderGrViz({
    grViz("digraph {

                      # Initialize graph
                      graph [layout = dot, rankdir = LR]
                      
                      # Define the nodes and global styles of the nodes. We can override these in box if we wish
                      node [shape = circle, style = filled, fixedsize = true,fontcolor=white]
                      
                      #Edges in graphs must be between nodes, therefore dummy nodes are needed for signaling edges and input/output edges
                      
                      S [fillcolor = dimgray]

                      E [fillcolor = blue]
                      X [fillcolor = blue]
                      fEX [style=invis, shape=point, width=0] #Dummy node
                      rEX [style=invis, shape=point, width=0] #Dummy node
                      XE [fillcolor = blue]
                      iXP [style=invis, shape=point, width=0] #Dummy node
                      XP [fillcolor = blue]
                      
                      Y [fillcolor = red]
                      fYXP [style=invis, shape=point, width=0] #Dummy node
                      rYXP [style=invis, shape=point, width=0] #Dummy node
                      YXP [fillcolor = red]
                      iYP [style=invis, shape=point, width=0] #Dummy node
                      YP [fillcolor = red]

                      Z [fillcolor = ForestGreen]
                      fZYP [style=invis, shape=point, width=0] #Dummy node
                      rZYP [style=invis, shape=point, width=0] #Dummy node
                      ZYP [fillcolor = ForestGreen]
                      iZP [style=invis, shape=point, width=0] #Dummy node
                      ZP [fillcolor = ForestGreen]
          

                      # edge definitions with the node IDs
                      edge[fontsize = 10]
                      # Edge weights are not necessary, but can be used to get edges in straight line (higher weight edges have more emphasis)
                      
                      ranksep=0.2;
                      
                      S -> fEX [label = '+', tailport = n, style = dashed, fontcolor = green, color = green,weight=0]

                      E -> fEX [dir = none, weight=1]
                      X -> fEX [dir = none, weight=2]
                      fEX -> XE [label = 'k@_{on}', weight=2]
                      rEX -> X [weight = 2]
                      rEX -> E [weight = 1]
                      XE -> rEX [weight = 2, dir = none, label = 'k@_{off}']
                      XE -> iXP [label = 'k@_{cat}', weight = 2, dir = none]
                      iXP -> XP [weight = 2]
                      iXP -> E [weight = 1]
                      XP -> X [label = 'k@_{dp}',weight = 1]
                      
                      XP -> fYXP [dir = none, weight=1]
                      Y -> fYXP [dir = none, weight=2]
                      fYXP -> YXP [label = 'k@_{on}', weight=2]
                      rYXP -> Y [weight = 2]
                      rYXP -> XP [weight = 1]
                      YXP -> rYXP [weight = 2, dir = none, label = 'k@_{off}']
                      YXP -> iYP [label = 'k@_{cat}', weight = 2, dir = none]
                      iYP -> YP [weight = 2]
                      iYP -> XP [weight = 1]
                      YP -> Y [label = 'k@_{dp}',weight = 1]

                      YP -> fZYP [dir = none, weight=1]
                      Z -> fZYP [dir = none, weight=2]
                      fZYP -> ZYP [label = 'k@_{on}', weight=2]
                      rZYP -> Z [weight = 2]
                      rZYP -> YP [weight = 1]
                      ZYP -> rZYP [weight = 2, dir = none, label = 'k@_{off}']
                      ZYP -> iZP [label = 'k@_{cat}', weight = 2, dir = none]
                      iZP -> ZP [weight = 2]
                      iZP -> YP [weight = 1]
                      ZP -> Z [label = 'k@_{dp}',weight = 1]
                  }")
  })
  
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "mapk_report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "mapk_report.Rmd")
      file.copy("mapk_report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(input=input,output=outModel())
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )

}
shinyApp(ui = ui, server = server)