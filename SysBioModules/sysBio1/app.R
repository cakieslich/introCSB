library(deSolve)
library(shiny)

runBlueSky <- function(input){
  ## Chaos in the atmosphere
  blueSky <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dX <- Y
      dY <- a * X - b * X^3 - c * Y + d * sin(e*t)
      list(c(dX, dY))
    })
  }
  parameters <- c(a = input$parA, b = input$parB,c = input$parC, d = input$parD, e = input$parE)
  state <- c(X = input$x0, Y = input$y0)
  times <- seq(0, 500, by = 0.1)
  out <- ode(y = state, times = times, func = blueSky, parms = parameters)
  
  return(out)
}

runSIR <- function(input){
  ## Chaos in the atmosphere
  SIRQ <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- rBirth + rS*R - rI*S*I + rRI*Q
      dI <- rI*S*I - rR*I - rDeath*I - rQ*I
      dR <- rR*I - rS*R
      dQ <- rQ*I - rRI*Q
      list(c(dS, dI, dR, dQ))
    })
  }
  parameters <- c(rBirth = input$rBirth, rDeath = input$rDeath,rI = input$rI, rR = input$rR, rS = input$rS, rQ = input$rQ, rRI = input$rRI)
  state <- c(S = 1000 - input$I0, I = input$I0, R = 0, Q = 0)
  times <- seq(0,input$time, by = 1)
  out <- ode(y = state, times = times, func = SIRQ, parms = parameters)
  
  parameters <- c(rBirth = input$rBirth, rDeath = input$rDeath,rI = input$rI, rR = input$rR, rS = input$rS, rQ = 0, rRI = input$rRI)
  out2 <- ode(y = state, times = times, func = SIRQ, parms = parameters)
  
  return(list(out,out2))
}

runSignal <- function(input){
  ## Chaos in the atmosphere
  signalCas <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      if(t>timeOn&t<timeOff){
        SS <- S
      } else{
        SS <- 0
      }
      dP0 <- rdP * P1 - rP * SS * P0
      dP1 <- rdP * P2 - rP * P1 - rdP * P1 + rP * SS * P0
      dP2 <- rP * P1 - rdP * P2
      list(c(dP0, dP1, dP2))
    })
  }
  parameters <- c(rP = input$rP, rdP = input$rdP,S = input$S, timeOn = input$timeOn, timeOff = input$timeOff)
  state <- c(P0 = 10, P1 = 0, P2 = 0)
  times <- seq(0, 100, by = 0.1)
  out <- ode(y = state, times = times, func = signalCas, parms = parameters)
  
  return(out)
}

dat <<- read.table("expDat.txt",header=T)

runParEst <- function(input){
  ## Chaos in the atmosphere
  hillEq <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dX <- a*Y^g - Vmax*X^4/(KM^4 + X^4)
      dY <- Vmax*X^4/(KM^4 + X^4) - b*Y^h
      list(c(dX, dY))
    })
  }
  parameters <- c(Vmax = input$Vmax, KM = input$KM,a = input$a, b = input$b, g = input$g, h = input$h)
  state <- c(X = 0.1, Y = 0.8)
  times <- seq(0, 10, by = 0.1)
  out <- ode(y = state, times = times, func = hillEq, parms = parameters)
  
  return(out)
}

ui <- fluidPage(
    tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 1000px; /* or 950px */}
                                          body { color: #222; }")),
    br(),
    hr(),
    titlePanel(NULL,windowTitle = "BSM1: Intro to Modeling"),
    h3("Biological Systems Modeling - Module 1: Intro to Modeling"),
    br(),
    tabsetPanel(
      tabPanel("Overview",
               br(),
               p("This is the first in a series of modules that are intended to provide simple interfaces
                  to example systems models and to provide extra practice problems. The modules are developed
                  using R Shiny applications, and include \"widgets\" that allow easy interaction with 
                  various models. These modules are only intended as a supplement and are not a sufficient 
                  replacement of a complete reading of the textbook."),
               p(code("Module 1: Intro to Modeling"),"provides interfaces to selected examples from",
                  strong("Chapter 2: Introduction to Mathematical Modeling"),"of",em("A First Course in Systems 
                  Biology."),"The expected outcomes of this module are  (i) learn the parts of a model and (ii) 
                  practice predicting how perturbations to model inputs and parameters affect model response"),
               tags$hr(),
               h4("Examples"),
               tags$ol(
                 tags$li("Blue Sky Catastrophe (p. 31)"),
                 tags$ul(
                   tags$li("Interesting system that exhibits chaos, but is completely deterministic.")
                 ),
                 tags$li("Protein Phosphorylation (p. 32)"),
                 tags$ul(
                   tags$li("A simple signaling model that contains a single protein that can exist in three 
                           distinct states: unphosphorylated (P0), singly phosphorylated (P1), and doubly
                           phosphorylated. Phosphorylation of the protein is driven by the presence of a 
                           signal (S)."
                   )
                 ),
                 tags$li("SIR Model (p. 35)"),
                 tags$ul(
                   tags$li("A simple epidemiological compartment model that predicts the spread of an infectious 
                           disease. In addition to an infected population, the model also includes a recovered, 
                           or at least temporarily immune, population, as well as, a susceptible population of 
                           people who can catch the disease from an infected individual.")
                 )
               ),
               h4("Futher Practice"),
               tags$ol(
                 tags$li("Parameter Estimation Demo"),
                 tags$ul(
                   tags$li("Parameter estimation by trial and error. The exercise has two purposes. First, develop a 
                       feel for the effect of parameter values on the dynamics of a system. Second, describe in 
                       detail (step by step) how you select new parameter values."
                   )
                 )
               ),
               tags$hr(),
               h4("Citation"),
               p("Voit, E.O.:",em("A First Course in Systems Biology."),"Garland Science, New York, NY, 2017, 2nd edition"),
               p("**All page numbers reference the Second Edition**"),
               br()
      ),              
      tabPanel("Blue Sky Catastrophe",
        br(),
        h4("Blue Sky Catastrophe"),
        tags$ul(
           tags$li("Interesting system that exhibits chaos, but is completely deterministic."),
           tags$li("Change the initial conditions and parameters to observe their effects on the system.")
        ),
        br(),
        sidebarLayout(
        sidebarPanel(
          tabsetPanel(
            tabPanel("Initial Conditions",
              br(),
              sliderInput("x0", "X0", min = 0, max = 1,value = 0.9),
              br(),
              sliderInput("y0", "Y0", min = 0, max = 1,value = 0.4)
            ),
            tabPanel("Parameters",
              br(),
              sliderInput("parA", "A", min = 0, max = 1,value = 1),
              br(),
              sliderInput("parB", "B", min = 0, max = 1,value = 1),
              br(),
              sliderInput("parC", "C", min = 0, max = 1,value = 0.25),
              br(),
              sliderInput("parD", "D", min = 0, max = 1,value = 0.2645),
              br(),
              sliderInput("parE", "E", min = 0, max = 1,value = 1)
            )
          )
        ),
        mainPanel(
          plotOutput("plotBlueSky",height = 350),
          withMathJax(),
          br(),
          "$$dX = Y$$",
          "$$dY = A \\cdot X - B \\cdot X^3 - C \\cdot Y + D \\cdot sin(E \\cdot t)$$"
        )
      )
    ),
    tabPanel("Protein Phosphorylation",
             br(),
             h4("Protein Phosphorylation"),
             tags$ul(
               tags$li("A simple signaling model that contains a single protein that can exist in three 
                       distinct states: unphosphorylated (P0), singly phosphorylated (P1), and doubly
                       phosphorylated. Phosphorylation of the protein is driven by the presence of a 
                       signal (S)."
               ),
               tags$li("Explore how the amount and duration of the signal effect the phosphorylation state."
               ),
               tags$li("What needs to changed in the model so that all of protein is in the doubly phosphorylated
                       state?"
               )
             ),
             br(),
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(
                   tabPanel("Signal",
                            br(),
                            sliderInput("S", "Signal (S)", min = 0, max = 5,value = 0,step = 0.1),
                            br(),
                            sliderInput("timeOn", "Time when Signal turns ON", min = 0, max = 100,value = 10,step = 0.1),
                            br(),
                            sliderInput("timeOff", "Time when Signal turns OFF", min = 0, max = 100,value = 30,step = 0.1)
                   ),
                   tabPanel("Parameters",
                            br(),
                            sliderInput("rP", "Rate of phosphorylation", min = 0, max = 1,value = 0.1),
                            br(),
                            sliderInput("rdP", "Rate of dephosphorylation", min = 0, max = 1,value = 0.1)
                            
                   )
                 )
               ),
               mainPanel(
                 plotOutput("plotSignal",height = 350),
                 withMathJax(),
                 br(),
                 "$$dP0 = rdP \\cdot P1 - rP \\cdot S \\cdot P0$$",
                 "$$dP1 = rdP \\cdot P2 - rP \\cdot P1 - rdP \\cdot P1 + rP \\cdot S \\cdot P0$$",
                 "$$dP2 = rP \\cdot P1 - rdP \\cdot P2$$"
               )
             )
    ),
    tabPanel("SIR Models",
             br(),
             h4("SIR Model"),
            tags$ul(
              tags$li("A simple epidemiological compartment model that predicts the spread of an infectious 
                            disease. In addition to an infected population, the model also includes a recovered, 
                            or at least temporarily immune, population, as well as, a susceptible population of 
                            people who can catch the disease from an infected individual."),
              tags$li("Explore the effects of the various rates on the size of the infected population."),
              tags$li("Can you eradicate the disease using quarantine?"),
              tags$li("What is the effect of the initial number of infected individuals on the effectiveness 
                           of quarantine?")
              ),
            br(),
            sidebarLayout(
            sidebarPanel(
              tabsetPanel(
                tabPanel("Time/Init. Cond.",
                  br(),
                  sliderInput("time", "Final Time", min = 0, max = 730,value = 100, step = 1),
                  br(),
                  sliderInput("I0", "Initial # of Infected Individuals (I0)", min = 0, max = 1000,value = 10, step = 1),
                  br()
                ),
                tabPanel("Parameters",
                  br(),
                  sliderInput("rBirth", "Birth Rate (rBirth)", min = 0, max = 100,value = 3, step = 1),
                  br(),
                  sliderInput("rDeath", "Death Rate (rDeath)", min = 0, max = 1,value = 0.02, step = 0.01),
                  br(),
                  sliderInput("rI", "Infection Rate (rI)", min = 0, max = 0.001,value = 0.0005, step = 0.00001),
                  br(),
                  sliderInput("rR", "Immunity Rate (rR)", min = 0, max = 1,value = 0.05, step = 0.01),
                  br(),
                  sliderInput("rS", "Susceptibility Rate (rS)", min = 0, max = 1,value = 0.01, step = 0.01),
                  br(),
                  sliderInput("rQ", "Quarantine Rate (rQ)", min = 0, max = 1,value = 0, step = 0.01),
                  br(),
                  sliderInput("rRI", "Reintroduction Rate (rRI)", min = 0, max = 1,value = 0, step = 0.01)
                )
               )
               ),
               mainPanel(
                 plotOutput("plotSIR",height = 350),
                 br(),
                 h5("Solid lines - With quarantine; Dashed lines - Without quarantine",align="center"),
                 "$$dS = rBirth + rS \\cdot  R - rI \\cdot  S \\cdot  I + rRI \\cdot  Q$$",
                 "$$dI = rI \\cdot  S \\cdot  I - rR \\cdot  I - rDeath \\cdot  I - rQ \\cdot  I$$",
                 "$$dR = rR \\cdot  I - rS \\cdot  R$$",
                 "$$dQ = rQ \\cdot  I - rRI \\cdot  Q$$"
               )
             )
    ),
    tabPanel("Parameter Estimation",
             br(),
             h4("Parameter Estimation"),
             tags$ul(
               tags$li("Parameter estimation by trial and error. The exercise has two purposes. First, develop a 
                       feel for the effect of parameter values on the dynamics of a system. Second, describe in 
                       detail (step by step) how you select new parameter values."
               ),
               tags$li("Enter values for each of the 6 parameters, press \"Try Parameters\", and the app will 
                       run a simulation to calculate the error between the model (lines) and the provided data (points)."
               ),
               tags$li("Identify parameters that minimize the error in as few trials as possible." 
               )
             ),
             br(),
             sidebarLayout(
               sidebarPanel(
                 h4("Parameters"),
                 numericInput("Vmax", "Vmax",value=0),
                 numericInput("KM", "KM",value=0),
                 numericInput("a", "a",value=0),
                 numericInput("b", "b",value=0),
                 numericInput("g", "g",value=0),
                 numericInput("h", "h",value=0),
                 br(),
                 actionButton("run", "Try Parameters")
               ),
               mainPanel(
                 plotOutput("plotParEst"),
                 withMathJax(),
                 "$$dX = a Y^g - \\frac{V_{max} \\cdot X^4}{K_M^4 + X^4}$$",
                 "$$dY = \\frac{V_{max} \\cdot X^4}{K_M^4 + X^4} - b Y^h$$"
               )
             )
    )
  )
)
server <- function(input, output, session) {
  
  outBlueSky <- reactive({
    runBlueSky(input)
  })
  
  output$plotBlueSky <- renderPlot({
    out <-  outBlueSky()
    
    par(mar=c(4,4,1,1))
    plot(out[,1],out[,2],xlab="Time",ylab="X",col=4,lwd=2,type="l",cex.lab=1.3,cex.axis=1.3)
    
  })
  
  outSignal <- reactive({
    runSignal(input)
  })
  
  output$plotSignal <- renderPlot({
    out <- outSignal()
    
    par(mar=c(4,4,1,1))
    plot(out[,1],out[,2],ylim=c(0,10),xlab="Time",ylab="Concentration",col=4,lwd=2,type="l",cex.lab=1.3,cex.axis=1.3)
    lines(out[,1],out[,3],xlab="Time",ylab="",col=2,lwd=2,type="l")
    lines(out[,1],out[,4],xlab="Time",ylab="",col=3,lwd=2,type="l")
    legend("topright",c("P0","P1","P2"),lwd = 2, col=c(4,2,3))
    
  })
  
  outSIR <- reactive({
    runSIR(input)
  })
  
  output$plotSIR <- renderPlot({
    out <- outSIR()
      
    par(mar=c(4,4,1,1))
    plot(out[[1]][,1],out[[1]][,2],ylim=c(0,1000),xlab="Time",ylab="Number of individuals",col=4,lwd=2,type="l",cex.lab=1.3,cex.axis=1.3)
    lines(out[[1]][,1],out[[1]][,3],xlab="Time",ylab="",col=2,lwd=2,type="l")
    lines(out[[1]][,1],out[[1]][,4],xlab="Time",ylab="",col=3,lwd=2,type="l")
    lines(out[[1]][,1],out[[1]][,5],xlab="Time",ylab="",col=1,lwd=2,type="l")
    
    lines(out[[2]][,1],out[[2]][,2],col=4,lwd=2,type="l",lty=2)
    lines(out[[2]][,1],out[[2]][,3],xlab="Time",ylab="",col=2,lwd=2,type="l",lty=2)
    lines(out[[2]][,1],out[[2]][,4],xlab="Time",ylab="",col=3,lwd=2,type="l",lty=2)
    legend("topright",c("S","I","R","Q"),lwd = 2, col=c(4,2,3,1))
    
  })
  
  outParEst <- eventReactive(input$run,{
    runParEst(input)
  })
  
  output$plotParEst <- renderPlot({
    out <- outParEst()
    
    par(mar=c(4,4,1,1))
    plot(dat$t,dat$X,ylim=c(0,2),pch=19,col=2,cex=2,xlab="Time",ylab="",cex.lab=1.3,cex.axis=1.3)
    points(dat$t,dat$Y,pch=19,col=4,cex=2)
    lines(out[,1],out[,2],col=2,lwd=2)
    lines(out[,1],out[,3],col=4,lwd=2)
    legend("topright",c("X","Y"),lwd = 2, col=c(2,4))
    error <- sum(abs(out[out[,1]%in%dat$t,2]-dat$X))+sum(abs(out[out[,1]%in%dat$t,3]-dat$Y))
    text(0,1.8,paste("Number of trys:",input$run),cex=1.5,pos=4)
    text(0,1.6,sprintf("Absolute Error: %.3f",error),cex=1.5,pos=4)
  })
  
  
}
shinyApp(ui = ui, server = server)