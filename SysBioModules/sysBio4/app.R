library(deSolve)
library(shiny)
library(plotly)

runLinear <- function(input){
  ## Chaos in the atmosphere
  linMod <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dX <- rxin + rxx*X + rxy*Y + rxz*Z
      dY <- ryin + ryx*X + ryy*Y + ryz*Z
      dZ <- rzin + rzx*X + rzy*Y + rzz*Z
    
      list(c(dX, dY, dZ))
    })
  }
  parameters <- c(rxin = input$rxin,rxx = input$rxx,rxy = input$rxy,rxz = input$rxz, 
                  ryin = input$ryin,ryx = input$ryx,ryy = input$ryy,ryz = input$ryz,
                  rzin = input$rzin,rzx = input$rzx,rzy = input$rzy,rzz = input$rzz
                  )
  
  state <- c(X = input$x0, Y = input$y0, Z = input$z0)
  times <- seq(0, 300, by = 0.1)
  out <- ode(y = state, times = times, func = linMod, parms = parameters)
  
  return(data.frame(out))
}

runMA <- function(input){
  ## Chaos in the atmosphere
  maMod <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dX <- rxin + rxx*X + rxy*Y + rxz*Z + rxxx*X*X + rxxy*X*Y + rxxz*X*Z + rxxyz*X*Y*Z
      dY <- ryin + ryx*X + ryy*Y + ryz*Z + ryxx*X*X + ryxy*X*Y + ryxz*X*Z + ryxyz*X*Y*Z
      dZ <- rzin + rzx*X + rzy*Y + rzz*Z + rzxx*X*X + rzxy*X*Y + rzxz*X*Z + rzxyz*X*Y*Z
      
      list(c(dX, dY, dZ))
    })
  }
  parameters <- c(rxin = input$mrxin,rxx = input$mrxx,rxy = input$mrxy,rxz = input$mrxz,rxxx = input$mrxxx,rxxy = input$mrxxy,rxxz = input$mrxxz,rxxyz = input$mrxxyz, 
                  ryin = input$mryin,ryx = input$mryx,ryy = input$mryy,ryz = input$mryz,ryxx = input$mryxx,ryxy = input$mryxy,ryxz = input$mryxz,ryxyz = input$mryxyz,
                  rzin = input$mrzin,rzx = input$mrzx,rzy = input$mrzy,rzz = input$mrzz,rzxx = input$mrzxx,rzxy = input$mrzxy,rzxz = input$mrzxz,rzxyz = input$mrzxyz
  )
  
  state <- c(X = input$mx0, Y = input$my0, Z = input$mz0)
  times <- seq(0, 300, by = 0.1)
  out <- ode(y = state, times = times, func = maMod, parms = parameters)
  
  return(data.frame(out))
}

runLV <- function(input){
  ## Chaos in the atmosphere
  lvMod <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dX <- rX*X + rXX*X*X + rXY*X*Y + rXZ*X*Z
      dY <- rY*Y + rYY*Y*Y + rYX*Y*X + rYZ*Y*Z
      dZ <- rZ*Z + rZZ*Z*Z + rZX*Z*X + rZY*Z*Y
      
      list(c(dX, dY, dZ))
    })
  }
  parameters <- c(rX = input$rX, rXX = input$rXX,rXY = input$rXY,rXZ = input$rXZ, 
                  rY = input$rY, rYY = input$rYY,rYX = input$rYX,rYZ = input$rYZ,
                  rZ = input$rZ, rZZ = input$rZZ,rZX = input$rZX,rZY = input$rZY)
  
  state <- c(X = input$X0, Y = input$Y0, Z = input$Z0)
  times <- seq(0, 50, by = 0.01)
  out <- ode(y = state, times = times, func = lvMod, parms = parameters)
  
  return(data.frame(out))
}

runSSys <- function(input){
  ## Chaos in the atmosphere
  ssysMod <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dX <- a1*X^g11*Y^g12*Z^g13 - b1*X^h11*Y^h12*Z^h13
      dY <- a2*X^g21*Y^g22*Z^g23 - b2*X^h21*Y^h32*Z^h23
      dZ <- a3*X^g31*Y^g32*Z^g33 - b3*X^h31*Y^h32*Z^h33
      
      list(c(dX, dY, dZ))
    })
  }
  parameters <- c(a1 = input$a1, b1 = input$b1, g11 = input$g11,g12 = input$g12,g13 = input$g13, h11 = input$h11,h12 = input$h12,h13 = input$h13, 
                  a2 = input$a2, b2 = input$b2, g21 = input$g21,g22 = input$g22,g23 = input$g23, h21 = input$h21,h22 = input$h22,h23 = input$h23,
                  a3 = input$a3, b3 = input$b3, g31 = input$g31,g32 = input$g32,g33 = input$g33, h31 = input$h31,h32 = input$h32,h33 = input$h33)
  
  state <- c(X = input$sX0, Y = input$sY0, Z = input$sZ0)
  times <- seq(0, 20, by = 0.01)
  out <- ode(y = state, times = times, func = ssysMod, parms = parameters)
  
  return(data.frame(out))
}

ui <- fluidPage(
  tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 1000px; /* or 950px */}")),
  titlePanel("Biomedical Systems and Modeling - Module 4: Canonical Models"),
  tabsetPanel(
    tabPanel("Overview",
             br(),
             p("This is the fourth of the Biomedical Systems and Modeling modules that are intended to provide simple interfaces
                  to example systems models and to provide extra practice problems, and is the follow up to ",code("Module 3: Discrete, Recursive Models")),
             p(code("Module 4: Canonical Models"),"provides interfaces to example models based on the canonical models described in pg. 102-105 of ",
               strong("Chapter 4: The Mathematics of Biological Systems"),"of",em("A First Course in Systems 
                  Biology."),"The expected outcomes of this module are to learn the basic structure of (i) linear, (ii) Lotka-Voltera, (iii) mass action, and (iv) S-system models."),
             tags$hr(),
             h4("Examples"),
             tags$ol(
               tags$li("Linear Model"),
               tags$ul(
                 tags$li("Interface to a dynamic model with 3 dependent variables composed of linear ordinary 
                       differential equations.")
               ),
               tags$li("Mass Action Model"),
               tags$ul(
                 tags$li("Interface to a dynamic model with 3 dependent variables composed of mass action equations.")
                ),
               tags$li("Lotka-Volterra Model"),
               tags$ul(
                 tags$li("Interface to a dynamic model with 3 dependent variables composed of Lotka-Volterra equations."
                 )
               ),
               tags$li("S-system Model"),
               tags$ul(
                 tags$li("Interface to a dynamic model with 3 dependent variables composed of S-system equations.")
                 )
             ),
             tags$hr(),
             h4("Citation"),
             p("Voit, E.O.:",em("A First Course in Systems Biology."),"Garland Science, New York, NY, 2017, 2nd edition"),
             p("**All page numbers reference the Second Edition**"),
             br()
               ),              
    tabPanel("Linear Model",
             h3("Linear Model"),
             tags$ul(
               tags$li("Interface to a dynamic model with 3 dependent variables composed of linear ordinary 
                       differential equations. Simulation results are presented as time series, as well as,
                        in a phase plane plot, where the line color changes from purple to yellow over time."),
               tags$li("Change the initial conditions and parameters to observe their effects on the system.")
             ),
             br(),
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(
                   tabPanel("I.C.",
                            br(),
                            sliderInput("x0", "X0", min = -10, max = 10,value = -5,step=1),
                            br(),
                            sliderInput("y0", "Y0", min = -10, max = 10,value = -2,step=1),
                            br(),
                            sliderInput("z0", "Z0", min = -10, max = 10,value = 2,step=1)
                   ),
                   tabPanel("X",
                            br(),
                            sliderInput("rxin", "rxin", min = -2, max = 2,value = 2,step=0.1),
                            br(),
                            sliderInput("rxx", "rxx", min = -2, max = 2,value = -0.2,step=0.1),
                            br(),
                            sliderInput("rxy", "rxy", min = -2, max = 2,value = -0.5,step=0.1),
                            br(),
                            sliderInput("rxz", "rxz", min = -2, max = 2,value = -0.3,step=0.1)
                   ),
                   tabPanel("Y",
                            br(),
                            sliderInput("ryin", "ryin", min = -2, max = 2,value = 1.2,step=0.1),
                            br(),
                            sliderInput("ryx", "ryx", min = -2, max = 2,value = 0.3,step=0.1),
                            br(),
                            sliderInput("ryy", "ryy", min = -2, max = 2,value = 0.1,step=0.1),
                            br(),
                            sliderInput("ryz", "ryz", min = -2, max = 2,value = 0,step=0.1)
                   ),
                   tabPanel("Z",
                            br(),
                            sliderInput("rzin", "rzin", min = -2, max = 2,value = -1.2,step=0.1),
                            br(),
                            sliderInput("rzx", "rzx", min = -2, max = 2,value = -0.2,step=0.1),
                            br(),
                            sliderInput("rzy", "rzy", min = -2, max = 2,value = 0,step=0.1),
                            br(),
                            sliderInput("rzz", "rzz", min = -2, max = 2,value = 0,step=0.1)
                   )
                 )
               ),
               mainPanel(
                 fluidRow(
                   column(6,plotlyOutput("plotLinear")),
                   column(6,plotlyOutput("plotLinPP"))
                 ),
                 withMathJax(),
                 br(),
                 "$$dX = r_{X,in} + r_{XX} \\cdot X + r_{XY} \\cdot Y + r_{XZ} \\cdot  Z$$",
                 "$$dY = r_{Y,in} + r_{YX} \\cdot X + r_{YY} \\cdot Y + r_{YZ} \\cdot  Z$$",
                 "$$dZ = r_{Z,in} + r_{ZX} \\cdot X + r_{ZY} \\cdot Y + r_{ZZ} \\cdot  Z$$"
               )
             )
    ),
    tabPanel("Mass Action Model",
             h3("Mass Action Model"),
             tags$ul(
               tags$li("Interface to a dynamic model with 3 dependent variables composed of mass action 
                        equations. Simulation results are presented as time series, as well as,
                       in a phase plane plot, where the line color changes from purple to yellow over time."),
               tags$li("Change the initial conditions and parameters to observe their effects on the system.")
               ),
             br(),
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(
                   tabPanel("I.C.",
                            br(),
                            sliderInput("mx0", "X0", min = -10, max = 10,value = 0,step=1),
                            br(),
                            sliderInput("my0", "Y0", min = -10, max = 10,value = 0,step=1),
                            br(),
                            sliderInput("mz0", "Z0", min = -10, max = 10,value = 0,step=1)
                   ),
                   tabPanel("X",
                            br(),
                            sliderInput("mrxin", "rxin", min = -2, max = 2,value = 0.4,step=0.1),
                            br(),
                            sliderInput("mrxx", "rxx", min = -2, max = 2,value = -0.1,step=0.1),
                            br(),
                            sliderInput("mrxy", "rxy", min = -2, max = 2,value = -0.3,step=0.1),
                            br(),
                            sliderInput("mrxz", "rxz", min = -2, max = 2,value = 0,step=0.1),
                            br(),
                            sliderInput("mrxxx", "rxxx", min = -0.2, max = 0.2,value = 0.01,step=0.01),
                            br(),
                            sliderInput("mrxxy", "rxxy", min = -0.2, max = 0.2,value = -0.06,step=0.01),
                            br(),
                            sliderInput("mrxxz", "rxxz", min = -0.2, max = 0.2,value = -0.04,step=0.01),
                            br(),
                            sliderInput("mrxxyz", "rxxyz", min = -0.02, max = 0.02,value = 0.004,step=0.001)
                   ),
                   tabPanel("Y",
                            br(),
                            sliderInput("mryin", "ryin", min = -2, max = 2,value = 0.7,step=0.1),
                            br(),
                            sliderInput("mryx", "ryx", min = -2, max = 2,value = -0.1,step=0.1),
                            br(),
                            sliderInput("mryy", "ryy", min = -2, max = 2,value = -1.1,step=0.1),
                            br(),
                            sliderInput("mryz", "ryz", min = -2, max = 2,value = 0,step=0.1),
                            br(),
                            sliderInput("mryxx", "ryxx", min = -0.2, max = 0.2,value = 0,step=0.01),
                            br(),
                            sliderInput("mryxy", "ryxy", min = -0.2, max = 0.2,value = 0,step=0.01),
                            br(),
                            sliderInput("mryxz", "ryxz", min = -0.2, max = 0.2,value = -0.04,step=0.01),
                            br(),
                            sliderInput("mryxyz", "ryxyz", min = -0.02, max = 0.02,value = -0.02,step=0.001)
                   ),
                   tabPanel("Z",
                            br(),
                            sliderInput("mrzin", "rzin", min = -2, max = 2,value = -0.2,step=0.1),
                            br(),
                            sliderInput("mrzx", "rzx", min = -2, max = 2,value = 0,step=0.1),
                            br(),
                            sliderInput("mrzy", "rzy", min = -2, max = 2,value = 0.2,step=0.1),
                            br(),
                            sliderInput("mrzz", "rzz", min = -2, max = 2,value = -0.5,step=0.1),
                            br(),
                            sliderInput("mrzxx", "rzxx", min = -0.2, max = 0.2,value = -0.01,step=0.01),
                            br(),
                            sliderInput("mrzxy", "rzxy", min = -0.2, max = 0.2,value = -0.07,step=0.01),
                            br(),
                            sliderInput("mrzxz", "rzxz", min = -0.2, max = 0.2,value = -0.02,step=0.01),
                            br(),
                            sliderInput("mrzxyz", "rzxyz", min = -0.02, max = 0.02,value = 0.005,step=0.001)
                   )
                 )
               ),
               mainPanel(
                 fluidRow(
                   column(6,plotlyOutput("plotMA")),
                   column(6,plotlyOutput("plotMAPP"))
                 ),
                 withMathJax(),
                 br(),
                 "$$dX = r_{X,in} + r_{XX} \\cdot X + r_{XY} \\cdot Y + r_{XZ} \\cdot  Z + r_{XXX} \\cdot X \\cdot X + r_{XXY} \\cdot X \\cdot Y + r_{XXZ} \\cdot X \\cdot  Z + r_{XXYZ} \\cdot X \\cdot  Y \\cdot  Z$$",
                 "$$dY = r_{Y,in} + r_{YX} \\cdot X + r_{YY} \\cdot Y + r_{YZ} \\cdot  Z + r_{YXX} \\cdot Y \\cdot X + r_{YXY} \\cdot X \\cdot Y + r_{YXZ} \\cdot X \\cdot  Z + r_{YXYZ} \\cdot X \\cdot  Y \\cdot  Z$$",
                 "$$dZ = r_{Z,in} + r_{ZX} \\cdot X + r_{ZY} \\cdot Y + r_{ZZ} \\cdot  Z + r_{ZXX} \\cdot Z \\cdot X + r_{ZXY} \\cdot X \\cdot Y + r_{ZXZ} \\cdot X \\cdot  Z + r_{ZXYZ} \\cdot X \\cdot  Y \\cdot  Z$$"
               )
             )
             ),
    tabPanel("LV Model",
             h3("Lotka-Volterra Model"),
             tags$ul(
               tags$li("Interface to a dynamic model with 3 dependent variables composed using Lotka-Volterra equations. 
                        Simulation results are presented as time series, as well as, in a phase plane plot, where the 
                        line color changes from purple to yellow over time."),
               tags$li("Change the initial conditions and parameters to observe their effects on the system.")
               ),
             br(),
             sidebarLayout(
               sidebarPanel(
                 tabsetPanel(
                   tabPanel("I.C.",
                            br(),
                            sliderInput("X0", "X0", min = 0, max = 4,value = 1.2,step=0.1),
                            br(),
                            sliderInput("Y0", "Y0", min = 0, max = 4,value = 0.6,step=0.1),
                            br(),
                            sliderInput("Z0", "Z0", min = 0, max = 4,value = 2.5,step=0.1)
                   ),
                   tabPanel("X",
                            br(),
                            sliderInput("rX", "rX", min = -3, max = 3,value = 0.84,step=0.01),
                            br(),
                            sliderInput("rXX", "rXX", min = -3, max = 3,value = 0.02,step=0.01),
                            br(),
                            sliderInput("rXY", "rXY", min = -3, max = 3,value = 0.19,step=0.01),
                            br(),
                            sliderInput("rXZ", "rXZ", min = -3, max = 3,value = -0.34,step=0.01)
                   ),
                   tabPanel("Y",
                            br(),
                            sliderInput("rY", "rY", min = -3, max = 3,value = -0.82,step=0.01),
                            br(),
                            sliderInput("rYY", "rYY", min = -3, max = 3,value = -0.04,step=0.01),
                            br(),
                            sliderInput("rYX", "rYX", min = -3, max = 3,value = 1.06,step=0.01),
                            br(),
                            sliderInput("rYZ", "rYZ", min = -3, max = 3,value = -0.07,step=0.01)
                   ),
                   tabPanel("Z",
                            br(),
                            sliderInput("rZ", "rZ", min = -3, max = 3,value = -3,step=0.01),
                            br(),
                            sliderInput("rZZ", "rZZ", min = -3, max = 3,value = -0.08,step=0.01),
                            br(),
                            sliderInput("rZX", "rZX", min = -3, max = 3,value = 2.23,step=0.01),
                            br(),
                            sliderInput("rZY", "rZY", min = -3, max = 3,value = 0,step=0.01)
                   )
                 )
               ),
               mainPanel(
                 fluidRow(
                   column(6,plotlyOutput("plotLV")),
                   column(6,plotlyOutput("plotLVPP"))
                 ),
                 withMathJax(),
                 br(),
                 "$$dX = r_{X} \\cdot X + r_{XX} \\cdot X \\cdot X + r_{XY} \\cdot X \\cdot Y + r_{XZ} \\cdot X \\cdot Z$$",
                 "$$dY = r_{Y} \\cdot Y + r_{YX} \\cdot Y \\cdot X + r_{YY} \\cdot Y \\cdot Y + r_{YZ} \\cdot Y \\cdot Z$$",
                 "$$dZ = r_{Z} \\cdot Z + r_{ZX} \\cdot Z \\cdot X + r_{ZY} \\cdot Z \\cdot Y + r_{ZZ} \\cdot Z \\cdot Z$$"
               )
             )
        ),
    tabPanel("S-system Model",
            h3("S-system Model"),
            tags$ul(
              tags$li("Interface to a dynamic model with 3 dependent variables composed of S-system 
                      equations. Simulation results are presented as time series, as well as, 
                      in a phase plane plot, where the line color changes from purple to yellow over time."),
              tags$li("Change the initial conditions and parameters to observe their effects on the system.")
              ),
            br(),
            sidebarLayout(
              sidebarPanel(
                tabsetPanel(
                  tabPanel("I.C.",
                           br(),
                           sliderInput("sX0", "X0", min = 0, max = 4,value = 0.2,step=0.1),
                           br(),
                           sliderInput("sY0", "Y0", min = 0, max = 4,value = 1,step=0.1),
                           br(),
                           sliderInput("sZ0", "Z0", min = 0, max = 4,value = 1,step=0.1)
                  ),
                  tabPanel("X", 
                           br(),
                           sliderInput("a1", "a1", min = 0, max = 3,value = 0.7,step=0.1),
                           br(),
                           sliderInput("g11", "g11", min = -3, max = 3,value = -1.1,step=0.1),
                           br(),
                           sliderInput("g12", "g12", min = -3, max = 3,value = 2,step=0.1),
                           br(),
                           sliderInput("g13", "g13", min = -3, max = 3,value = -2.3,step=0.1),
                           br(),
                           sliderInput("b1", "b1", min = 0, max = 3,value = 0.2,step=0.1),
                           br(),
                           sliderInput("h11", "h11", min = -3, max = 3,value = 1,step=0.1),
                           br(),
                           sliderInput("h12", "h12", min = -3, max = 3,value = -2,step=0.1),
                           br(),
                           sliderInput("h13", "h13", min = -3, max = 3,value = 0,step=0.1)
                  ),
                  tabPanel("Y",
                           br(),
                           sliderInput("a2", "a2", min = 0, max = 3,value = 1.5,step=0.1),
                           br(),
                           sliderInput("g21", "g21", min = -3, max = 3,value = 1,step=0.1),
                           br(),
                           sliderInput("g22", "g22", min = -3, max = 3,value = -1,step=0.1),
                           br(),
                           sliderInput("g23", "g23", min = -3, max = 3,value = -2.8,step=0.1),
                           br(),
                           sliderInput("b2", "b2", min = 0, max = 3,value = 1.1,step=0.1),
                           br(),
                           sliderInput("h21", "h21", min = -3, max = 3,value = 0,step=0.1),
                           br(),
                           sliderInput("h22", "h22", min = -3, max = 3,value = -2,step=0.1),
                           br(),
                           sliderInput("h23", "h23", min = -3, max = 3,value = 2,step=0.1)
                  ),
                  tabPanel("Z",
                           br(),
                           sliderInput("a3", "a3", min = 0, max = 3,value = 0.1,step=0.1),
                           br(),
                           sliderInput("g31", "g31", min = -3, max = 3,value = 0.5,step=0.1),
                           br(),
                           sliderInput("g32", "g32", min = -3, max = 3,value = 2,step=0.1),
                           br(),
                           sliderInput("g33", "g33", min = -3, max = 3,value = 0.2,step=0.1),
                           br(),
                           sliderInput("b3", "b3", min = 0, max = 3,value = 0.9,step=0.1),
                           br(),
                           sliderInput("h31", "h31", min = -3, max = 3,value = 1,step=0.1),
                           br(),
                           sliderInput("h32", "h32", min = -3, max = 3,value = 0.1,step=0.1),
                           br(),
                           sliderInput("h33", "h33", min = -3, max = 3,value = 0,step=0.1)
                  )
                )
              ),
              mainPanel(
                fluidRow(
                  column(6,plotlyOutput("plotSSys")),
                  column(6,plotlyOutput("plotSSysPP"))
                ),
                withMathJax(),
                br(),
                "$$dX = \\alpha_1 \\cdot X^{g_{11}} \\cdot Y^{g_{12}} \\cdot Z^{g_{13}} - \\beta_1 \\cdot X^{h_{11}} \\cdot Y^{h_{12}} \\cdot Z^{h_{13}}$$",
                "$$dY = \\alpha_2 \\cdot X^{g_{21}} \\cdot Y^{g_{22}} \\cdot Z^{g_{23}} - \\beta_2 \\cdot X^{h_{21}} \\cdot Y^{h_{22}} \\cdot Z^{h_{23}}$$",
                "$$dZ = \\alpha_3 \\cdot X^{g_{31}} \\cdot Y^{g_{32}} \\cdot Z^{g_{33}} - \\beta_3 \\cdot X^{h_{31}} \\cdot Y^{h_{32}} \\cdot Z^{h_{33}}$$"
              )
            )
          )
  )
)

server <- function(input, output, session) {
  
  outLinear <- reactive({
    runLinear(input)
  })
  
  output$plotLinear <- renderPlotly({
    out <-  outLinear()
    
    p <- plot_ly(x = out[,"time"], y = out[,"X"], name = 'X', type = 'scatter', mode = 'lines',line = list(width = 1,color = "blue")) %>%
      add_lines(x = out[,"time"],y = out[,"Y"], name = 'Y', mode = 'lines',line = list(width = 1,color="green")) %>%
      add_lines(x = out[,"time"],y = out[,"Z"], name = 'Z', mode = 'lines',line = list(width = 1,color="red")) %>%
      layout(title="Time series",legend = list(x = 100, y = 0.5))
    
    p
  })
  
  output$plotLinPP <- renderPlotly({
    out <-  outLinear()
    
    p <- plot_ly(out, x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
                 line = list(width = 3, color = ~time, colorscale = 'Viridis')) %>%
    layout(title="Phase Plane",scene = list(camera = list(eye = list(x = 1.5, y = 2.25, z = 1.5))))
    
    p
  })
  
  outMA <- reactive({
    runMA(input)
  })
  
  output$plotMA <- renderPlotly({
    out <-  outMA()
    
    p <- plot_ly(x = out[,"time"], y = out[,"X"], name = 'X', type = 'scatter', mode = 'lines',line = list(width = 1,color = "blue")) %>%
      add_lines(x = out[,"time"],y = out[,"Y"], name = 'Y', mode = 'lines',line = list(width = 1,color="green")) %>%
      add_lines(x = out[,"time"],y = out[,"Z"], name = 'Z', mode = 'lines',line = list(width = 1,color="red")) %>%
      layout(title="Time series",legend = list(x = 100, y = 0.5))
    
    p
  })
  
  output$plotMAPP <- renderPlotly({
    out <-  outMA()
    
    p <- plot_ly(out, x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
                 line = list(width = 3, color = ~time, colorscale = 'Viridis')) %>%
      layout(title="Phase Plane",scene = list(camera = list(eye = list(x = 1.5, y = 2.25, z = 1.5))))
    
    p
  })
  
  outLV <- reactive({
    runLV(input)
  })
  
  output$plotLV <- renderPlotly({
    out <-  outLV()
    
    p <- plot_ly(x = out[,"time"], y = out[,"X"], name = 'X', type = 'scatter', mode = 'lines',line = list(width = 1,color = "blue")) %>%
      add_lines(x = out[,"time"],y = out[,"Y"], name = 'Y', mode = 'lines',line = list(width = 1,color="green")) %>%
      add_lines(x = out[,"time"],y = out[,"Z"], name = 'Z', mode = 'lines',line = list(width = 1,color="red")) %>%
      layout(title="Time series",legend = list(x = 100, y = 0.5))
    
    p
  })
  
  output$plotLVPP <- renderPlotly({
    out <-  outLV()
    
    p <- plot_ly(out, x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
                 line = list(width = 3, color = ~time, colorscale = 'Viridis')) %>%
      layout(title="Phase Plane",scene = list(camera = list(eye = list(x = 1.5, y = -2, z = 1.5))))
    
    p
  })
  
  outSSys <- reactive({
    runSSys(input)
  })
  
  output$plotSSys <- renderPlotly({
    out <-  outSSys()
    
    p <- plot_ly(x = out[,"time"], y = out[,"X"], name = 'X', type = 'scatter', mode = 'lines',line = list(width = 1,color = "blue")) %>%
      add_lines(x = out[,"time"],y = out[,"Y"], name = 'Y', mode = 'lines',line = list(width = 1,color="green")) %>%
      add_lines(x = out[,"time"],y = out[,"Z"], name = 'Z', mode = 'lines',line = list(width = 1,color="red")) %>%
      layout(title="Time series",legend = list(x = 100, y = 0.5))
    
    p
  })
  
  output$plotSSysPP <- renderPlotly({
    out <-  outSSys()
    
    p <- plot_ly(out, x = ~X, y = ~Y, z = ~Z, type = 'scatter3d', mode = 'lines',
                 line = list(width = 3, color = ~time, colorscale = 'Viridis')) %>%
      layout(title="Phase Plane",scene = list(camera = list(eye = list(x = 1.5, y = -2, z = 1.5))))
    
    p
  })
  
}
shinyApp(ui = ui, server = server)
