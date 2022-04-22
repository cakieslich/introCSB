library(deSolve)
library(shiny)
library(shinyjs)
library(igraph)

runSIR <- function(input){
  ## Chaos in the atmosphere
  SIR <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- rBirth + rS*R - rI*S*I 
      dI <- rI*S*I - rR*I - rDeath*I
      dR <- rR*I - rS*R
      list(c(dS, dI, dR))
    })
  }
  parameters <- c(rBirth = 3, rDeath = 0.02,rI = 0.0005, rR = 0.05, rS = 0.01)
  state <- c(S = 1000 - input$I0, I = input$I0, R = 0)
  times <- seq(0,100, by = 1)
  out <- ode(y = state, times = times, func = SIR, parms = parameters)
  
  out2 <- data.frame(S = 1000 - input$I0, I = input$I0, R = 0)
  timeD <- seq(0,100,input$timeStep) 
  while(nrow(out2)<length(timeD)){
    out2 <- rbind(out2,with(tail(out2,1),
                            data.frame(S=S+input$rBirth+input$rS*R-input$rI*S*I,
                                       I=I+input$rI*S*I-input$rDeath*I-input$rR*I,
                                       R=R+input$rR*I-input$rS*R)
                            )
                  )
  }
  out2 <- cbind(time=timeD,out2)
  return(list(out,out2))
}

runPop <- function(input){
  ## Chaos in the atmosphere
  Pop <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dP <- rG*P
      list(c(dP))
    })
  }
  parameters <- c(rG = input$rG)
  state <- c(P = 1)
  times <- seq(0,input$nIter, by = 0.1)
  out <- ode(y = state, times = times, func = Pop, parms = parameters)
  
  times <- seq(0,input$nIter, by = 1)
  out2 <- data.frame(P = 1)
  while(nrow(out2)<length(times)){
    out2 <- rbind(out2,with(tail(out2,1),
                            data.frame(P=2*P)
                            )
                  )
  }
  out2 <- cbind(time=times,out2)
  return(list(out,out2))
}

runRBC <- function(input){
  
  out <- data.frame(R = input$rbcR0, M = input$rbcM0)
  timeD <- seq(0,20,1) 
  while(nrow(out)<length(timeD)){
    out <- rbind(out,with(tail(out,1),
                            data.frame(R=(1-input$rbcf)*R+M,
                                       M=input$rbcg*input$rbcf*R)
    )
    )
  }
  out <- cbind(time=timeD,out)
  return(out)
}

runLeslie <- function(input){
  
  out <- data.frame(P1 = input$P10, P2 = input$P20, P3 = input$P30, P4 = input$P40)
  timeD <- seq(0,10,1) 
  while(nrow(out)<length(timeD)){
    out <- rbind(out,with(tail(out,1),
                          data.frame(P1=input$alpha1*P1+input$alpha2*P2+input$alpha3*P3+input$alpha4*P4,
                                     P2=input$sigma1*P1,
                                     P3=input$sigma2*P2,
                                     P4=input$sigma3*P3)
                          )
                )
  }
  
  out <- cbind(time=timeD,out)
  return(out)
}

ui <- fluidPage(
    tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 1000px; /* or 950px */}; .checkbox{ margin-top: 0px; }")),
    titlePanel("Biomedical Systems and Modeling - Module 3: Discrete Models"),
    tabsetPanel(
      tabPanel("Overview",
               br(),
               p("This is the third of the Biomedical Systems and Modeling modules that are intended to provide simple interfaces
                  to example systems models and to provide extra practice problems, and is the follow up to ",code("Module 2: Static Network Models")),
               p(code("Module 3: Discrete, Recursive Models"),"provides interfaces to example discrete and recursive models, including examples from",
                  strong("Chapter 4: The Mathematics of Biological Systems"),"of",em("A First Course in Systems 
                  Biology."),"The expected outcomes of this module are  (i) learn the basic principles of discrete, recursive models; (ii) 
                  explore examples of dynamic discrete recrusive models; (iii) practice generating Markov matrices."),
               tags$hr(),
               h4("Examples"),
               tags$ol(
                 tags$li("Cell Growth Models (p. 85)"),
                 tags$ul(
                   tags$li("Interface for comparing simple models for cell growth.")
                 ),
                 tags$li("SIR Models Revisited"),
                 tags$ul(
                   tags$li("A comparison of continuous and discrete SIR models that predict the spread of an infectious 
                            disease."
                   )
                 ),
                 tags$li("Red Blood Cell Dynamics (p. 86)"),
                 tags$ul(
                   tags$li("A simple interface for a discrete, recursive model of red blood cell dynamics."
                   )
                 ),
                  tags$li("Growth of Age-structured Populations (p. 289)"),
                 tags$ul(
                   tags$li("Interface to a simple model for the growth of a population stratified by age 
                            (equal age brackets)."
                   )
                 )
               ),
               h4("Futher Practice"),
               tags$ol(
                 tags$li("Markov Matrices (p. 88)"),
                 tags$ul(
                   tags$li("Every time this app is loaded a new 3 node Markov model is generated. For each Markov model, the app 
                           generates a plot of the graph with the edges labeled by the probability of a given transition and the 
                           corresponding Markov. The graph and Markov matrix only appear after clicking the associated button."
                   )
                 )
               ),
               tags$hr(),
               h4("Citation"),
               p("Voit, E.O.:",em("A First Course in Systems Biology."),"Garland Science, New York, NY, 2017, 2nd edition"),
               p("**All page numbers reference the Second Edition**"),
               br()
      ),              
      tabPanel("Cell Growth",
               h3("Cell Growth Models"),
               tags$ul(
                 tags$li("Interface for comparing simple models for cell growth."),
                 tags$li("Identify the continuous growth rate (rG) that results in the best agreement between 
                         the continuous and discrete models."),
                 tags$li("Is there any significance to this value?")
                 ),
               br(),
               sidebarLayout(
                 sidebarPanel(
                   sliderInput("nIter", "Max. # of divsions", min = 0, max = 20,value = 6, step = 1),
                   br(),
                   sliderInput("rG", "Growth rate (rG)", min = 0.1, max = 2,value = 2, step = 0.01)
                 ),
                 mainPanel(
                   plotOutput("plotPop",height = 350),
                   br(),
                   withMathJax(),
                   fluidRow(
                     column(6,h4("ODE model",style="text-align:center;"),
                            "$$dP = r_G \\cdot  P$$"
                     ),
                     column(6,h4("Recursive model",style="text-align:center;"),
                            "$$P_{t+\\tau} = 2 \\cdot  P_t$$"
                     )
                   )
                 )
               )
      ),
      tabPanel("SIR Models",
               h3("SIR Model"),
               tags$ul(
                 tags$li("A comparison of continuous and discrete SIR models that predict the spread of an infectious 
                          disease. The interface allows for manipulation of the parameters of the discrete model, while 
                          the parameters of the dynamic model
                         are fixed to the initial values assigned for the discrete model. Agreement between the discrete and 
                         continuous models is reported in terms of average percent error."),
                 tags$li("What are the differences in the equations for the two models?"),
                 tags$li("Explore the effects of the discrete time step (tau) on the agreement between the models."),
                 tags$li("Do you need to change the parameters of the discrete model along with time step to maximize 
                         agreement between the models?"),
                 tags$li("What changes to the model provide the best agreement (minimal error) between the models?")
                 ),
               br(),
               sidebarLayout(
                 sidebarPanel(
                   tabsetPanel(
                     tabPanel("Time/Init. Cond.",
                              br(),
                              sliderInput("timeStep", "Time Step", min = 0.1, max = 10,value = 10, step = 0.1),
                              br(),
                              sliderInput("I0", "Initial # of Infected Individuals (I0)", min = 0, max = 1000,value = 10, step = 1),
                              br()
                     ),
                     tabPanel("Parameters",
                              br(),
                              sliderInput("rBirth", "Birth Rate (rBirth)", min = 0.3, max = 30,value = 3, step = 0.3),
                              br(),
                              sliderInput("rDeath", "Death Rate (rDeath)", min = 0.002, max = 0.2,value = 0.02, step = 0.002),
                              br(),
                              sliderInput("rI", "Infection Rate (rI)", min = 0.00005, max = 0.005,value = 0.0005, step = 0.00005),
                              br(),
                              sliderInput("rR", "Immunity Rate (rR)", min = 0.005, max = 0.5,value = 0.05, step = 0.005),
                              br(),
                              sliderInput("rS", "Susceptibility Rate (rS)", min = 0.001, max = 0.1,value = 0.01, step = 0.001)
                     )
                   )
                 ),
                 mainPanel(
                   plotOutput("plotSIR",height = 350),
                   br(),
                   fluidRow(
                     column(6,h4("ODE model",style="text-align:center;"),
                            "$$dS = r_{Birth} + r_S \\cdot  R - r_I \\cdot  S \\cdot  I$$",
                            "$$dI = r_I \\cdot  S \\cdot  I - r_R \\cdot  I - r_{Death} \\cdot  I$$",
                            "$$dR = r_R \\cdot  I - r_S \\cdot  R$$"
                            ),
                     column(6,h4("Recursive model",style="text-align:center;"),
                            "$$S_{t+\\tau} = S_t + r_{Birth} + r_S \\cdot  R_t - r_I \\cdot  S_t \\cdot  I_t$$",
                            "$$I_{t+\\tau} = I_t + r_I \\cdot  S_t \\cdot  I_t - r_R \\cdot  I - r_{Death} \\cdot  I_t$$",
                            "$$R_{t+\\tau} = R_t + r_R \\cdot  I_t - r_S \\cdot  R_t$$"
                     )
                   )
                 )
               )
            ),
      tabPanel("Red Blood Cells",
               h3("Red Blood Cell Dynamics"),
               tags$ul(
                 tags$li("A simple interface for a discrete, recursive model of red blood cell (RBC) dynamics. The
                         model assumes that RBCs are constantly formed by the bone marrow and that RBCs are constantly 
                         destroyed by the spleen. In normal conditions, the total number of RBCs is approximately constant from 
                         day to day. Rn represents the number of red blood cells in circulation on day n, while Mn is 
                         the number of red blood cells produced by the bone marrow on day n."),
                 tags$li("Explore the effects of initial conditions and the parameter values."),
                 tags$li("What is the meaning of g lesser/greater than 1? What is the meaning of g equals 1?"),
                 tags$li("Can you find a steady state for this system? What determines to observed steady state 
                         values of Rn anf Mn?")
                 ),
               br(),
               sidebarLayout(
                 sidebarPanel(
                   sliderInput("rbcR0", "Initial # of RBCs in circulation", min = 0, max = 1500,value = 1000),
                   br(),
                   sliderInput("rbcM0", "Initial # of RBCs produced by marrow", min = 0, max = 20,value = 10, step = 1),
                   br(),
                   sliderInput("rbcf", "Fraction removed by spleen", min = 0, max = 0.2,value = 0.01, step = 0.01),
                   br(),
                   sliderInput("rbcg", "Relative growth rate", min = 0, max = 3,value = 2, step = 0.1)
                 ),
                 mainPanel(
                   fluidRow(
                     column(6,plotOutput("plotRBCR",height = 350)
                     ),
                     column(6,plotOutput("plotRBCM",height = 350)
                     )
                   ),
                   br(),
                   h4("Model Equations",style="text-align:center;"),
                            "$$R_n = (1-f) \\cdot  R_{n-1} + M_{n-1}$$",
                            "$$M_n = g \\cdot f \\cdot  R_{n-1}$$"
                 )
               )
          ),
      tabPanel("Age-structured Populations",
               h3("Growth of Age-structured Populations"),
               tags$ul(
                 tags$li("Interface to a simple model for the growth of a population stratified by age 
                         (equal age brackets). Age groups have specific proliferation rates, alpha, and all 
                         individuals have the same potential maximal life span, but the survival rates (sigma)  
                         change with age. This model assumes that the maximum age is 120 years. The model 
                          contains 4 age classes, and as a result each iteration represents 30 years."),
                 tags$li("Explore the effects of the sruvival and proliferation rates, as well as, the initial 
                         populations on the population size."
                         )
               ),
               br(),
               sidebarLayout(
                 sidebarPanel(
                   tabsetPanel(
                     tabPanel("Init. Pop.",
                       br(),
                       sliderInput("P10", "Initial # of people aged 0 to 30", min = 0, max = 20,value = 8, step = 1),
                       br(),
                       sliderInput("P20", "Initial # of people aged 31 to 60", min = 0, max = 20,value = 4, step = 1),
                       br(),
                       sliderInput("P30", "Initial # of people aged 61 to 90", min = 0, max = 20,value = 2, step = 1),
                       br(),
                       sliderInput("P40", "Initial # of people aged 91 to 120", min = 0, max = 20,value = 1, step = 1)
                     ),
                     tabPanel("Parameters",
                       br(),
                       sliderInput("alpha1", "Proliferation rate for ages 0 to 30", min = 0, max = 1,value = 0.7, step = 0.05),
                       br(),
                       sliderInput("sigma1", "Survival rate for ages 0 to 30", min = 0, max = 1,value = 0.95, step = 0.05),
                       br(),
                       sliderInput("alpha2", "Proliferation rate for ages 31 to 60", min = 0, max = 1,value = 0.4, step = 0.05),
                       br(),
                       sliderInput("sigma2", "Survival rate for ages 31 to 60", min = 0, max = 1,value = 0.9, step = 0.05),
                       br(),
                       sliderInput("alpha3", "Proliferation rate for ages 61 to 90", min = 0, max = 1,value = 0, step = 0.05),
                       br(),
                       sliderInput("sigma3", "Survival rate for ages 61 to 90", min = 0, max = 1,value = 0.1, step = 0.05),
                       br(),
                       sliderInput("alpha4", "Proliferation rate for ages 91 to 120", min = 0, max = 1,value = 0, step = 0.05),
                       br()
                     )
                   )
                 ),
                 mainPanel(
                   fluidRow(
                     column(6,plotOutput("plotAgePop",height = 350)
                     ),
                     column(6,plotOutput("plotTotPop",height = 350)
                     )
                   ),
                   br(),
                   h4("Model Equations",style="text-align:center;"),
                   "$$P_{1,t+\\tau} = \\alpha_1 \\cdot  P_{1,t} + \\alpha_2 \\cdot  P_{2,t} + \\alpha_3 \\cdot  P_{3,t} + \\alpha_4 \\cdot  P_{4,t}$$",
                   "$$P_{2,t+\\tau} = \\sigma_1 \\cdot  P_{1,t}$$",
                   "$$P_{3,t+\\tau} = \\sigma_2 \\cdot  P_{2,t}$$",
                   "$$P_{4,t+\\tau} = \\sigma_3 \\cdot  P_{3,t}$$"
                 )
               )
            ),
      tabPanel("Markov Matrices",
               h3("Markov Matrices"),
               tags$ul(
                 tags$li("Every time this app is loaded a new 3 node Markov model is generated. For each Markov model, the app 
                          generates a plot of the graph with the edges labeled by the probability of a given transition and the 
                          corresponding Markov. The graph and Markov matrix only appear after clicking the associated button."),
                 tags$li("Reveal the graph. Determine the Markov matrix and check your answer on completion. Compute powers 
                         of the Markov matrix assuming the initial position is node 2. Compute M^2 and M^3, then compute M^4 and M^6, 
                         as well as higher powers. Document results. Do the matrices and the vector converge? What do the results 
                         mean? Are there absorbing states in the Markov model?"),
                 tags$li("Reveal the Markov matrix. Draw a graph with the appropriate edge weights and check your answer on completion.")
                 ),
               br(),
               fluidRow(
                 column(6,actionButton("showGraph", "Show Graph"),
                        plotOutput("plotSys")
                 ),
                 column(1,""),
                 column(4,actionButton("showMat", "Show Markov Matrix"),
                        htmlOutput("adjMat"))
               ),
               fluidRow(HTML("</br></br>"))
      )
    )
)

server <- function(input, output, session) {
  
  outSIR <- reactive({
    runSIR(input)
  })
  
  output$plotSIR <- renderPlot({
    out <- outSIR()
    
    par(mar=c(4,4,2,1))
    plot(out[[1]][,1],out[[1]][,2],ylim=c(0,1000),xlab="Time",ylab="Number of individuals",col=4,lwd=2,type="l",cex.lab=1.3,cex.axis=1.3)
    lines(out[[1]][,1],out[[1]][,3],xlab="Time",ylab="",col=2,lwd=2,type="l")
    lines(out[[1]][,1],out[[1]][,4],xlab="Time",ylab="",col=3,lwd=2,type="l")
    
    points(out[[2]][,1],out[[2]][,2],col=4,lwd=2)
    points(out[[2]][,1],out[[2]][,3],xlab="Time",ylab="",col=2,lwd=2)
    points(out[[2]][,1],out[[2]][,4],xlab="Time",ylab="",col=3,lwd=2)
    legend("topright",c("S","I","R"),lwd = 2, col=c(4,2,3))
    error <- abs(out[[1]][out[[1]][,1]%in%out[[2]][,1],-1]-out[[2]][out[[2]][,1]%in%out[[1]][,1],-1])/out[[1]][out[[1]][,1]%in%out[[2]][,1],-1]
    mtext(sprintf("Average Percent Error: %.0f%%",100*mean(unlist(error),na.rm = T)),cex=1.5)
    
  })
  
  outPop <- reactive({
    runPop(input)
  })
  
  output$plotPop <- renderPlot({
    out <- outPop()
    
    par(mar=c(4,4,2,1))
    plot(out[[1]][,1],out[[1]][,2],ylim=c(0,max(out[[2]][,2])),xlab="Time",ylab="Number of cells",col=4,lwd=2,type="l",cex.lab=1.3,cex.axis=1.3)
    points(out[[2]][,1],out[[2]][,2],col=4,lwd=2,pch=19)
    
  })
  
  outRBC <- reactive({
    runRBC(input)
  })
  
  output$plotRBCR <- renderPlot({
    out <- outRBC()
    
    par(mar=c(4,4,2,1))
    plot(out$time,out$R,xlab="Number of iterations",ylab="Number of RBCs in circulation",col=4,lwd=2,type="p",cex.lab=1.3,cex.axis=1.3,pch=19)
    
  })
  
  output$plotRBCM <- renderPlot({
    out <- outRBC()
    
    par(mar=c(4,4,2,1))
    plot(out$time,out$M,xlab="Number of iterations",ylab="Number of RBCs produced in marrow",col=2,lwd=2,type="p",cex.lab=1.3,cex.axis=1.3,pch=19)
    
  })
  
  outLeslie <- reactive({
    runLeslie(input)
  })
  
  output$plotAgePop <- renderPlot({
    out <- outLeslie()
    
    par(mar=c(4,4,2,1))
    barplot(t(out[,-1]),ylab="Population size",col=1:4,lwd=2,cex.lab=1.3,
            cex.axis=1.3,beside = T,legend=c("0-30 yrs","31-60 yrs","61-90 yrs","91-120 yrs"),
            args.legend = list(x="topleft"),names.arg=as.character(seq(0,300,30)))
    
  })
  
  output$plotTotPop <- renderPlot({
    out <- outLeslie()
    
    par(mar=c(4,4,2,1))
    plot(out$time*30,apply(out[,-1],1,sum),xlab="Years",ylab="Total population size",col="dimgray",lwd=2,type="p",cex.lab=1.3,cex.axis=1.3,pch=19)
    
  })
  
  m <<- 3
  
  mm1 <- round(runif(m),digits = 2)
  mm2 <- round((1-mm1)*runif(m),digits = 2)
  mm3 <- 1-mm1-mm2
  mm <- cbind(mm1,mm2,mm3)
  
  output$adjMat <- renderUI({
    if(!input$showMat) return()     
    HTML(paste("<h3>Markov Matrix</h3>",
               paste("<h3>",
                     paste(apply(mm,1,function(x) paste(sprintf("%.2f",x),collapse="&nbsp;&nbsp;")),collapse = "<br/><br/>"),
                     "</h3>"))
    )
  })
  
  observeEvent(input$showGraph, {
    #show("plotSys")
  })
  
  observeEvent(input$showMat, {
    #show("adjMat")
  })
  
  am <- sapply(1:m,function(x) rep(1,m))
  
  output$plotSys <- renderPlot({
    if(!input$showGraph) return()      
    net.bg <- graph_from_adjacency_matrix(t(am))
    V(net.bg)$size <- 32
    V(net.bg)$frame.color <- "black"
    V(net.bg)$color <- 1:m
    #V(net.bg)$label <- "" 
    E(net.bg)$arrow.mode <- 1
    E(net.bg)$label <- sprintf("%.2f",as.vector(mm))
    #l <- layout_in_circle(net.bg)
    l <- layout.fruchterman.reingold(net.bg)
    plot(net.bg,edge.label.color="black",edge.arrow.size	= 0.6,edge.curved=rep(0.2, m^2),layout=l,edge.label.cex=1.4,vertex.label.cex=1.4)
  })
  
  
}
shinyApp(ui = ui, server = server)