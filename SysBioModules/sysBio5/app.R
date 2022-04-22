library(deSolve)
library(shiny)

runNL1d <- function(input,x0){
  ## Chaos in the atmosphere
  SS <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dX <- d*(X - a)*(X - b)*(X - c)
      list(c(dX))
    })
  }
  parameters <- c(a = input$a, b = input$b,c = input$c, d = input$d)
  times <- seq(0,1, by = 0.01)
  
  out <- list()
  for(i in 1:length(x0)){
    state <- c(X = x0[i])
    out[[i]] <- ode(y = state, times = times, func = SS, parms = parameters)
  }
  
  return(out)
}

runSpat <- function(input,C){
  
  A <- rbind(c(input$a11,input$a12),c(input$a21,input$a22))
  #A <- t(A)
  
  traceA <- sum(diag(A))
  detA <- det(A)
  desA <- traceA^2 - 4*det(A)
  
  
  time <- seq(0, 10, 0.01)
  eigenA <- eigen(A)
  
  traj <- list()
  for(i in 1:nrow(C)){
    const <- solve(eigenA$vectors,C[i,])
    
    x1 <- eigenA$vectors[1,1]*const[1]*exp(eigenA$values[1]*time)+eigenA$vectors[1,2]*const[2]*exp(eigenA$values[2]*time)
    x2 <- eigenA$vectors[2,1]*const[1]*exp(eigenA$values[1]*time)+eigenA$vectors[2,2]*const[2]*exp(eigenA$values[2]*time)
    traj[[i]] <- data.frame(x1,x2)
    
  }
  
  return(list(C=C,traceA=traceA,detA=detA,desA=desA,traj=traj,eigenA=eigenA))
}

runEstab <- function(){

  A <- round(runif(4)*10-8)
  dim(A) <- c(2,2)
  
  traceA <- sum(diag(A))
  detA <- det(A)
  desA <- traceA^2 - 4*det(A)
  
  
  time <- seq(0, 10, 0.01)
  eigenA <- eigen(A)
  
  C <- cbind(c(10,10,5,-10,-10,-5),c(10,5,10,-10,-5,-10))
  
  traj <- list()
  for(i in 1:nrow(C)){
    const <- solve(eigenA$vectors,C[i,])
    
    x1 <- eigenA$vectors[1,1]*const[1]*exp(eigenA$values[1]*time)+eigenA$vectors[1,2]*const[2]*exp(eigenA$values[2]*time)
    x2 <- eigenA$vectors[2,1]*const[1]*exp(eigenA$values[1]*time)+eigenA$vectors[2,2]*const[2]*exp(eigenA$values[2]*time)
    traj[[i]] <- data.frame(x1,x2)
    
  }
  
  return(list(A=A,C=C,traceA=traceA,detA=detA,desA=desA,traj=traj,eigenA=eigenA))
}

zprintf <- function(z){
  if(is.complex(z)){
    isign <- c("+","")[(sign(Im(z))==-1)+1]
    sprintf('%0.3f %s %0.3fi', Re(z),isign, Im(z))  
  } else{
    sprintf('%0.3f', z)
  }
}

ui <- fluidPage(
  tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 1000px; /* or 950px */}")),
  titlePanel("Systems Biology - Module 5: Analysis of Dynamic Models"),
  tabsetPanel(
    tabPanel("Overview",
             br(),
             p("This is the fourth of the Biomedical Systems and Modeling modules that are intended to provide simple interfaces
                  to example systems models and to provide extra practice problems, and is the follow up to ",code("Module 4: Canonical Models")),
             p(code("Module 5: Analysis of Dynamic Models"),"provides interfaces to example models described in ",
               strong("Chapter 4: The Mathematics of Biological Systems"),"of",em("A First Course in Systems 
                  Biology."),"The expected outcomes of this module are to learn the fundamentals of stability analysis including: (i) linearization of nonlinear models, (ii) eigenvalue analysis."),
             tags$hr(),
             h4("Examples"),
             tags$ol(
               tags$li("1D-Nonlinear System"),
               tags$ul(
                 tags$li("Interface to a simple 1D-nonlinear system for exploring possible steady-states.")
               ),
               tags$li("Stability Patterns"),
               tags$ul(
                 tags$li("Interface for exploring stability patterns.")
                )
             ),
             h4("Futher Practice"),
             tags$ol(
               tags$li("Eigenvalues"),
               tags$ul(
                 tags$li("Every time this app is loaded a 2X2 matrix is randomly generated and the results of the 
                         corresponding eigenvalue analysis are revealed when prompted."
                 )
              ),
              tags$li("Linear Systems"),
              tags$ul(
                tags$li("Every time this app is loaded a random 2-variable linear system is generated and the results of the 
                         corresponding stability analysis are revealed when prompted."
                )
              )
             ),
             tags$hr(),
             h4("Citation"),
             p("Voit, E.O.:",em("A First Course in Systems Biology."),"Garland Science, New York, NY, 2017, 2nd edition"),
             p("**All page numbers reference the Second Edition**"),
             br()
               ),              
    tabPanel("1D-Nonlinear System",
             h3("1D-Nonlinear System"),
             tags$ul(
               tags$li("Interface to a simple 1D-nonlinear system for exploring possible steady-states.")
             ),
             br(),
             sidebarLayout(
               sidebarPanel(
                 br(),
                 sliderInput("nTraj", "Number of trajectories", min = 1, max = 25,value = 1, step = 1),
                 br(),
                 checkboxInput("showSSlines", label = "Show Steady States", value = F),
                 br(),
                 h4("Parameters"),
                 sliderInput("a", "a", min = 0, max = 5,value = 1, step = 1),
                 br(),
                 sliderInput("b", "b", min = 0, max = 5,value = 3, step = 1),
                 br(),
                 sliderInput("c", "c", min = 0, max = 5,value = 5, step = 1),
                 br(),
                 sliderInput("d", "d", min = -2, max = 2,value = -0.4, step = 0.1),
                 br()
               ),
               mainPanel(
                 fluidRow(
                   plotOutput("plotNL1d")
                 ),
                 withMathJax(),
                 br(),
                 "$$dX = d \\cdot (X - a) \\cdot (X - b) \\cdot (X - c)$$"
               )
             )
    ),
    tabPanel("Stability Patterns",
             h3("Stability Patterns and Trace-Determinant Plots"),
             tags$ul(
               tags$li("Interface for exploring stability patterns.")
             ),
             br(),
             sidebarLayout(
               sidebarPanel(
                 br(),
                 h4("Parameters"),
                 sliderInput("a11", "a11", min = -6, max = 6,value = -3, step = 1),
                 br(),
                 sliderInput("a12", "a12", min = -6, max = 6,value = -6, step = 1),
                 br(),
                 sliderInput("a21", "a21", min = -6, max = 6,value = 2, step = 1),
                 br(),
                 sliderInput("a22", "a22", min = -6, max = 6,value = 1, step = 1),
                 br()
               ),
               mainPanel(
                 fluidRow(
                   column(6,plotOutput("plotSpat")),
                   column(6,plotOutput("plotSpat2"))
                 ),
                 withMathJax(),
                 br(),
                 "$$ A = \\left( \\begin{array}{cc} a_{11} &  a_{12} \\\\ a_{21} &  a_{22} \\end{array} \\right)$$"
               )
             )
      ),
    tabPanel("Eigenvalues",
             h3("Eigenvalue Analysis"),
             tags$ul(
               tags$li("Every time this app is loaded a 2X2 matrix is randomly generated and the results of the 
                       corresponding eigenvalue analysis are revealed when prompted."
               )
             ),
             br(),
             
             fluidRow(column(12, align="center",htmlOutput("sysMat"),HTML("<br/><br/>"))),
             actionButton("showEigen", "Show Eigen Values"),
             h3("Trace"),
             verbatimTextOutput("trace"),
             h3("Determinant"),
             verbatimTextOutput("det"),
             h3("Discriminant"),
             verbatimTextOutput("des"),
             h3("Eigen values"),
             verbatimTextOutput("eigen"),
             br(),
             fluidRow(
               column(6,actionButton("showTrace", "Show Trace-Determinant"),plotOutput("plotEstab2")),
               column(6,actionButton("showPhase", "Show Phase Plane"),plotOutput("plotEstab"))
             ),
             br(),br(),br()
      ),
    tabPanel("Linear Systems",
             h3("Analysis of Linear Systems"),
             tags$ul(
               tags$li("Every time this app is loaded a random 2-variable linear system is generated and the results of the 
                       corresponding stability analysis are revealed when prompted."
               )
             ),
             br(),
             
             fluidRow(
               withMathJax(),
               column(5,uiOutput("eq1"),
                      uiOutput("eq2"),
                      actionButton("showSS", "Show Steady States"),
                      uiOutput("equSS"),
                      uiOutput("SS"),
                      actionButton("showPlot", "Show Vector Field"),
                      HTML("</br></br>The point is at the SS, and the color indicates stability: green - stable; magenta - unstable; cyan - saddle. Red and blue arrows are the eigenvectors."),
                      plotOutput("plotSys")
               ),
               column(1,""),
               column(6,actionButton("showJac", "Show Jacobian Matrix"),
                      uiOutput("jac"),
                      actionButton("showStab", "Show Stability Analysis"),
                      HTML("</br></br>"),
                      uiOutput("stab")
               )
             ),
             br(),br(),br()
        )
  )
)

server <- function(input, output, session) {
  x0 <- sample(seq(-6,6,0.5))
  
  sysNL1d <- reactive({
    runNL1d(input,x0)
  })
  
  output$plotNL1d <- renderPlot({
    out <- sysNL1d()
    plot(out[[1]][,1],out[[1]][,2],ylim=c(-6,6),xlim=c(0,1),xlab="Time",ylab="X",col=1,lwd=2,type="l")
    if(input$nTraj > 1){
      for(i in 1:(input$nTraj-1)){
        lines(out[[i+1]][,1],out[[i+1]][,2],col=1,lwd=2)
      }  
    }
    
    if(input$showSSlines){
      times <- seq(0,1, by = 0.01)
      lines(times,rep(input$a,length(times)),col=2,lwd=2)
      lines(times,rep(input$b,length(times)),col=3,lwd=2)
      lines(times,rep(input$c,length(times)),col=4,lwd=2)
      legend("topright",c("a","b","c"),lwd = 2, col=2:4)  
    }
    
  })
  
  C <- cbind(runif(10)*6-3,runif(10)*6-3)
  
  sysSpat <- reactive({
    runSpat(input,C)
  })
  
  output$plotSpat <- renderPlot({
    
    out <- sysSpat()
    
    x1r <- sapply(1:4, function(x) range(Re(out$traj[[x]]$x1)))
    x2r <- sapply(1:4, function(x) range(Re(out$traj[[x]]$x2)))
    
    plot(Re(out$traj[[1]]$x1),Re(out$traj[[1]]$x2),xlim=c(max(min(x1r[1,]),-20),min(max(x1r[2,]),20)),ylim=c(max(min(x2r[1,]),-20),min(max(x2r[2,]),20)),type="l",xlab="x1",ylab="x2",lwd=2)
    for(i in 2:length(out$traj)){
      lines(Re(out$traj[[i]]$x1),Re(out$traj[[i]]$x2),col=i,lwd=2)
    }
  })
  
  output$plotSpat2 <- renderPlot({
    out <- sysSpat()
    
    mt <- ifelse(out$traceA==0,1,abs(out$traceA))
    
    xax <- seq(-20,20,1)
    yax <- seq(min(out$detA-10,0),max((3*mt)^2/4,out$detA+10),1)
    
    plot(seq(-3*mt,3*mt,0.01),seq(-3*mt,3*mt,0.01)^2/4,xlab="tr A",ylab="det A",type="l",lwd=2,ylim=c(min(out$detA-10,0),max((3*mt)^2/4,out$detA+10)),cex.lab=1.25,cex.axis=1.25)
    #plot(seq(-10*mt,10*mt,0.01),seq(-10*mt,10*mt,0.01)^2 - 4*out$detA,xlab="tr A",ylab="det A",type="l",lwd=2)
    lines(xax,rep(0,length(xax)),lwd=3,col="darkgray")
    lines(rep(0,length(yax)),yax,lwd=3,col="darkgray")
    points(out$traceA,out$detA,pch=19,col=2,cex=2)
    
  })
  
  out <- runEstab()
  
  output$sysMat <- renderUI({
    HTML(paste("<br/><h3>Linear System</h3>","<h3>",
               paste(apply(out$A,1,function(x) paste(sprintf("%.0f",x),collapse="&nbsp;&nbsp;")),collapse = "<br/><br/>"),
               "</h3>"))
  })
  
  output$trace <- renderText({
    if(!input$showEigen) return()
    out$traceA
  })
  
  output$det <- renderText({
    if(!input$showEigen) return()
    out$detA
  })
  
  output$des <- renderText({
    if(!input$showEigen) return()
    out$desA
  })
  
  output$eigen <- renderText({
    if(!input$showEigen) return()
    out$eigenA$value
  })
  
  
  output$plotEstab <- renderPlot({
    if(!input$showPhase) return()
    x1r <- sapply(1:4, function(x) range(Re(out$traj[[x]]$x1)))
    x2r <- sapply(1:4, function(x) range(Re(out$traj[[x]]$x2)))
    
    plot(Re(out$traj[[1]]$x1),Re(out$traj[[1]]$x2),xlim=c(max(min(x1r[1,]),-20),min(max(x1r[2,]),20)),ylim=c(max(min(x2r[1,]),-20),min(max(x2r[2,]),20)),type="l",xlab="x1",ylab="x2",lwd=2,cex.lab=1.25,cex.axis=1.25)
    for(i in 2:length(out$traj)){
      lines(Re(out$traj[[i]]$x1),Re(out$traj[[i]]$x2),col=i,lwd=2)
    }
  })
  
  output$plotEstab2 <- renderPlot({
    if(!input$showTrace) return()
    mt <- ifelse(out$traceA==0,1,abs(out$traceA))
    plot(seq(-3*mt,3*mt,0.01),seq(-3*mt,3*mt,0.01)^2/4,xlab="tr A",ylab="det A",type="l",lwd=2,ylim=c(min(out$detA-10,0),max((3*mt)^2/4,out$detA+10)),cex.lab=1.25,cex.axis=1.25)
    #plot(seq(-10*mt,10*mt,0.01),seq(-10*mt,10*mt,0.01)^2 - 4*out$detA,xlab="tr A",ylab="det A",type="l",lwd=2,ylim=c(min(out$detA-10,0),(10*mt)^2),cex.lab=1.25,cex.axis=1.25)
    points(out$traceA,out$detA,pch=19,col=2,cex=2)
    
  })
  
  sSys <- function(x,r) {
    dX1 <- r[1] + r[2]*x[1] + r[3]*x[2] 
    dX2 <- r[4] + r[5]*x[1] + r[6]*x[2] 
    c(dX1,dX2)
  }
  
  jac <- function(x,r){
    J <- rbind(c(r[2], r[3]),
               c(r[5],r[6]))
    
    return(J)
  }
  
  kr <- 5
  
  r <- round(runif(6)*2*kr-kr)
  
  ssA <- rbind(c(r[2],r[3]),c(r[5],r[6]))
  ssB <- c(-r[1],-r[4])
  while(sum(r[-c(1,4)])==0|det(ssA)==0){
    r <- round(runif(6)*2*kr-kr)
    ssA <- rbind(c(r[2],r[3]),c(r[5],r[6]))
    ssB <- c(-r[1],-r[4])
  }
  
  SS <- solve(ssA,ssB)
  SS[abs(SS)<1e-9] <- 0
  ssOrigin <- T #sum(!is.finite(sSys(t(c(0,0)),r,g)))>0
  
  J1 <- jac(SS,r)
  J1[abs(J1) < 1e-9] <- 0
  eig1 <- eigen(J1)
  eig1$values[abs(eig1$values) < 1e-9] <- 0
  
  sterms <- c(ifelse(r[1]==0,"",as.character(r[1])),
              ifelse(r[2]==0,"",ifelse(r[2]==1,"X",paste(r[2],"X"))),
              ifelse(r[3]==0,"",ifelse(r[3]==1,"Y",paste(r[3],"Y"))),
              ifelse(r[4]==0,"",as.character(r[4])),
              ifelse(r[5]==0,"",ifelse(r[5]==1,"X",paste(r[5],"X"))),
              ifelse(r[6]==0,"",ifelse(r[6]==1,"Y",paste(r[6],"Y"))))
  
  jequs <- c(as.character(r[2]),
             as.character(r[3]),
             as.character(r[5]),
             as.character(r[6]))
  jequs[jequs==""] <- "0"
  
  output$eq1 <- renderUI({ 
    withMathJax(paste0("$$X' = ",sterms[1],ifelse(sterms[2]=="","",c("","+")[(r[2]>0)+1]),sterms[2],ifelse(sterms[3]=="","",c("","+")[(r[3]>0)+1]),sterms[3],"$$"))
  })
  
  output$eq2 <- renderUI({ 
    withMathJax(paste0("$$Y' = ",sterms[4],ifelse(sterms[5]=="","",c("","+")[(r[5]>0)+1]),sterms[5],ifelse(sterms[6]=="","",c("","+")[(r[6]>0)+1]),sterms[6],"$$"))
  })
  
  observeEvent(input$showJac, {
    show("jac")
  })
  
  output$jac <- renderUI({
    if(!input$showJac) return()
    withMathJax(paste("$$J=\\left( \\begin{array}{cc}",
                      jequs[1],"&",jequs[2],"\\\\",
                      jequs[3],"&",jequs[4],
                      "\\end{array} \\right)$$"))
  })
  
  observeEvent(input$showSS, {
    show("SS")
  })
  
  output$SS <- renderUI({
    if(!input$showSS) return()
    withMathJax(sprintf("$$X=%0.3f, Y=%0.3f$$",SS[1],SS[2]))
  })
  
  output$equSS <- renderUI({
    if(!input$showSS) return()
    withMathJax(paste("$$ \\left( \\begin{array}{cc}",
                      r[2],"&",r[3],"\\\\",
                      r[5],"&",r[6],
                      "\\end{array} \\right) \\left( \\begin{array}{c} X \\\\ Y \\end{array} \\right) = ",
                      "\\left( \\begin{array}{c}",-r[1],"\\\\",-r[4],"\\end{array} \\right) $$"))
  })
  
  observeEvent(input$showStab, {
    show("stab")
  })
  
  output$stab <- renderUI({
    if(!input$showStab) return()
    withMathJax(paste(sprintf("SS1 $$X=%0.3f, Y=%0.3f \\\\ J=\\left( \\begin{array}{cc}",SS[1],SS[2]),
                      zprintf(J1[1,1]),"&",zprintf(J1[1,2]),"\\\\",
                      zprintf(J1[2,1]),"&",zprintf(J1[2,2]),
                      "\\end{array} \\right) \\\\",
                      sprintf("\\lambda = %s, %s $$",zprintf(eig1$values[1]),zprintf(eig1$values[2]))))
  })
  
  observeEvent(input$showGraph, {
    show("plotSys")
  })
  
  output$plotSys <- renderPlot({
    if(!input$showPlot) return()
    
    X1 <- seq(-15,15,length.out = 20)
    X2 <- seq(-15,15,length.out = 20)
    x <- expand.grid(X1,X2)
    
    par(cex.lab=1.5,cex.axis=1.5)
    plot(NULL,NULL,type = "l",las = 1,xlim = range(X1),ylim = range(X2),axes=F,col=2,lwd=5,xlab="X",ylab="Y")
    #lines(X1,nc2,col=4,lwd=5)
    
    evs <- 2
    #EV <- rbind(c(SS,Re(eig1$vectors[,1]*eig1$values[1])+SS),c(SS,Re(eig1$vectors[,2]*eig1$values[2])+SS))
    EV <- rbind(c(SS,Re(eig1$vectors[,1]*eig1$values[1])*evs+SS),c(SS,Re(eig1$vectors[,2]*eig1$values[2])*evs+SS))
    EV[Re(eig1$values)<0,] <- cbind(EV[Re(eig1$values)<0,3:4],EV[Re(eig1$values)<0,1:2])
    
    dX <- t(apply(x,1,function(xx) sSys(xx,r)))
    arrows(x[,1], x[,2],(dX[,1]/20+x[,1]),(dX[,2]/20+x[,2]),length=0.08,col="dark gray",lwd=2)
    axis(1, pos=0,lwd=3)
    axis(2, pos=0,lwd=3)
    points(SS[1],SS[2],lwd=4,cex=2,col=c("magenta","cyan","green")[sum(Re(eig1$value)<0)+1])
    arrows(EV[,1], EV[,2],EV[,3],EV[,4],length=0.12,col=c(2,4),lwd=3)
    
  })
  
}
shinyApp(ui = ui, server = server)
