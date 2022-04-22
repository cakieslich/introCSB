library(shiny)
library(igraph)
library(pracma)

runDegDist <- function(input){
  ## Chaos in the atmosphere
  
  m <- 10
  
  edges <- NULL
  for(i in 1:(m-1)){
    for(id in paste0(letters[(i+1):m],i)){
      #print(id)
      edges <- c(edges,input[[id]])
    }
  }
  edges <- edges+0
  
  adjMat <- matrix(0,nrow=m,ncol=m)
  adjMat[lower.tri(adjMat)] <- edges
  adjMat[upper.tri(adjMat)] <- t(adjMat)[upper.tri(adjMat)]
  #adjMat <- t(adjMat)
  
  return(adjMat)
}

genStoichMat <- function(am,nc){
  sm <- array(0,dim = c(nc,sum(am)))
  ncc <- round(runif(nc)*2)+1
  
  nn <- 1
  for(n in 1:(2*nc)){
    outn <- which(am[n,]==1)
    if(length(outn)>0){
      for(on in outn){
        if(n <= nc & on <= nc){
          if(ncc[n] == ncc[on]){
            sm[on,nn] <- 1
            sm[n,nn] <- -1
            nn <- nn + 1
          } else{
            sm[on,nn] <- ncc[n]
            sm[n,nn] <- -ncc[on]
            nn <- nn + 1
          }
        } else{
          if(on > nc){
            sm[n,nn] <- -1
            nn <- nn + 1
          } else if(n > nc){
            sm[on,nn] <- 1
            nn <- nn + 1
          }
        }
      }
    }
  }
  
  smHTML <- sapply(1:nc,function(x) paste("<font color=\"#515151\">",x,"</font>&nbsp;&nbsp;",paste(sm[x,],collapse="&nbsp;&nbsp;")))
  smHTML <- gsub("&nbsp;&nbsp;-1","&nbsp;-1",smHTML)
  smHTML <- gsub("&nbsp;&nbsp;-2","&nbsp;-2",smHTML)
  smHTML <- gsub("&nbsp;&nbsp;-3","&nbsp;-3",smHTML)
  smHTML <- gsub("; 3",";&nbsp;&nbsp; 3",smHTML)
  smHTML <- gsub("; 2",";&nbsp;&nbsp; 2",smHTML)
  smHTML <- gsub("; 1",";&nbsp;&nbsp; 1",smHTML)
  smHTML <- gsub("; 0",";&nbsp;&nbsp; 0",smHTML)
  
  return(list(sm=sm,html=smHTML,ncc=ncc))
}

runStoichSys <- function(input,selF,sm){
  nc <- 5
  
  mfv <- NULL
  for(m in 1:8){
    eval(parse(text=paste0("mfv <- c(mfv,input$flux",m,")")))
  }
  mfv <- as.numeric(mfv)
  
  of <- c(6,8,10,12,13)
  ford <- c(of,(1:13)[-of][mfv>0])

  fo <- c(selF,ford[ford!=selF])
  
  rrm <- rref(sm[,fo])
  cfi <- which(apply(abs(rrm),2,sum)==1 & apply(rrm,2,sum)==1)[1:5]
  
  cf <- fo[cfi]
  mf <- fo[!(fo%in%cf)]
  
  cfv <- -rrm[,-cfi]%*%mfv[mfv>0]
  
  fvals <- rep(0,13)
  for(i in 1:13){
    if(sum(cf==i)>0){
      fvals[i] <- cfv[cf==i]
    }
    
    if(sum(mf==i)>0){
      fvals[i] <- mfv[mfv>0][mf==i]
    }
  }
  
  return(fvals)
}

runStoichSysP <- function(sm){
  nc <- 5
  
  mfv <- round(runif(ncol(sm)-5)*9) + 1 
  
  of <- sample(1:ncol(sm),5)
  fo <- c(of,(1:ncol(sm))[-of])
  
  rrm <- rref(sm[,fo])
  cfi <- which(apply(abs(rrm),2,sum)==1 & apply(rrm,2,sum)==1)[1:5]
  
  cf <- fo[cfi]
  mf <- fo[!(fo%in%cf)]
  
  cfv <- -rrm[,-cfi]%*%mfv
  
  count <- 0
  while(sum(cfv<0)>0 & count < 100){
    mfv <- round(runif(ncol(sm)-5)*9) + 1
    
    of <- sample(1:ncol(sm),5)
    fo <- c(of,(1:ncol(sm))[-of])
    
    rrm <- rref(sm[,fo])
    cfi <- which(apply(abs(rrm),2,sum)==1 & apply(rrm,2,sum)==1)[1:5]
    
    cf <- fo[cfi]
    mf <- fo[!(fo%in%cf)]
    
    cfv <- -rrm[,-cfi]%*%mfv
    
    count <- count + 1
  }

  fvals <- rep(0,ncol(sm))
  for(i in 1:ncol(sm)){
    if(sum(cf==i)>0){
      fvals[i] <- cfv[cf==i]
    }
    
    if(sum(mf==i)>0){
      fvals[i] <- mfv[mf==i]
    }
  }
  
  fType <- rep("Measured",ncol(sm))
  fType[cf] <- "Unknown"
  
  ftab <- data.frame(Flux=1:ncol(sm),Value=fvals,Type=fType)
  ftab <- ftab[order(ftab$Flux),]
  
  return(ftab)
}

getStoichNet <- function(){
  nc <- 5
  nflux <- 12
  
  inf <- sample(1:nc,round(nc*0.4))
  outf <- (1:nc)[-inf]
  
  AM <- array(0,dim=c(nc*2,nc*2))
  
  tmp <- rep(0,nc^2)
  tmp[sample(1:nc^2,round(1/3*nc^2))] <- 1 
  if(sum(tmp)>nflux){
    tmp[sample(which(tmp==1),sum(tmp)-nflux)] <- 0
  }
  
  dim(tmp) <- c(nc,nc)
  diag(tmp) <- rep(0,nc)
  AM[1:nc,1:nc] <- tmp
  
  nn2 <- 1
  for(i in inf){
    AM[nc+nn2,i] <- 1
    nn2 <- nn2 + 1
  }
  
  for(o in outf){
    AM[o,nc+nn2] <- 1
    nn2 <- nn2 + 1
  }
  
  noIn <- which(apply(AM[,1:nc],2,sum)==0)
  for(ni in noIn){
    sn <- sample(which(AM[ni,]==0),1)
    AM[ni,sn] <- 0
    AM[sn,ni] <- 1
  }
  
  for(j in 1:nc){ AM[j,(AM[j,1:nc]+AM[1:nc,j])==2] <- 0;}
  AM[(1:(2*nc)^2)[(1:(2*nc)+(1:(2*nc)-1)*2*nc)]] <- 0
  
  for(n in 1:nc){
    if(sum(AM[,n])==0){
      AM[sample((1:(2*nc))[-c(n,which(AM[n,]==1))],1),n] <- 1
    }
    if(sum(AM[n,])==0){
      AM[n,sample((1:(2*nc))[-c(n,which(AM[,n]==1))],1)] <- 1
    }
  }
  
  while(sum(AM)<nflux){
    pn <- as.vector(sapply(1:nc,function(x) (1:nc)+(x-1)*2*nc))
    pn <- pn[!(pn%in%(1:(nc)+(1:(nc)-1)*2*nc))]
    pn <- pn[!(pn%in%which(AM==1))]
    
    AM[sample(pn,nflux-sum(AM))] <- 1
    for(j in 1:nc){ AM[j,(AM[j,1:nc]+AM[1:nc,j])==2] <- 0;}
  }
  return(AM)
}

ui <- fluidPage(
    tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 1000px; /* or 950px */}
                                            body { color: #222; }")),
    br(),
    hr(),
    titlePanel(NULL,windowTitle = "BSM2: Static Network Models"),
    h3("Biological Systems Modeling - Module 2: Static Network Models"),
    br(),
    tabsetPanel(
      tabPanel("Overview",
               br(),
               p("This is the second of the Biomedical Systems and Modeling modules that are intended to provide simple interfaces
                  to example systems models and to provide extra practice problems, and is the follow up to ",code("Module 1: Intro to Modeling")),
               p(code("Module 2: Static Network Models"),"provides interfaces to example static network models in line with",
                  strong("Chapter 3: Static Network Models"),"of",em("A First Course in Systems 
                  Biology."),"The expected outcomes of this module are  (i) learn the basics of working with graphs; (ii) 
                  practice computing statistics based on graphs; (iii) practice using stoichiometric network models."),
               tags$hr(),
               h4("Examples"),
               tags$ol(
                 tags$li("Undirected Networks"),
                 tags$ul(
                   tags$li("Simple interface to generate/manipulate a 10 node undirected network.")
                 ),
                 tags$li("Stoichiometric Networks"),
                 tags$ul(
                   tags$li("A stoichiometric network model containing 5 species and 8 chemical reactions. Toggles allow chemical 
                           reactions to be manipulated to maximize outputs."
                   )
                 )
               ),
               h4("Futher Practice"),
               tags$ol(
                 tags$li("Undirected Graphs"),
                 tags$ul(
                   tags$li("Every time this app is loaded a new 6 node undirected graph is generated. For each undirected graph, the 
                           app generates a plot of the graph, the adjacency matrix of the graph, and calculates graph statistics which 
                           includes node degree, node/graph clustering coefficient, shortest path between each pair of nodes, the 
                           diameter of the graph, and the number of paths of length 3 between each pair of nodes. The graph, adjacency 
                           matrix, and statistics only appear after clicking the associated button."
                   )
                 ),
                 tags$li("Directed Graphs"),
                 tags$ul(
                   tags$li("Every time this app is loaded a new 6 node directed graph is generated. For each directed graph, the app 
                            generates a plot of the graph, the adjacency matrix of the graph, and calculates graph statistics which 
                            includes node degree (in and out), node/graph clustering coefficient, shortest path between each pair of 
                            nodes, the diameter of the graph, and the number of paths of length 3 between each pair of nodes. The graph, 
                            adjacency matrix, and statistics only appear after clicking the associated button." 
                   )
                 ),
                 tags$li("Stoichiometric Matrices"),
                 tags$ul(
                   tags$li("Every time this app is loaded a new stoichiometric network is generated. Each network contains 5 species, 
                            with 1-3 carbons each, and approximately 12 fluxes. For each network, the app generates a plot of the directed 
                            graph representation of the network with the node color indicating the number of carbons in each species 
                            (blue = 1, yellow = 2, red = 3). Based on the directed graph and number of carbons per species the app also 
                            generates the corresponding stoichiometric matrix. The network and stoichiometric matrix only appear after 
                            clicking the associated button."
                   )
                 )
               ),
               tags$hr(),
               h4("Citation"),
               p("Voit, E.O.:",em("A First Course in Systems Biology."),"Garland Science, New York, NY, 2017, 2nd edition"),
               br()
      ),              
      tabPanel("Undirected Networks",
        h4("Undirected Networks"),
        tags$ul(
          tags$li("Simple interface to generate/manipulate a 10 node undirected network."),
           tags$li("This app initially generates a random 10 node undirected graph, and 
                   plots the corresponding degree distribution on a log-log plot. The edges 
                   of the graph can be edited by exchanging 0s and 1s in the adjacency matrix."
           ),
          tags$li("Try to make the graph into a small-world network by editing the edges in the network. 
                   Judge the result visually by the regression line in the degree graph."
                  )
        ),
        br(),
        sidebarLayout(
        sidebarPanel(
              h4("Adjacency Matrix"),
              tags$style(".shiny-input-container {margin-bottom: 0px} .checkbox { margin-top: 0px; margin-bottom: 0px;}"),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",1),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",2),
                div(style="display: inline-block;",checkboxInput("b1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",3),
                div(style="display: inline-block;",checkboxInput("c1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("c2", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",4),
                div(style="display: inline-block;",checkboxInput("d1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("d2", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("d3", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",5),
                div(style="display: inline-block;",checkboxInput("e1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("e2", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("e3", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("e4", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",6),
                div(style="display: inline-block;",checkboxInput("f1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("f2", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("f3", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("f4", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("f5", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",7),
                div(style="display: inline-block;",checkboxInput("g1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("g2", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("g3", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("g4", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("g5", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("g6", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",8),
                div(style="display: inline-block;",checkboxInput("h1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("h2", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("h3", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("h4", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("h5", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("h6", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("h7", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",9),
                div(style="display: inline-block;",checkboxInput("i1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("i2", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("i3", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("i4", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("i5", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("i6", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("i7", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("i8", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px; vertical-align:top;",10),
                div(style="display: inline-block;",checkboxInput("j1", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("j2", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("j3", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("j4", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("j5", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("j6", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("j7", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("j8", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block;",checkboxInput("j9", label = "", value = runif(1)>0.5,width=20)),
                div(style="display: inline-block; vertical-align:top;","X")
              ),
              div(
                div(style="display: inline-block; width: 20px;"),
                div(style="display: inline-block; width: 20px;",1),
                div(style="display: inline-block; width: 20px;",2),
                div(style="display: inline-block; width: 20px;",3),
                div(style="display: inline-block; width: 20px;",4),
                div(style="display: inline-block; width: 20px;",5),
                div(style="display: inline-block; width: 20px;",6),
                div(style="display: inline-block; width: 20px;",7),
                div(style="display: inline-block; width: 20px;",8),
                div(style="display: inline-block; width: 20px;",9),
                div(style="display: inline-block; width: 20px;",10)
              )
        ),
        mainPanel(
          column(6,
                 plotOutput("plotSys")
          ),
          column(6,
                 plotOutput("plotSys2")
          ),
          column(8,
                HTML("<strong>Degree of nodes</strong>"),
                htmlOutput("nodeDeg"),
                HTML("<strong>Clustering coefficient</strong>"),
                htmlOutput("clustCoef")
          ),
          column(6,
                 HTML("<strong>Graph clustering coefficient</strong>"),
                 htmlOutput("clustCoefG")
          ),
          column(6,
                HTML("<strong>Diameter</strong>"),
                htmlOutput("diameter")
          ),
          br()
        )
      )
    ),
    tabPanel("Stoichiometric Networks",
             h4("Stoichiometric Networks"),
             tags$ul(
               tags$li("A stoichiometric network model containing 5 species and 8 chemical reactions. Toggles allow chemical 
                           reactions to be manipulated to maximize outputs."),
               tags$li("Change the reaction reaction fluxes to maximize the magenta flux.")
             ),
             br(),
             sidebarLayout(
               sidebarPanel(h4("Chemical Reaction Fluxes"),
                            sliderInput("flux1", "F1", min = 0, max = 10,value = round(runif(1)*100)/10),
                            sliderInput("flux2", "F2", min = 0, max = 10,value = round(runif(1)*100)/10),
                            sliderInput("flux3", "F3", min = 0, max = 10,value = round(runif(1)*100)/10),
                            sliderInput("flux4", "F4", min = 0, max = 10,value = round(runif(1)*100)/10),
                            sliderInput("flux5", "F5", min = 0, max = 10,value = round(runif(1)*100)/10),
                            sliderInput("flux6", "F7", min = 0, max = 10,value = round(runif(1)*100)/10),
                            sliderInput("flux7", "F9", min = 0, max = 10,value = round(runif(1)*100)/10),
                            sliderInput("flux8", "F11", min = 0, max = 10,value = round(runif(1)*100)/10)
                  
               ),
               mainPanel(
                 column(6, align="center",
                        plotOutput("plotSsys"),
                        fluidRow(htmlOutput("stoMat"))
                        ),
                 column(6, align="center",
                        plotOutput("plotFlux"),
                        tableOutput("stoichTab")
                        )
               )
             )
    ),
    tabPanel("Undirected Graphs",
             h4("Undirected Graphs"),
             tags$ul(
               tags$li("Every time this app is loaded a new 6 node undirected graph is generated. For each undirected graph, the 
                     app generates a plot of the graph, the adjacency matrix of the graph, and calculates graph statistics which 
                       includes node degree, node/graph clustering coefficient, shortest path between each pair of nodes, the 
                       diameter of the graph, and the number of paths of length 3 between each pair of nodes. The graph, adjacency 
                       matrix, and statistics only appear after clicking the associated button."),
               tags$li("Reveal the graph. Determine the adjacency matrix and check your answer on completion. Based on the 
                     graph and adjacency network calculate the statistics. Check your answers."),
               tags$li("Reveal the adjacency matrix. Draw the graph and check your answer on completion."),
               tags$li("Reveal the statistics. Determine the adjacency matrix and draw the graph. Check your answers.")
             ),
             br(),
             fluidRow(
               column(4,actionButton("showGraphUD", "Show Graph"),
                      plotOutput("plotSysUD")
               ),
               column(1,""),
               column(4,actionButton("showMatUD", "Show Adjacency Matrix"),
                      htmlOutput("adjMatUD"))
             ),
             actionButton("showStatsUD", "Show Graph Statistics"),
             fluidRow(column(6,HTML("<h3>Degree of nodes</h3>"),
                             htmlOutput("nodeDegUD"),
                             HTML("<h3>Clustering coefficient</h3>"),
                             htmlOutput("clustCoefUD"),
                             HTML("<h3>Clustering coefficient of graph</h3>"),
                             htmlOutput("clustCoefGUD")),
                      column(3,HTML("<h3>Shortest paths</h3>"),
                             htmlOutput("shortPathUD"),
                             HTML("<h3>Diameter</h3>"),
                             htmlOutput("diameterUD")),
                      column(3,HTML("<h3>Paths of length 3</h3>"),
                             htmlOutput("pathL3UD"))
             ),
             fluidRow(HTML("</br></br>"))
    ),
    tabPanel("Directed Graphs",
           h4("Directed Graphs"),
           tags$ul(
             tags$li("Every time this app is loaded a new 6 node directed graph is generated. For each directed graph, the app 
                     generates a plot of the graph, the adjacency matrix of the graph, and calculates graph statistics which 
                     includes node degree (in and out), node/graph clustering coefficient, shortest path between each pair of 
                     nodes, the diameter of the graph, and the number of paths of length 3 between each pair of nodes. The graph, 
                     adjacency matrix, and statistics only appear after clicking the associated button."),
             tags$li("Reveal the graph. Determine the adjacency matrix and check your answer on completion. Based on the 
                     graph and adjacency network calculate the statistics. Check your answers."),
             tags$li("Reveal the adjacency matrix. Draw the graph and check your answer on completion."),
             tags$li("Reveal the statistics. Determine the adjacency matrix and draw the graph. Check your answers.")
           ),
           br(),
           fluidRow(
             column(4,actionButton("showGraphD", "Show Graph"),
                    plotOutput("plotSysD")
             ),
             column(1,""),
             column(4,actionButton("showMatD", "Show Adjacency Matrix"),
                    htmlOutput("adjMatD"))
           ),
           actionButton("showStatsD", "Show Graph Statistics"),
           fluidRow(column(6,HTML("<h3>In-degree</h3>"),
                           htmlOutput("nodeDegIn"),
                           HTML("<h3>Out-degree</h3>"),
                           htmlOutput("nodeDegOut"),
                           HTML("<h3>Clustering coefficient</h3>"),
                           htmlOutput("clustCoefD"),
                           HTML("<h3>Clustering coefficient of graph</h3>"),
                           htmlOutput("clustCoefGD")),
                    column(3,HTML("<h3>Shortest paths</h3>"),
                           htmlOutput("shortPathD"),
                           HTML("<h3>Diameter</h3>"),
                           htmlOutput("diameterD")),
                    column(3,HTML("<h3>Paths of length 3</h3>"),
                           htmlOutput("pathL3D"))
           ),
           fluidRow(HTML("</br></br>"))
  
    ),
    tabPanel("Stoichiometric Matrices",
             h4("Stoichiometric Matrices"),
             tags$ul(
               tags$li("Every time this app is loaded a new stoichiometric network is generated. Each network contains 5 species, 
                        with 1-3 carbons each, and approximately 12 fluxes. For each network, the app generates a plot of the directed 
                        graph representation of the network with the node color indicating the number of carbons in each species 
                        (blue = 1, yellow = 2, red = 3). Based on the directed graph and number of carbons per species the app also 
                        generates the corresponding stoichiometric matrix. The network and stoichiometric matrix only appear after 
                        clicking the associated button."),
               tags$li("Reveal the network. Determine the stoichiometric matrix. Check your answer on completion. How does the
                       stoichiometric matrix compare to the adjacency matrix that would represent the same network? "),
               tags$li("Based on the stoichiometric matrices, setup a system of equations to determine the unknown flux values using 
                       the given measured values."),
               tags$li("Reveal the stoichiometric matrix. Draw the directed graph and determine the number of carbons for 
                       each species. Check your answer on completion.")
             ),
             br(),
             fluidRow(
               column(4,actionButton("showNetP", "Show Network"),
                      plotOutput("plotSSysP")),
               column(3,HTML("<br/><br/>"),tableOutput("stoichTabP")),
               column(5,actionButton("showMatP", "Show Matrix"),
                      htmlOutput("stoMatP"),
                      br(),br(),br(),br(),br(),
                      actionButton("showVals", "Show Unknown Fluxes"),
                      tableOutput("ftab"))
             ),
             
             fluidRow(HTML("</br></br>"))
    )
    )
)

server <- function(input, output, session) {
  
  sysOut <- reactive({
    runDegDist(input)
  })
  
  output$plotSys <- renderPlot({
    
    net.bg <- graph_from_adjacency_matrix(sysOut())
    V(net.bg)$size <- 18
    V(net.bg)$frame.color <- "black"
    V(net.bg)$color <- "orange"
    #V(net.bg)$label <- "" 
    E(net.bg)$arrow.mode <- 0
    #l <- layout_in_circle(net.bg)
    plot(net.bg)
  })
  
  output$plotSys2 <- renderPlot({
    adjMat <- sysOut()
    deg <- diag(adjMat%*%adjMat)
    logD <- log(sort(unique(deg)))
    logN <- log(sapply(sort(unique(deg)),function(x) sum(deg==x))/nrow(adjMat))
    plot(exp(logD),exp(logN),log="xy",pch = 16, cex = 2.5,cex.lab=1.5,cex.axis=1.5,col=1,xlab="Degree of node",ylab="Fraction of nodes")
    points(exp(logD),exp(logN),pch = 16, cex = 2,col=3)
    model <- lm(logN~logD)
    lines(exp(logD), exp(predict(model, newdata=list(logD=logD))))
    
  })
  
  output$nodeDeg <- renderUI({
    am <- sysOut()
    d <- apply(am,1,sum)
    
    HTML(paste("<p>",paste(d,collapse = "&nbsp;&nbsp;"),"</p>"))
  })
  
  output$clustCoef <- renderUI({ 
    am <- sysOut()
    d <- apply(am,1,sum)
    cc <- NULL
    for(i in 1:10){
      neigh <- which(am[i,]==1|am[,i]==1)
      nMat <- am[neigh,neigh]
      
      cc <- c(cc,ifelse(d[i]>1,(sum(nMat)/2)/(d[i]*(d[i]-1)/2),0))
    }
    
    HTML(paste("<p>",paste(sprintf("%.2f",cc),collapse = "&nbsp;&nbsp;"),"</p>"))
  })
  
  output$clustCoefG <- renderUI({ 
    am <- sysOut()
    d <- apply(am,1,sum)
    cc <- NULL
    for(i in 1:nrow(am)){
      neigh <- which(am[i,]==1|am[,i]==1)
      nMat <- am[neigh,neigh]
      
      cc <- c(cc,ifelse(d[i]>1,(sum(nMat)/2)/(d[i]*(d[i]-1)/2),0))
    }
    
    HTML(paste("<p>",sprintf("%.2f",mean(cc)),"</p>"))
  })
  
  output$diameter <- renderUI({
    am <- sysOut()
    sd <- matrix(0,nrow = 10,ncol = 10)
    ad <- am
    
    count <- 1
    while(sum(sd==0)>0 & count < nrow(am)){
      sd[ad!=0 & sd==0] <- count
      ad <- ad%*%am
      count <- count + 1
    }

    sd[(lower.tri(sd)|upper.tri(sd)) & sd==0] <- Inf
    diag(sd) <- 0
    
    HTML(paste("<p>",max(sd),"</p>"))
  })
  
  inf <- 1
  outf <- 2:5
  nc <- 5
  nflux <- 12
  
  am <- array(0,dim=c(10,10))
  am[1,1:5] <- c(0,1,1,1,1)
  am[2,1:5] <- c(0,0,1,0,0)
  am[3,1:5] <- c(0,0,0,1,0)
  am[4,1:5] <- c(0,0,0,0,1)
  am[5,1:5] <- c(0,1,0,0,0)
    
  nn2 <- 1
  for(i in inf){
    am[nc+nn2,i] <- 1
    nn2 <- nn2 + 1
  }
  
  for(o in outf){
    am[o,nc+nn2] <- 1
    nn2 <- nn2 + 1
  }
  
  sMat <- genStoichMat(am,nc)
  selF <- sample(c(6,8,10,12),1)
  
  output$plotSsys <- renderPlot({
    
    net.bg <- graph_from_adjacency_matrix(am,mode="directed")
    l <- layout_with_kk(net.bg)
    
    colors <- c("royalblue1", "gold","tomato", "white")
    V(net.bg)$size <- 16
    V(net.bg)$frame.color <- c(rep("black",nc),rep("white",5))
    clrs <- colors[c(sMat$ncc,rep(4,5))]
    V(net.bg)$color <- clrs
    E(net.bg)$edge.label <- as.character(1:sum(am))
    E(net.bg)$color <- rep("gray",13)
    E(net.bg)$color[selF] <- "magenta"
    
    ew <- sysStoichNet()
    ew[ew<=0] <- 0
    ew <- ew + 0.5
    ew <- 5*ew/max(ew)

    plot(net.bg,edge.arrow.size=.4,layout=l,
         vertex.label.color= c(rep("black",nc),rep("white",5)),
         edge.label=as.character(1:sum(am)),
         edge.label.color="black",edge.label.font=2,edge.width = ew,edge.arrow.size	= 5,edge.arrow.width	= ew/2+1)
    par(mar=rep(0.5,4))
  })
  
  output$stoMat <- renderUI({
    
    HTML(paste("<strong>Stoichiometric Matrix</strong>",
               paste("<p>","&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",
                     "<font color=\"#D35400\">",paste(sapply(1:13, function(x) paste0(x,ifelse(x<9,"&nbsp;&nbsp;","&nbsp;"))),collapse=""),
                     "</font>","<br/>",paste(sMat$html,collapse="<br/>"),"</p>"))
    )
  })
  
  output$plotFlux <- renderPlot({
    cols <- rep("gray",13)
    cols[selF] <- "magenta"
    fVals <- sysStoichNet()
    names(fVals) <- 1:13
    barplot(fVals,ylab = "Flux values",xlab="Fluxes",col=cols)
  })
  
  output$stoichTab <- renderTable({
    data.frame(Species=1:nc,N_carbons = sMat$ncc)
  })

  sysStoichNet <- reactive({
    runStoichSys(input,selF,sMat$sm)
  })
  
  amUD <- matrix(0,nrow = 6,ncol = 6)
  amUD[upper.tri(amUD)] <- (runif((6^2 - 6)/2) < (0.1*runif(1)+0.5)) + 0
  amUD[lower.tri(amUD)] <- t(amUD)[lower.tri(amUD)]
  
  output$adjMatUD <- renderUI({
    if(!input$showMatUD) return()     
    
    HTML(paste("<h3>Adjacency Matrix</h3>",
               paste("<h3>",
                     paste(paste("<font color=\"red\">0</font>",amUD[1,2],amUD[1,3],amUD[1,4],amUD[1,5],amUD[1,6],sep="&nbsp;&nbsp;"),
                           paste(amUD[2,1],"<font color=\"red\">0</font>",amUD[2,3],amUD[2,4],amUD[2,5],amUD[2,6],sep="&nbsp;&nbsp;"),
                           paste(amUD[3,1],amUD[3,2],"<font color=\"red\">0</font>",amUD[3,4],amUD[3,5],amUD[3,6],sep="&nbsp;&nbsp;"),
                           paste(amUD[4,1],amUD[4,2],amUD[4,3],"<font color=\"red\">0</font>",amUD[4,5],amUD[4,6],sep="&nbsp;&nbsp;"),
                           paste(amUD[5,1],amUD[5,2],amUD[5,3],amUD[5,4],"<font color=\"red\">0</font>",amUD[5,6],sep="&nbsp;&nbsp;"),
                           paste(amUD[6,1],amUD[6,2],amUD[6,3],amUD[6,4],amUD[6,5],"<font color=\"red\">0</font>",sep="&nbsp;&nbsp;")
                           ,sep="<br/>"),"</h3>"))
    )

  })
  
  
  output$nodeLabels <- renderUI({ 
    HTML(paste("<br/><br/><br/><h3>",paste(1:10,collapse="<br/><br/>"),"</h3>"))    
  })
  
  output$nodeDegUD <- renderUI({
    if(!input$showStatsUD) return()
    
    d <- apply(amUD,1,sum)
    
    HTML(paste("<h4>",paste(d,collapse = "&nbsp;&nbsp;"),"</h4>"))
  })
  
  output$clustCoefUD <- renderUI({ 
    if(!input$showStatsUD) return()
    
    d <- apply(amUD,1,sum)
    cc <- NULL
    for(i in 1:nrow(amUD)){
      neigh <- which(amUD[i,]==1|amUD[,i]==1)
      nMat <- amUD[neigh,neigh]
      
      cc <- c(cc,ifelse(d[i]>1,(sum(nMat)/2)/(d[i]*(d[i]-1)/2),0))
    }
    
    HTML(paste("<h4>",paste(sprintf("%.2f",cc),collapse = "&nbsp;&nbsp;"),"</h4>"))
  })
  
  output$clustCoefGUD <- renderUI({ 
    if(!input$showStatsUD) return()
    
    d <- apply(amUD,1,sum)
    cc <- NULL
    for(i in 1:nrow(amUD)){
      neigh <- which(amUD[i,]==1|amUD[,i]==1)
      nMat <- amUD[neigh,neigh]
      
      cc <- c(cc,ifelse(d[i]>1,(sum(nMat)/2)/(d[i]*(d[i]-1)/2),0))
    }
    
    HTML(paste("<h4>",sprintf("%.2f",mean(cc)),"</h4>"))
  })
  
  output$shortPathUD <- renderUI({ 
    if(!input$showStatsUD) return()
    
    sd <- matrix(0,nrow = nrow(amUD),ncol = nrow(amUD))
    ad <- amUD
    count <- 1
    while(sum(sd==0)>0 & count < nrow(amUD)){
      sd[ad!=0 & sd==0] <- count
      
      ad <- ad%*%amUD
      count <- count + 1
    }
    sd[(lower.tri(sd)|upper.tri(sd)) & sd==0] <- Inf
    
    sdHTML <- paste(paste("<font color=\"red\">0</font>",sd[1,2],sd[1,3],sd[1,4],sd[1,5],sd[1,6],sep="&nbsp;&nbsp;"),
                    paste(sd[2,1],"<font color=\"red\">0</font>",sd[2,3],sd[2,4],sd[2,5],sd[2,6],sep="&nbsp;&nbsp;"),
                    paste(sd[3,1],sd[3,2],"<font color=\"red\">0</font>",sd[3,4],sd[3,5],sd[3,6],sep="&nbsp;&nbsp;"),
                    paste(sd[4,1],sd[4,2],sd[4,3],"<font color=\"red\">0</font>",sd[4,5],sd[4,6],sep="&nbsp;&nbsp;"),
                    paste(sd[5,1],sd[5,2],sd[5,3],sd[5,4],"<font color=\"red\">0</font>",sd[5,6],sep="&nbsp;&nbsp;"),
                    paste(sd[6,1],sd[6,2],sd[6,3],sd[6,4],sd[6,5],"<font color=\"red\">0</font>",sep="&nbsp;&nbsp;")
                    ,sep="<br/>")
    
    sdHTML <- gsub("Inf","&nbsp;I",sdHTML)
    HTML(paste("<h4>",sdHTML,"</h4>"))
  })
  
  output$pathL3UD <- renderUI({ 
    if(!input$showStatsUD) return()
    
    sd <- amUD
    for(i in 1:2){
      sd <- sd%*%amUD
    }
    
    sdHTML <- paste(paste(sd[1,1],sd[1,2],sd[1,3],sd[1,4],sd[1,5],sd[1,6],sep="&nbsp;&nbsp;"),
                    paste(sd[2,1],sd[2,2],sd[2,3],sd[2,4],sd[2,5],sd[2,6],sep="&nbsp;&nbsp;"),
                    paste(sd[3,1],sd[3,2],sd[3,3],sd[3,4],sd[3,5],sd[3,6],sep="&nbsp;&nbsp;"),
                    paste(sd[4,1],sd[4,2],sd[4,3],sd[4,4],sd[4,5],sd[4,6],sep="&nbsp;&nbsp;"),
                    paste(sd[5,1],sd[5,2],sd[5,3],sd[5,4],sd[5,5],sd[5,6],sep="&nbsp;&nbsp;"),
                    paste(sd[6,1],sd[6,2],sd[6,3],sd[6,4],sd[6,5],sd[6,6],sep="&nbsp;&nbsp;")
                    ,sep="<br/>")
    
    HTML(paste("<h4>",sdHTML,"</h4>"))
  })
  
  output$diameterUD <- renderUI({ 
    if(!input$showStatsUD) return()
    
    sd <- matrix(0,nrow = nrow(amUD),ncol = nrow(amUD))
    ad <- amUD
    count <- 1
    while(sum(sd==0)>0 & count < nrow(amUD)){
      sd[ad!=0 & sd==0] <- count
      
      ad <- ad%*%amUD
      count <- count + 1
    }
    sd[(lower.tri(sd)|upper.tri(sd)) & sd==0] <- Inf
    
    HTML(paste("<h4>",max(sd),"</h4>"))
  })
  
  observeEvent(input$showGraphUD, {
    #show("plotSys")
  })
  
  observeEvent(input$showMatUD, {
    #show("adjMat")
  })
  
  observeEvent(input$showStatsUD, {
    #show("nodeDegreeUD")
    #show("clustCoefUD")
    #show("clustCoefGUD")
    #show("shortPathUD")
    #show("diameterUD")
  })
  
  output$plotSysUD <- renderPlot({
    if(!input$showGraphUD) return() 
    
    net.bg <- graph_from_adjacency_matrix(amUD)
    V(net.bg)$size <- 18
    V(net.bg)$frame.color <- "black"
    V(net.bg)$color <- "orange"
    
    E(net.bg)$arrow.mode <- 0
    
    plot(net.bg)
  })
  
  amD <- matrix(0,nrow = 6,ncol = 6)
  amD[upper.tri(amD)] <- (runif((6^2 - 6)/2) < (0.1*runif(1)+0.5)) + 0
  amD[lower.tri(amD)] <- (runif((6^2 - 6)/2) < (1-0.2*t(amD)[lower.tri(amD)])*(0.1*runif(1)+0.5)) + 0
  
  output$adjMatD <- renderUI({
    if(!input$showMatD) return()
    
    HTML(paste("<h3>Adjacency Matrix</h3>",
               paste("<h3>",
                     paste(paste("<font color=\"red\">0</font>",amD[1,2],amD[1,3],amD[1,4],amD[1,5],amD[1,6],sep="&nbsp;&nbsp;"),
                           paste(amD[2,1],"<font color=\"red\">0</font>",amD[2,3],amD[2,4],amD[2,5],amD[2,6],sep="&nbsp;&nbsp;"),
                           paste(amD[3,1],amD[3,2],"<font color=\"red\">0</font>",amD[3,4],amD[3,5],amD[3,6],sep="&nbsp;&nbsp;"),
                           paste(amD[4,1],amD[4,2],amD[4,3],"<font color=\"red\">0</font>",amD[4,5],amD[4,6],sep="&nbsp;&nbsp;"),
                           paste(amD[5,1],amD[5,2],amD[5,3],amD[5,4],"<font color=\"red\">0</font>",amD[5,6],sep="&nbsp;&nbsp;"),
                           paste(amD[6,1],amD[6,2],amD[6,3],amD[6,4],amD[6,5],"<font color=\"red\">0</font>",sep="&nbsp;&nbsp;")
                           ,sep="<br/>"),"</h3>"))
    )
  })
  
  output$nodeDegOut <- renderUI({
    if(!input$showStatsD) return()
    
    d <- apply(amD,1,sum)
    
    HTML(paste("<h4>",paste(d,collapse = "&nbsp;&nbsp;"),"</h4>"))
  })
  
  output$nodeDegIn <- renderUI({
    if(!input$showStatsD) return()
    
    d <- apply(amD,2,sum)
    
    HTML(paste("<h4>",paste(d,collapse = "&nbsp;&nbsp;"),"</h4>"))
  })
  
  output$clustCoefD <- renderUI({ 
    if(!input$showStatsD) return()
    
    nNeigh <- sapply(1:nrow(amD), function(x) sum(amD[x,]==1|amD[,x]==1))
    cc <- NULL
    for(i in 1:nrow(amD)){
      neigh <- which(amD[i,]==1|amD[,i]==1)
      nMat <- amD[neigh,neigh]
      
      cc <- c(cc,ifelse(nNeigh[i]>1,(sum(nMat))/(nNeigh[i]*(nNeigh[i]-1)),0))
    }
    
    HTML(paste("<h4>",paste(sprintf("%.2f",cc),collapse = "&nbsp;&nbsp;"),"</h4>"))
  })
  
  output$clustCoefGD <- renderUI({ 
    if(!input$showStatsD) return()
    
    nNeigh <- sapply(1:nrow(amD), function(x) sum(amD[x,]==1|amD[,x]==1))
    cc <- NULL
    for(i in 1:nrow(amD)){
      neigh <- which(amD[i,]==1|amD[,i]==1)
      nMat <- am[neigh,neigh]
      
      cc <- c(cc,ifelse(nNeigh[i]>1,(sum(nMat))/(nNeigh[i]*(nNeigh[i]-1)),0))
    }
    
    HTML(paste("<h4>",sprintf("%.2f",mean(cc)),"</h4>"))
  })
  
  output$shortPathD <- renderUI({ 
    if(!input$showStatsD) return()
    
    sd <- matrix(0,nrow = nrow(amD),ncol = nrow(amD))
    ad <- t(amD)
    count <- 1
    while(sum(sd==0)>0 & count < nrow(amD)){
      sd[ad!=0 & sd==0] <- count
      
      ad <- ad%*%t(amD)
      count <- count + 1
    }
    sd <- t(sd)
    sd[(lower.tri(sd)|upper.tri(sd)) & sd==0] <- Inf
    
    sdHTML <- paste(paste("<font color=\"red\">0</font>",sd[1,2],sd[1,3],sd[1,4],sd[1,5],sd[1,6],sep="&nbsp;&nbsp;"),
                    paste(sd[2,1],"<font color=\"red\">0</font>",sd[2,3],sd[2,4],sd[2,5],sd[2,6],sep="&nbsp;&nbsp;"),
                    paste(sd[3,1],sd[3,2],"<font color=\"red\">0</font>",sd[3,4],sd[3,5],sd[3,6],sep="&nbsp;&nbsp;"),
                    paste(sd[4,1],sd[4,2],sd[4,3],"<font color=\"red\">0</font>",sd[4,5],sd[4,6],sep="&nbsp;&nbsp;"),
                    paste(sd[5,1],sd[5,2],sd[5,3],sd[5,4],"<font color=\"red\">0</font>",sd[5,6],sep="&nbsp;&nbsp;"),
                    paste(sd[6,1],sd[6,2],sd[6,3],sd[6,4],sd[6,5],"<font color=\"red\">0</font>",sep="&nbsp;&nbsp;")
                    ,sep="<br/>")
    
    sdHTML <- gsub("Inf","&nbsp;I",sdHTML)
    HTML(paste("<h4>",sdHTML,"</h4>"))
  })
  
  output$pathL3D <- renderUI({ 
    if(!input$showStatsD) return()
    
    sd <- amD
    for(i in 1:2){
      sd <- sd%*%amD
    }
    
    sdHTML <- paste(paste(sd[1,1],sd[1,2],sd[1,3],sd[1,4],sd[1,5],sd[1,6],sep="&nbsp;&nbsp;"),
                    paste(sd[2,1],sd[2,2],sd[2,3],sd[2,4],sd[2,5],sd[2,6],sep="&nbsp;&nbsp;"),
                    paste(sd[3,1],sd[3,2],sd[3,3],sd[3,4],sd[3,5],sd[3,6],sep="&nbsp;&nbsp;"),
                    paste(sd[4,1],sd[4,2],sd[4,3],sd[4,4],sd[4,5],sd[4,6],sep="&nbsp;&nbsp;"),
                    paste(sd[5,1],sd[5,2],sd[5,3],sd[5,4],sd[5,5],sd[5,6],sep="&nbsp;&nbsp;"),
                    paste(sd[6,1],sd[6,2],sd[6,3],sd[6,4],sd[6,5],sd[6,6],sep="&nbsp;&nbsp;")
                    ,sep="<br/>")
    
    HTML(paste("<h4>",sdHTML,"</h4>"))
  })
  
  output$diameterD <- renderUI({ 
    if(!input$showStatsD) return()
    
    sd <- matrix(0,nrow = nrow(amD),ncol = nrow(amD))
    ad <- t(amD)
    count <- 1
    while(sum(sd==0)>0 & count < nrow(amD)){
      sd[ad!=0 & sd==0] <- count
      
      ad <- ad%*%t(amD)
      count <- count + 1
    }
    sd[(lower.tri(sd)|upper.tri(sd)) & sd==0] <- Inf
    
    HTML(paste("<h4>",max(sd),"</h4>"))
  })
  
  observeEvent(input$showGraphD, {
    #show("plotSys")
  })
  
  observeEvent(input$showMatD, {
    #show("adjMat")
  })
  
  observeEvent(input$showStatsD, {
    #show("nodeDegree")
    #show("clustCoef")
    #show("clustCoefG")
    #show("shortPath")
    #show("diameter")
  })
  
  output$plotSysD <- renderPlot({
    if(!input$showGraphD) return()      
    net.bg <- graph_from_adjacency_matrix(t(amD))
    V(net.bg)$size <- 18
    V(net.bg)$frame.color <- "black"
    V(net.bg)$color <- "orange"
    E(net.bg)$arrow.mode <- 1
    
    plot(net.bg,edge.arrow.size	= 0.6)
  })
  
  output$stoMatP <- renderUI({
    if(!input$showMatP) return()
    HTML(paste("<h4>Stoichiometric Matrix</h4>",
               paste("<h4>","&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;",
                     "<font color=\"#D35400\">",paste(sapply(1:sum(AM), function(x) paste0(x,ifelse(x<9,"&nbsp;&nbsp;","&nbsp;"))),collapse=""),
                     "</font>","<br/>",paste(sMatP$html,collapse="<br/>"),"</h4>"))
    )
  })
  
  output$plotSSysP <- renderPlot({
    if(!input$showNetP) return()
    
    net.bg <- graph_from_adjacency_matrix(AM,mode="directed")
    l <- layout_with_kk(net.bg)
    
    colors <- c("royalblue1", "gold","tomato", "white")
    V(net.bg)$size <- 16
    V(net.bg)$frame.color <- c(rep("black",nc),rep("white",5))
    clrs <- colors[c(sMatP$ncc,rep(4,5))]
    V(net.bg)$color <- clrs
    E(net.bg)$edge.label <- as.character(1:sum(AM))
    
    ew <- 3
    plot(net.bg,edge.arrow.size=.4,layout=l,
         vertex.label.color= c(rep("black",nc),rep("white",5)),
         edge.label=as.character(1:sum(AM)),
         edge.label.color="black",edge.label.font=2,edge.width = ew,edge.arrow.size	= ew*3,edge.arrow.width	= ew/2+1)
    par(mar=rep(0.5,4))
  })
  
  observeEvent(input$showMatP, {
    #show("stoMat")
  })
  
  observeEvent(input$showNetP, {
    #show("plotSys")
  })
  
  output$stoichTabP <- renderTable({
    if(!input$showNetP) return()
    data.frame(Species=1:nc,N_carbons = sMatP$ncc)
  })
  
  AM <- getStoichNet()
  sMatP <- genStoichMat(AM,nc)
  fTabP <- runStoichSysP(sMatP$sm)
  while(sum(fTabP$Value<0)>0){
    AM <- getStoichNet()
    sMatP <- genStoichMat(AM,nc)
    fTabP <- runStoichSysP(sMatP$sm)
  }
  
  
  output$ftab <- renderTable({
    ftab <- fTabP
    if(!input$showVals){
      ftab$Value[ftab$Type=="Unknown"] <- ""
    }
    return(ftab)
  })
  
}
shinyApp(ui = ui, server = server)