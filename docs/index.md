## R Resources for Learning Comptutational Systems Biology
On this page you will find a variety of resources for teaching/learning Computational Systems Biology that have been developed using the R statistical language. These resources include interactive modules that are generally applicable to learning the basics of systems biology modeling, as well as, coding examples for using R to develop and analyze various types of systems biology models. The coding examples also include R *shiny* examples, which provide a template for developing interactive interfaces for systems biology models. Source codes for all interfaces and coding examples are also available.

### R Shiny Interactive Modules 
This is a series of modules that are intended to provide simple interfaces to example systems models and to provide extra practice problems and are intended as a supplement the textbook: A First Course in Systems Biology . The textbook is not needed to be able to use these applications, but may be useful for further explanation.

1. [Intro to Modeling](https://kieslich.shinyapps.io/sysBio1/): Provides interfaces to selected examples from Chapter 2: Introduction to Mathematical Modeling of A First Course in Systems  Biology. The expected outcomes of this module are (i) learn the parts of a model and (ii) practice predicting how perturbations to model inputs and parameters affect the model response. ([Source](https://github.com/cakieslich/introCSB/tree/main/SysBioModules/sysBio1))
2. [Static Network Models](https://kieslich.shinyapps.io/sysBio2/): Provides interfaces to example static network models in line with Chapter 3: Static Network Models of A First Course in Systems Biology. The expected outcomes of this module are (i) learn the basics of working with graphs; (ii) practice computing statistics based on graphs; (iii) practice using stoichiometric network models. ([Source](https://github.com/cakieslich/introCSB/tree/main/SysBioModules/sysBio2))
3. [Discrete, Recursive Models](https://kieslich.shinyapps.io/sysBio3/): Provides interfaces to example discrete and recursive models, including examples from Chapter 4: The Mathematics of Biological Systems of A First Course in Systems Biology. The expected outcomes of this module are (i) learn the basic principles of discrete, recursive models; (ii) explore examples of dynamic discrete recrusive models; (iii) practice generating Markov matrices. ([Source](https://github.com/cakieslich/introCSB/tree/main/SysBioModules/sysBio3))
4. [Canonical Models](https://kieslich.shinyapps.io/sysBio4/): Provides interfaces to example models based on the canonical models described in pg. 102-105 of Chapter 4: The Mathematics of Biological Systems of A First Course in Systems Biology. The expected outcomes of this module are to learn the basic structure of (i) linear, (ii) Lotka-Voltera, (iii) mass action, and (iv) S-system models. ([Source](https://github.com/cakieslich/introCSB/tree/main/SysBioModules/sysBio4))
5. [Analysis of Dynamic Models](https://kieslich.shinyapps.io/sysBio5/): Provides interfaces to example models described in Chapter 4: The Mathematics of Biological Systems of A First Course in Systems Biology. The expected outcomes of this module are to learn the fundamentals of stability analysis including: (i) linearization of nonlinear models, (ii) eigenvalue analysis. ([Source](https://github.com/cakieslich/introCSB/tree/main/SysBioModules/sysBio5))

### R Coding Examples by Topic
Below is a list of links to HTML pages with coding examples for a range of topics related to developing and analysiing systems biology models. The examples cover developing systems diagrams and developing/analyzing ODE-based systems models.      

1. [Intro to Modeling](https://cakieslich.github.io/introCSB/IntroToModelling.html) 
2. [Properties of Undirected Graphs](https://cakieslich.github.io/introCSB/PropertiesOfUndirGraphs.html) 
3. [Properties of Directed Graphs](https://cakieslich.github.io/introCSB/PropertiesOfDirGraphs.html) 
4. [Discrete, Recursive Models](https://cakieslich.github.io/introCSB/DiscreteModelsSIR.html) 
5. [Markov Models](https://cakieslich.github.io/introCSB/MarkovModels.html) 
6. [Parameter Estimation](https://cakieslich.github.io/introCSB/ParameterEstimation.html) 
7. [Stability of Linear Models](https://cakieslich.github.io/introCSB/StabilityLinearModels.html) 
8. [Linear Approximation of Nonlinear Models](https://cakieslich.github.io/introCSB/LinearApproximation.html) 
9. [Stability of Lotka-Volterra Models](https://cakieslich.github.io/introCSB/StabilityNonLinearModels.html) 
10. [Stability of S-system Models](https://cakieslich.github.io/introCSB/StabilitySsystemModels.html) 
11. [Nullcline Analysis](https://cakieslich.github.io/introCSB/nullclines.html) 
12. [Analyzing Hysteresis](https://cakieslich.github.io/introCSB/hysteresis.html)   

### R *shiny* Examples
*shiny* is an R package that enables the efficient development of web-based applications. The following examples demonstrate how R shiny can be used for developing interactive interfaces for systems biology models.

1. [Interface for ODE based model](https://github.com/cakieslich/introCSB/blob/main/shinyEx/fishing.R): Simple example of how to use R shiny to develop interactive interfaces for ODE-based models. ([Demo](https://kieslich.shinyapps.io/fishing/))
2. [Interface for ODE based model with report export](https://github.com/cakieslich/introCSB/tree/main/shinyEx/mapk): Advanced example of how to use R shiny to develop interactive interfaces for ODE-based models. Includes the ability to export analysis from the model using an R markdown derived report in PDF format. ([Demo](https://kieslich.shinyapps.io/mapk/))


***Citation***: Voit, E.O.: A First Course in Systems Biology. Garland Science, New York, NY, 2017, 2nd edition
