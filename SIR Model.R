#Niall Martin - 4BCT - 12301341 - Assignment 4

#include libraries
library(ggplot2)
library(deSolve)
library(scales)

#Outer sensitivity function
mySystem <- function(numRun)
{
  #Uniform distribution
  minInfected <- 1.0;
  maxInfected <- 20.0;
  minRecoveryDelay <- 2.0; 
  maxRecoveryDelay <- 10.0;
  minCe <- 0.75; 
  maxCe <- 8.0;
  
  #Random inital values for variables within uniform distribution thresholds and beta and gamma values for use in SIR model
  population <- 10000
  infected <- runif(1, minInfected, maxInfected)
  recovered <- 0
  recoveryDelay <- runif(1, minRecoveryDelay, maxRecoveryDelay)
  gammaValue <- 1/recoveryDelay
  Ce <- runif(1, minCe, maxCe)
  betaValue <- Ce/population

  #Initial Stock values
  stock <- c(S = population - infected, I = infected, R = 0)
  
  #parameter values
  parameters <- c(beta = betaValue, gamma = gammaValue)
  
  #Time values
  time <- seq(0, 60, by=0.1)
  
  #SIR model function
  SIR <- function(time, state, parameters) 
  {
    with(as.list(c(state, parameters)), 
    {
      #Formulae for the Susceptible, Infected and Recovered at a given time for SIR model
      #Susceptible
      dS <- -beta * S * I
      #Infected
      dI <- beta * S * I - gamma * I
      #Recovered
      dR <- gamma * I
      
      return(list(c(dS, dI, dR)))
    })
  }
  
  #solve with ode(Ordinary Differential Equation solver)
  result <- as.data.frame(ode(y = stock, times = time, func = SIR, parms = parameters, method="euler"))
  
  #get infected results, take out S and R, not needed
  infectedResult <- result
  infectedResult$S <- NULL
  infectedResult$R <- NULL
  
  #add Ce and RecoveryDelay to infectedResult (for use in correlation coefficient calculations)
  infectedResult$Ce <- Ce
  infectedResult$RD <- recoveryDelay
  
  #Add run number as a label to result vector
  infectedResult$numRun <- numRun
  return(infectedResult)
}

#run initially to create results matrix
results <- mySystem(1)

#loop 50 times and add to results matrix
for(i in 2:50)
{
  results <- rbind(results, mySystem(i))
}

#Get max, min and average peaks of epidemics over the 50 runs
maxPeak <- max(results["I"])
minPeak <- min(results["I"])
averagePeak <- colMeans(results["I"], na.rm = FALSE, dims = 1)

#plot infected over time for the 50 simulations
infectionGraph<-ggplot(results,aes(x=time,y=I,color=numRun,group=numRun)) + geom_line() + ylab("No. Infected") +  xlab("Time")  + guides(color=FALSE) + ggtitle("SIR Model Infections - 50 Simulations")
#print graph
print(infectionGraph)

#Split results into groups of their corresponding time steps
finalResults<-split(results,results$time)

#find Correlation between Infected and Ce
corCE <- sapply(finalResults,function(l){cor(l$I, l$Ce )})

#find correlation between Infected and Recovery Delay
corRD <- sapply(finalResults,function(l){cor(l$I, l$RD )})

#Time vector same as SIR model to be used in graphing correlation coefficients
simTime <- seq(0, 60, by=0.1)

#Graph correlation coefficients for Infected
CCgraph<-ggplot() + geom_line(aes(simTime,corRD,color="Recovery Delay")) + geom_line(aes(simTime,corCE,color="Ce")) + scale_y_continuous(labels = comma) + ylab("Infected Correlation Coefficients") + xlab("Time") + labs(color="") + theme(legend.position="top")
#print graph
print(CCgraph)

#show max min and average infected peak values in colsole
maxPeak
minPeak
averagePeak






