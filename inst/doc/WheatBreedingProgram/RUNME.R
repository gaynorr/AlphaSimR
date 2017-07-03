# This is the first script of the simulation
# All other scripts are called by this script
library(AlphaSimR)

#Path to MaCS
macsPath = "/Users/rgaynor/Documents/testPBS/macs"

# Load global parameters
source("GlobalParameters.R")

# Create initial parents
source("CreateParents.R")

# Fill breeding pipeline with unique individuals from initial parents
source("FillPipeline.R")

# Initialize variables for results
parMean = meanG(Parents)
parVar = varG(Parents)

# Cycle years to make more advanced parents
for(year in 1:20){ #Change to any number of desired years
  w = runif(1) #Sample environmental covariate
  cat("Working on year:",year,"\n")
  source("UpdateParents.R") #Pick parents
  source("AdvanceYear.R") #Advances yield trials by a year
  # Update results
  parMean = c(parMean,meanG(Parents))
  parVar = c(parVar,varG(Parents))
}


#View results
plot(0:20,parMean,type="l",xlab="Year",ylab="Mean GV",
     main="Mean genetic value of parents")
plot(0:20,parVar,type="l",xlab="Year",ylab="Var GV",
     main="Genetic variance of parents")
