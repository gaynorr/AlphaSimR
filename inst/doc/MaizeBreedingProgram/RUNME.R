# This is the first script of the simulation
# All other scripts are called by this script
library(AlphaSimR)

#Path to MaCS
macsPath = "/Users/rgaynor/Documents/testAS/macs"

# Load global parameters
source("GlobalParameters.R")

nReps = 10
burninYears = 20
futureYears = 20
# Initialize variables for results
hybridCorr = inbredMean = hybridMean = inbredVar = hybridVar = 
  rep(NA_real_,burninYears+futureYears)
output = list(inbredMean=NULL,
              inbredVar=NULL,
              hybridMean=NULL,
              hybridVar=NULL,
              hybridCorr=NULL)
saveRDS(output,"Scenario1.rds")
saveRDS(output,"Scenario2.rds")
for(rep in 1:nReps){
  # Create initial parents and set testers and hybrid parents
  source("CreateParents.R")
  # Fill breeding pipeline with unique individuals from initial parents
  source("FillPipeline.R")
  # Cycle years
  for(year in 1:burninYears){ #Change to any number of desired years
    cat("Working on year:",year,"\n")
    source("UpdateParents.R") #Pick new inbred parents
    source("UpdateTesters.R") #Pick new testers and hybrid parents
    source("AdvanceYear.R") #Advances yield trials by a year
    # Update results
    inbredMean[year] = (meanG(MaleParents)+meanG(FemaleParents))/2
    inbredVar[year] = (varG(MaleParents)+varG(FemaleParents))/2
    tmp = hybridCross(FemaleParents,MaleParents)
    hybridMean[year] = meanG(tmp)
    hybridVar[year] = varG(tmp)
    tmp = calcGCA(tmp,use="gv")
    hybridCorr[year] = cor(c(tmp$females[,2],tmp$males[,2]),
                           c(FemaleParents@gv[,1],MaleParents@gv[,1]))
    rm(tmp)
  }
  save.image("tmp.RData")
  
  #Scenario 1
  cat("Working on Scenario 1\n")
  for(year in (burninYears+1):(burninYears+futureYears)){
    cat("Working on year:",year,"\n")
    source("UpdateParents.R") #Pick new inbred parents
    source("UpdateTesters.R") #Pick new testers and hybrid parents
    source("AdvanceYear.R") #Advances yield trials by a year
    # Update results
    inbredMean[year] = (meanG(MaleParents)+meanG(FemaleParents))/2
    inbredVar[year] = (varG(MaleParents)+varG(FemaleParents))/2
    tmp = hybridCross(FemaleParents,MaleParents)
    hybridMean[year] = meanG(tmp)
    hybridVar[year] = varG(tmp)
    tmp = calcGCA(tmp,use="gv")
    hybridCorr[year] = cor(c(tmp$females[,2],tmp$males[,2]),
                           c(FemaleParents@gv[,1],MaleParents@gv[,1]))
    rm(tmp)
  }
  output = readRDS("Scenario1.rds")
  output = list(inbredMean=rbind(output$inbredMean,inbredMean),
                inbredVar=rbind(output$inbredVar,inbredVar),
                hybridMean=rbind(output$hybridMean,hybridMean),
                hybridVar=rbind(output$hybridVar,hybridVar),
                hybridCorr=rbind(output$hybridCorr,hybridCorr))
  saveRDS(output,"Scenario1.rds")
  #Scenario 2
  load("tmp.RData")
  cat("Working on Scenario 2\n")
  for(year in (burninYears+1):(burninYears+futureYears)){
    cat("Working on year:",year,"\n")
    source("UpdateParents.R") #Pick new inbred parents
    source("UpdateTesters.R") #Pick new testers and hybrid parents
    source("AdvanceYearAlt.R") #Advances yield trials by a year
    # Update results
    inbredMean[year] = (meanG(MaleParents)+meanG(FemaleParents))/2
    inbredVar[year] = (varG(MaleParents)+varG(FemaleParents))/2
    tmp = hybridCross(FemaleParents,MaleParents)
    hybridMean[year] = meanG(tmp)
    hybridVar[year] = varG(tmp)
    tmp = calcGCA(tmp,use="gv")
    hybridCorr[year] = cor(c(tmp$females[,2],tmp$males[,2]),
                           c(FemaleParents@gv[,1],MaleParents@gv[,1]))
    rm(tmp)
  }
  output = readRDS("Scenario2.rds")
  output = list(inbredMean=rbind(output$inbredMean,inbredMean),
                inbredVar=rbind(output$inbredVar,inbredVar),
                hybridMean=rbind(output$hybridMean,hybridMean),
                hybridVar=rbind(output$hybridVar,hybridVar),
                hybridCorr=rbind(output$hybridCorr,hybridCorr))
  saveRDS(output,"Scenario2.rds")
}

sc1 = readRDS("Scenario1.rds")
sc2 = readRDS("Scenario2.rds")

#Inbred Mean
plot(-19:20,colMeans(sc2$inbredMean),type="l",col=2,
     main="Inbred Yield",xlab="Year",ylab="Yield (bu/ac)")
lines(-19:20,colMeans(sc1$inbredMean))

#Hybrid Mean
plot(-19:20,colMeans(sc2$hybridMean),type="l",col=2,
     main="Hybrid Yield",xlab="Year",ylab="Yield (bu/ac)")
lines(-19:20,colMeans(sc1$hybridMean))

#Inbred Variance
plot(-19:20,sqrt(colMeans(sc2$inbredVar)),type="l",col=2,
     main="Inbred Standard Deviation",xlab="Year",ylab="Yield (bu/ac)")
lines(-19:20,sqrt(colMeans(sc1$inbredVar)))

#Hybrid Variance
plot(-19:20,sqrt(colMeans(sc2$hybridVar)),type="l",col=2,
     main="Hybrid Standard Deviation",xlab="Year",ylab="Yield (bu/ac)")
lines(-19:20,sqrt(colMeans(sc1$hybridVar)))

#Hybrid Corr
plot(-19:20,sqrt(colMeans(sc2$hybridCorr)),type="l",col=2,
     main="GCA vs per se",xlab="Year",ylab="Correlation")
lines(-19:20,sqrt(colMeans(sc1$hybridCorr)))

