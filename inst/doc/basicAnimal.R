# Basic animal breeding simulation
# Based on "historical breeding" portion of Jenko et al., 2015
# https://gsejournal.biomedcentral.com/articles/10.1186/s12711-015-0135-3

# Generate initial haplotypes, replace macsPath
macsPath = "/Users/rgaynor/Documents/testPBS/macs"
FOUNDERPOP = runMacs(macs=macsPath,nInd=1000,nChr=10,segSites=1000,
                     inbred=FALSE,species="TEST") #Using TEST instead of CATTLE for speed
rm(macsPath) #No longer needed

# Set simulation parameters
SIMPARAM = createSimulation(maxQtl=1000,maxSnp=0,gender="yes_sys")
# Add a trait with additive effects
SIMPARAM = addTraitA(nQtlPerChr=1000,meanG=0,varG=1)


calves = newPop(FOUNDERPOP)
popMean = meanG(calves)
popVar = varG(calves)
for(i in 1:20){
  bulls = selectMale(calves,25,useGv=TRUE)
  cows = selectFemale(calves,500,useGv=TRUE)
  calves = randCross2(cows,bulls,1000)
  popMean = c(popMean,meanG(calves))
  popVar = c(popVar,varG(bulls))
}
plot(0:20,popMean,type="l",xlab="Generation",ylab="GV",main="Genetic Gain")
plot(0:20,popVar,type="l",xlab="Generation",ylab="GV",main="Genetic Variance")


