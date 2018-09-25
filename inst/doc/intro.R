## ----eval=FALSE----------------------------------------------------------
#  founderPop = runMacs(nInd=1000, nChr=10, segSites=1000)

## ----eval=FALSE----------------------------------------------------------
#  SP = SimParam$new(founderPop)

## ----eval=FALSE----------------------------------------------------------
#  SP$addTraitA(nQtlPerChr=1000)

## ----eval=FALSE----------------------------------------------------------
#  SP$setGender("yes_sys")

## ----eval=FALSE----------------------------------------------------------
#  pop = newPop(founderPop)

## ----eval=FALSE----------------------------------------------------------
#  genMean = meanG(pop)

## ----eval=FALSE----------------------------------------------------------
#  for(generation in 1:20){
#    pop = selectCross(pop=pop, nFemale=500, nMale=25, use="gv", nCrosses=1000)
#    genMean = c(genMean, meanG(pop))
#  }

## ----eval=FALSE----------------------------------------------------------
#  plot(0:20, genMean, xlab="Generation", ylab="Mean Genetic Value")

## ------------------------------------------------------------------------
library(AlphaSimR)

# Creating Founder Haplotypes
founderPop = runMacs(nInd=1000, nChr=10, segSites=1000)

# Setting Simulation Parameters
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr=1000)
SP$setGender("yes_sys")

# Modeling the Breeding Program
pop = newPop(founderPop)
genMean = meanG(pop)
for(generation in 1:20){
  pop = selectCross(pop=pop, nFemale=500, nMale=25, use="gv", nCrosses=1000)
  genMean = c(genMean, meanG(pop))
}

# Examining the Results
plot(0:20, genMean, xlab="Generation", ylab="Mean Genetic Value")

