## ----eval=FALSE---------------------------------------------------------------
# founderPop = quickHaplo(nInd=1000, nChr=10, segSites=1000)

## ----eval=FALSE---------------------------------------------------------------
# SP = SimParam$new(founderPop)

## ----eval=FALSE---------------------------------------------------------------
# SP$addTraitA(nQtlPerChr=1000)

## ----eval=FALSE---------------------------------------------------------------
# SP$setSexes("yes_sys")

## ----eval=FALSE---------------------------------------------------------------
# pop = newPop(founderPop)

## ----eval=FALSE---------------------------------------------------------------
# genMean = meanG(pop)

## ----eval=FALSE---------------------------------------------------------------
# for(generation in 1:20){
#   pop = selectCross(pop=pop, nFemale=500, nMale=25, use="gv", nCrosses=1000)
#   genMean = c(genMean, meanG(pop))
# }

## ----eval=FALSE---------------------------------------------------------------
# plot(0:20, genMean, xlab="Generation", ylab="Mean Genetic Value", type="l")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(AlphaSimR)

## -----------------------------------------------------------------------------
# Creating Founder Haplotypes
founderPop = quickHaplo(nInd=1000, nChr=10, segSites=1000)

# Setting Simulation Parameters
SP = SimParam$new(founderPop)

## ----include=FALSE------------------------------------------------------------
SP$nThreads = 1L

## -----------------------------------------------------------------------------
SP$addTraitA(nQtlPerChr=1000)
SP$setSexes("yes_sys")

# Modeling the Breeding Program
pop = newPop(founderPop)
genMean = meanG(pop)
for(generation in 1:20){
  pop = selectCross(pop=pop, nFemale=500, nMale=25, use="gv", nCrosses=1000)
  genMean = c(genMean, meanG(pop))
}

# Examining the Results
plot(0:20, genMean, xlab="Generation", ylab="Mean Genetic Value", type="l")

