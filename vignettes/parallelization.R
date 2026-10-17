## -----------------------------------------------------------------------------
library(AlphaSimR)
getNumThreads()

## -----------------------------------------------------------------------------
founderPop = quickHaplo(nInd = 4, nChr = 2, segSites = 5)
SP = SimParam$new(founderPop)
SP$nThreads

## -----------------------------------------------------------------------------
SP$nThreads = 1L
pop = newPop(founderPop, simParam = SP)
pullSegSiteGeno(pop)
# Note that pullSegSiteGeno() uses internally SP$nThreads,
# when the nThreads argument is not specified.

## -----------------------------------------------------------------------------
pullSegSiteGeno(pop, nThreads = 1L)

## -----------------------------------------------------------------------------
system.time(runMacs(nInd = 50, nChr = 10, nThreads = 1L))
system.time(runMacs(nInd = 50, nChr = 10, nThreads = getNumThreads()))

