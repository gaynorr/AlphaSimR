## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(AlphaSimR)

## -----------------------------------------------------------------------------
founderPop = quickHaplo(nInd = 10, nChr = 2, segSites = 5, genLen = 1)
SP = SimParam$new(founderPop)

head(getGenMap(SP))

## -----------------------------------------------------------------------------
SP$setRecombRatio(2) # Twice as much recombination in females

sapply(SP$femaleMap, max)
sapply(SP$maleMap, max)

## -----------------------------------------------------------------------------
SP$v # Interference parameter, 2.6 approximates Kosambi
SP$p # Proportion of non-interfering crossovers

## -----------------------------------------------------------------------------
SP$p = 0.1 # 10% of crossovers from a non-interfering pathway

## -----------------------------------------------------------------------------
SP$quadProb # Probability of quadrivalent pairing, default is no quadrivalents

## -----------------------------------------------------------------------------
set.seed(5678)

founderPop = quickHaplo(nInd = 10, nChr = 1, segSites = 20, genLen = 1)

SP = SimParam$new(founderPop)
SP$setTrackRec(TRUE)

basePop = newPop(founderPop, simParam = SP)
pop = randCross(basePop, nCrosses = 4, simParam = SP)

## -----------------------------------------------------------------------------
SP$recHist[[pop@id[1]]][[1]][[1]] # First individual, chromosome 1, first copy

## -----------------------------------------------------------------------------
pullIbdHaplo(pop, simParam = SP)

## -----------------------------------------------------------------------------
dh = makeDH(pop, nDH = 1, simParam = SP)

pullIbdHaplo(dh, simParam = SP)

