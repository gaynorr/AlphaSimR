## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(AlphaSimR)

## -----------------------------------------------------------------------------
founderGenomes = runMacs(nInd = 10,
                         nChr = 3,
                         segSites = 4,
                         species = "CATTLE")
# Set simulation parameters
SP = SimParam$new(founderGenomes)
SP$setTrackRec(TRUE)
# Inspect the founder genomes
founderGenomes

## ----include=FALSE------------------------------------------------------------
SP$nThreads = 1L

## -----------------------------------------------------------------------------
pullSegSiteGeno(founderGenomes)

## -----------------------------------------------------------------------------
pullSegSiteHaplo(founderGenomes)

## -----------------------------------------------------------------------------
basePop = newPop(founderGenomes)

mutatedBasePop = mutate(basePop, mutRate = 0.1, simParam = SP)

pullSegSiteGeno(basePop)

pullSegSiteGeno(mutatedBasePop)

## -----------------------------------------------------------------------------
secondGenPop = randCross(pop = basePop, nCrosses = 1, nProgeny = 2, simParam = SP)

secondGenPop

# Collect the IDs of the mother and father
mother = secondGenPop@mother[1]
father = secondGenPop@father[1]

## -----------------------------------------------------------------------------
pullSegSiteHaplo(basePop[mother])

## -----------------------------------------------------------------------------
pullSegSiteHaplo(basePop[father])

## -----------------------------------------------------------------------------
pullSegSiteHaplo(secondGenPop)

## -----------------------------------------------------------------------------
pullIbdHaplo(basePop[mother])
pullIbdHaplo(basePop[father])
pullIbdHaplo(secondGenPop)

## -----------------------------------------------------------------------------
founderGenomes = runMacs(nInd = 10,
                        nChr = 3,
                        segSites = 4,
                        species = "CATTLE")
# Set simulation parameters
SP = SimParam$new(founderGenomes)

## ----include=FALSE------------------------------------------------------------
SP$nThreads = 1L

## -----------------------------------------------------------------------------
SP$addTraitA(nQtlPerChr = 3)
SP$addSnpChip(nSnpPerChr = 1)

## -----------------------------------------------------------------------------
basePop = newPop(founderGenomes, simParam = SP)

# Inspect basePop
basePop

## -----------------------------------------------------------------------------
# From SNP
SNP = pullSnpGeno(pop = basePop, snpChip =1)
print(SNP)

# From marker
location = colnames(SNP)
pullMarkerGeno(basePop, markers = location)

## -----------------------------------------------------------------------------
basePop@gv

## -----------------------------------------------------------------------------
getQtlMap(trait = 1, simParam = SP)

## -----------------------------------------------------------------------------
getGenMap(SP)
getSnpMap(snpChip = 1, simParam = SP)

