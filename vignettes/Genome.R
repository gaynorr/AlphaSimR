## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(AlphaSimR)

## -----------------------------------------------------------------------------
founderGenome = runMacs(nInd = 10, 
                        nChr = 3, 
                        segSites = 4, 
                        species = "CATTLE")
# Set simulation parameters
SP = SimParam$new(founderGenome)
SP$setTrackRec(TRUE)
# Inspect the founderGenome
founderGenome


## -----------------------------------------------------------------------------
pullSegSiteGeno(founderGenome)

## -----------------------------------------------------------------------------
pullSegSiteHaplo(founderGenome)

## -----------------------------------------------------------------------------
basePop = newPop(founderGenome)

mutatedBasePop = mutate(basePop, mutRate = 0.1, simParam = SP)

pullSegSiteGeno(basePop)

pullSegSiteGeno(mutatedBasePop)

## -----------------------------------------------------------------------------
secondGenPop = randCross(pop = basePop, nCrosses = 1, nProgeny = 2, simParam = SP)

str(secondGenPop)

# Collect the iid for mother and father
mother = as.integer(secondGenPop@mother[1])
father = as.integer(secondGenPop@father[1])

## -----------------------------------------------------------------------------
pullSegSiteHaplo(basePop[basePop@iid == mother,])

## -----------------------------------------------------------------------------
pullSegSiteHaplo(basePop[basePop@iid == father,])

## -----------------------------------------------------------------------------
pullSegSiteHaplo(secondGenPop)

## -----------------------------------------------------------------------------
pullIbdHaplo(basePop[basePop@iid == mother,])
pullIbdHaplo(basePop[basePop@iid == father,])
pullIbdHaplo(secondGenPop)

## -----------------------------------------------------------------------------
founderGenome = runMacs(nInd = 10, 
                        nChr = 3, 
                        segSites = 4, 
                        species = "CATTLE")
# Set simulation parameters
SP = SimParam$new(founderGenome)

## -----------------------------------------------------------------------------
SP$addTraitA(nQtlPerChr = 3)
SP$addSnpChip(nSnpPerChr = 1)

## -----------------------------------------------------------------------------
basePop = newPop(founderGenome, simParam = SP)

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

