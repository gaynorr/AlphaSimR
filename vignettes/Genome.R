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
founderGenomes = runMacs(nInd = 10,
                        nChr = 3,
                        segSites = 4,
                        species = "CATTLE")
# Set simulation parameters
SP = SimParam$new(founderGenomes)

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

