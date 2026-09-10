## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(AlphaSimR)

## -----------------------------------------------------------------------------
set.seed(1234)

founderPop = runMacs(nInd = 500, nChr = 10, segSites = 100)

SP = SimParam$new(founderPop)

## ----include=FALSE------------------------------------------------------------
SP$nThreads = 1L

## -----------------------------------------------------------------------------
# Keep the QTL and the SNP chip on separate segregating sites
SP$restrSegSites(minQtlPerChr = 20, minSnpPerChr = 50, overlap = FALSE)

SP$addTraitAD(nQtlPerChr = 20, meanDD = 0.4)
SP$addSnpChip(nSnpPerChr = 50)
SP$setVarE(h2 = 0.4)

pop = newPop(founderPop, simParam = SP)

## -----------------------------------------------------------------------------
# Two batches of records, phenotyped separately
batch1 = setPheno(pop[1:200], reps = 1, fixEff = 1L, simParam = SP)
batch2 = setPheno(pop[201:400], reps = 3, fixEff = 2L, simParam = SP)

trainPop = c(batch1, batch2)

# Held back for validation
predPop = pop[401:500]

nInd(trainPop)

## -----------------------------------------------------------------------------
table(trainPop@fixEff)

## -----------------------------------------------------------------------------
RRBLUPMemUse(nInd = 2000, nMarker = 10000)
RRBLUPMemUse(nInd = 2000, nMarker = 10000, model = "SCA")

## -----------------------------------------------------------------------------
ans = RRBLUP(trainPop, simParam = SP)

class(ans)

## -----------------------------------------------------------------------------
ans@Vu
ans@Ve

## -----------------------------------------------------------------------------
length(ans@gv)
length(ans@bv)
length(ans@female)

## -----------------------------------------------------------------------------
predPop = setEBV(predPop, ans, simParam = SP)

head(ebv(predPop))

## -----------------------------------------------------------------------------
# Breeding values relative to the prediction population
predPop = setEBV(predPop, ans, value = "bv", targetPop = predPop,
                 append = TRUE, simParam = SP)

cor(ebv(predPop)[, 1], ebv(predPop)[, 2])

## -----------------------------------------------------------------------------
# Accuracy in the prediction population
cor(ebv(predPop)[, 1], gv(predPop)[, 1])
cor(ebv(predPop)[, 1], bv(predPop, simParam = SP)[, 1])

## -----------------------------------------------------------------------------
trainPop = setEBV(trainPop, ans, simParam = SP)

# Inflated: these individuals trained the model
cor(ebv(trainPop)[, 1], gv(trainPop)[, 1])

## -----------------------------------------------------------------------------
selected = selectInd(trainPop, nInd = 50, use = "ebv", simParam = SP)

meanG(selected)
meanG(trainPop)

## -----------------------------------------------------------------------------
newPop = selectCross(trainPop, nInd = 50, nCrosses = 200,
                     use = "ebv", simParam = SP)

## -----------------------------------------------------------------------------
ans2 = RRBLUP2(trainPop, Vu = ans@Vu[1, 1], Ve = ans@Ve[1, 1],
               useEM = FALSE, simParam = SP)

predPop = setEBV(predPop, ans2, simParam = SP)

cor(ebv(predPop)[, 1], gv(predPop)[, 1])

## -----------------------------------------------------------------------------
geno = pullSnpGeno(trainPop, snpChip = 1, simParam = SP)

dim(geno)
geno[1:5, 1:5]

## -----------------------------------------------------------------------------
head(getSnpMap(snpChip = 1, simParam = SP))

## ----eval=FALSE---------------------------------------------------------------
# # 'pred' is a vector of predictions from an external package,
# # in the same order as the individuals in predPop
# predPop@ebv = as.matrix(pred)

