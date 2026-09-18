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
# 1. Training population: phenotyped and genotyped
trainPop = setPheno(pop[1:400], simParam = SP)

# Selection candidates: genotyped only
candidates = pop[401:500]

# 2. Fit
ans = RRBLUP(trainPop, simParam = SP)

# 3. Predict
candidates = setEBV(candidates, ans, simParam = SP)

# Accuracy, measured on individuals that did not train the model
cor(ebv(candidates)[, 1], gv(candidates)[, 1])

## -----------------------------------------------------------------------------
combined = c(setPheno(pop[1:200], reps = 1, simParam = SP),
             setPheno(pop[201:400], reps = 3, simParam = SP))

nInd(combined)

## -----------------------------------------------------------------------------
ans2 = RRBLUP2(trainPop, Vu = ans@Vu[1, 1], Ve = ans@Ve[1, 1],
               useEM = FALSE, simParam = SP)

candidates = setEBV(candidates, ans2, simParam = SP)

cor(ebv(candidates)[, 1], gv(candidates)[, 1])

## ----eval=FALSE---------------------------------------------------------------
# # Predict every DH and advance the best to the first yield trial
# DH = setEBV(DH, gsModel, simParam = SP)
# 
# PYT = selectWithinFam(DH, famMax, use = "ebv", simParam = SP)
# PYT = selectInd(PYT, nPYT, use = "ebv", simParam = SP)
# PYT = setPheno(PYT, varE = varE, reps = repPYT, simParam = SP)

## ----eval=FALSE---------------------------------------------------------------
# DH = setEBV(DH, gsModel, simParam = SP)
# 
# newParents = selectInd(DH, nNewParents, use = "ebv", simParam = SP)
# 
# # Replace the oldest parents with the newly predicted ones
# Parents = c(Parents[-(1:nNewParents)], newParents)

## ----eval=FALSE---------------------------------------------------------------
# if (year == startTP) {
#   trainPop = c(PYT, AYT, EYT)
# } else if (year <= nBurnin) {
#   # Accumulate
#   trainPop = c(trainPop, PYT, AYT, EYT)
# } else {
#   # Fixed size: drop the oldest year of records
#   nNew = c(PYT, AYT, EYT)@nInd
#   trainPop = c(trainPop[-(1:nNew)], PYT, AYT, EYT)
# }

## -----------------------------------------------------------------------------
gsPop = pop
gsTrainPop = setPheno(pop, simParam = SP)

for (cycle in 1:3) {
  # Fit to everything phenotyped so far
  gsModel = RRBLUP(gsTrainPop, simParam = SP)

  # Predict and select
  gsPop = setEBV(gsPop, gsModel, simParam = SP)
  gsPop = selectCross(gsPop, nInd = 50, nCrosses = 100, nProgeny = 2,
                      use = "ebv", simParam = SP)

  # Phenotype the new generation and add it to the training records
  gsTrainPop = c(gsTrainPop, setPheno(gsPop, simParam = SP))

  cat("cycle", cycle, " meanG =", round(meanG(gsPop), 3),
      " nTrain =", nInd(gsTrainPop), "\n")
}

## ----eval=FALSE---------------------------------------------------------------
# gsModel = RRBLUP_GCA(hybridTrainPop, simParam = SP)
# 
# maleDH = setEBV(maleDH, gsModel, value = "male", simParam = SP)
# femaleDH = setEBV(femaleDH, gsModel, value = "female", simParam = SP)
# 
# newMaleParents = selectInd(maleDH, 10, use = "ebv", simParam = SP)
# newFemaleParents = selectInd(femaleDH, 10, use = "ebv", simParam = SP)

## -----------------------------------------------------------------------------
candidates = setEBV(candidates, ans, value = "bv", targetPop = candidates,
                    append = TRUE, simParam = SP)

cor(ebv(candidates)[, 1], ebv(candidates)[, 2])

## -----------------------------------------------------------------------------
cor(ebv(candidates)[, 1], gv(candidates)[, 1])
cor(ebv(candidates)[, 1], bv(candidates, simParam = SP)[, 1])

## -----------------------------------------------------------------------------
trainPop = setEBV(trainPop, ans, simParam = SP)

# Inflated: these individuals trained the model
cor(ebv(trainPop)[, 1], gv(trainPop)[, 1])

## -----------------------------------------------------------------------------
RRBLUPMemUse(nInd = 5000, nMarker = 2000, model = "fastRRBLUP")
RRBLUPMemUse(nInd = 5000, nMarker = 2000, model = "RRBLUP")
RRBLUPMemUse(nInd = 5000, nMarker = 2000, model = "RRBLUP_SCA")

## -----------------------------------------------------------------------------
class(ans)

ans@Vu
ans@Ve

## -----------------------------------------------------------------------------
length(ans@gv)
length(ans@bv)
length(ans@female)

## -----------------------------------------------------------------------------
batch1 = setPheno(pop[1:200], fixEff = 1L, simParam = SP)
batch2 = setPheno(pop[201:400], fixEff = 2L, simParam = SP)

table(c(batch1, batch2)@fixEff)

## ----eval=FALSE---------------------------------------------------------------
# # Each trial phenotyped in its own environment, and labelled as such
# PYT = setPheno(PYT, reps = repPYT, p = pYear, fixEff = year, simParam = SP)
# AYT = setPheno(AYT, reps = repAYT, p = pYear, fixEff = year, simParam = SP)
# 
# trainPop = c(trainPop, PYT, AYT)

## -----------------------------------------------------------------------------
geno = pullSnpGeno(trainPop, snpChip = 1, simParam = SP)

dim(geno)
geno[1:5, 1:5]

## -----------------------------------------------------------------------------
head(getSnpMap(snpChip = 1, simParam = SP))

## ----eval=FALSE---------------------------------------------------------------
# # 'pred' is a vector of predictions from an external package,
# # in the same order as the individuals in candidates
# candidates@ebv = as.matrix(pred)

