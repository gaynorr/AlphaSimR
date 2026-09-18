## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(AlphaSimR)

## -----------------------------------------------------------------------------
founderPop = quickHaplo(nInd = 100, nChr = 2, segSites = 20)
SP = SimParam$new(founderPop)
SP$addTraitAD(nQtlPerChr = 10, meanDD = 0.5)
pop = newPop(founderPop, simParam = SP)

ans = genParam(pop, simParam = SP)

head(ans$alpha[[1]])    # Observed frequencies
head(ans$alpha_HW[[1]]) # Hardy-Weinberg frequencies

## -----------------------------------------------------------------------------
founderPop = quickHaplo(nInd = 100, nChr = 2, segSites = 20)
SP = SimParam$new(founderPop)
SP$addTraitADE(nQtlPerChr = 10, meanDD = 0.5, relAA = 0.2)
pop = newPop(founderPop, simParam = SP)

ans = genParam(pop, simParam = SP)

# G = mu + A + D + AA, for every individual
resid = ans$gv[,1] - (ans$mu[1] + ans$bv[,1] + ans$dd[,1] + ans$aa[,1])
range(resid)

## -----------------------------------------------------------------------------
# The full variance identity
varG = ans$varG[1,1]
parts = ans$varA[1,1] + ans$varD[1,1] + ans$varAA[1,1] +
  2*(ans$covAD_L[1] + ans$covAAA_L[1] + ans$covDAA_L[1])

c(varG = varG, parts = parts)

## -----------------------------------------------------------------------------
set.seed(1234)

founderPop = quickHaplo(nInd = 1000, nChr = 10, segSites = 20)
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = 10)
SP$setVarE(h2 = 0.5)
pop = newPop(founderPop, simParam = SP)

sel = selectInd(pop, nInd = 100, simParam = SP)
prog = randCross(sel, nCrosses = 1000, simParam = SP)

before = genParam(pop, simParam = SP)
after = genParam(prog, simParam = SP)

round(c(varA_before = before$varA[1,1],
        varA_after = after$varA[1,1],
        genic_before = before$genicVarA[1],
        genic_after = after$genicVarA[1],
        covL_before = before$covA_L[1],
        covL_after = after$covA_L[1]), 3)

## -----------------------------------------------------------------------------
founderPop = quickHaplo(nInd = 100, nChr = 2, segSites = 20, ploidy = 4L)
SP = SimParam$new(founderPop)
SP$addTraitAD(nQtlPerChr = 10, meanDD = 0.5)
pop = newPop(founderPop, simParam = SP)

ans4 = genParam(pop, simParam = SP)
round(c(varA = ans4$varA[1,1],
        varD = ans4$varD[1,1],
        genicVarA = ans4$genicVarA[1],
        genicVarD = ans4$genicVarD[1]), 3)

