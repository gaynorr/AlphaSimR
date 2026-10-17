## ----setup, include=FALSE-----------------------------------------------------
library(AlphaSimR)
knitr::opts_chunk$set(collapse=TRUE, comment="#>")

## ----parents-and-crosses------------------------------------------------------
set.seed(1942)
founders = quickHaplo(nInd=12, nChr=2, segSites=40, inbred=TRUE)
SPgca = SimParam$new(founders)
SPgca$nThreads = 1L
SPgca$addTraitAD(nQtlPerChr=20, mean=10, meanDD=0.8, varDD=0.2)
parents = newPop(founders, simParam=SPgca)
females = parents[1:6]
males = parents[7:12]

hybrids = hybridCross(females, males, crossPlan="testcross",
                      returnHybridPop=TRUE, simParam=SPgca)
combining = calcGCA(hybrids, use="gv")

## ----combining-means----------------------------------------------------------
combining$GCAf
combining$GCAm
head(combining$SCA)

## ----centered-combining-ability-----------------------------------------------
M = matrix(NA_real_, nrow=nInd(females), ncol=nInd(males),
           dimnames=list(females@id, males@id))
rowIndex = match(hybrids@mother, females@id)
colIndex = match(hybrids@father, males@id)
M[cbind(rowIndex, colIndex)] = gv(hybrids)[, 1]
stopifnot(!anyNA(M))

crossMean = mean(M)
gcaFemale = rowMeans(M) - crossMean
gcaMale = colMeans(M) - crossMean
expected = crossMean + outer(gcaFemale, gcaMale, "+")
sca = M - expected

round(gcaFemale, 3)
round(gcaMale, 3)
round(sca, 3)

# Verify centering and reconstruction before rounding.
stopifnot(max(abs(rowMeans(sca))) < 1e-10,
          max(abs(colMeans(sca))) < 1e-10,
          max(abs(M - (expected + sca))) < 1e-10)

# Check the relationship to calcGCA's uncentered output.
stopifnot(isTRUE(all.equal(
  unname(rowMeans(M)),
  combining$GCAf$Trait1[match(females@id, combining$GCAf$id)]
)))
stopifnot(isTRUE(all.equal(
  unname(colMeans(M)),
  combining$GCAm$Trait1[match(males@id, combining$GCAm$id)]
)))
crossIds = paste(hybrids@mother, hybrids@father, sep="_")
stopifnot(isTRUE(all.equal(
  gv(hybrids)[, 1],
  combining$SCA$Trait1[match(crossIds, combining$SCA$id)],
  check.attributes=FALSE
)))

## ----tester-comparison--------------------------------------------------------
allTesters = setPhenoGCA(females, testers=males, use="gv",
                         inbred=TRUE, onlyPheno=TRUE, simParam=SPgca)
fewerTesters = setPhenoGCA(females, testers=males[1:3], use="gv",
                           inbred=TRUE, onlyPheno=TRUE, simParam=SPgca)
stopifnot(max(abs(allTesters[, 1] - rowMeans(M))) < 1e-10,
          max(abs(fewerTesters[, 1] - rowMeans(M[, 1:3]))) < 1e-10)

round(cbind(
  mean_all=allTesters[, 1],
  mean_subset=fewerTesters[, 1],
  gca_all=allTesters[, 1] - mean(allTesters[, 1]),
  gca_subset=fewerTesters[, 1] - mean(fewerTesters[, 1])
), 3)

## ----combining-variances------------------------------------------------------
populationVariance = function(x) mean((x - mean(x))^2)
components = c(femaleGCA=populationVariance(gcaFemale),
               maleGCA=populationVariance(gcaMale),
               SCA=populationVariance(as.vector(sca)))
round(c(components, total=sum(components),
        observed=populationVariance(as.vector(M))), 4)
stopifnot(abs(sum(components) -
                populationVariance(as.vector(M))) < 1e-10)

