## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(AlphaSimR)

## -----------------------------------------------------------------------------
# Founders and a trait
founderPop = quickHaplo(nInd = 100, nChr = 1, segSites = 1)
SP = SimParam$new(founderPop)
SP$addTraitAD(nQtlPerChr = 1, mean = 10, meanDD = 0.75)

# Parameters of the genotypic parameterization
SP$traits[[1]]@intercept
SP$traits[[1]]@addEff
SP$traits[[1]]@domEff

## -----------------------------------------------------------------------------
# Simulate a population
pop = newPop(founderPop)

# QTL genotypes
popQtlGenos = pullQtlGeno(pop)
x_A = popQtlGenos - 1
x_D = popQtlGenos * (2 - popQtlGenos)

# Genotypic parameterization of genetic values
tmp = cbind(popQtlGenos, gv(pop),
            SP$traits[[1]]@intercept,
            x_A, x_A * SP$traits[[1]]@addEff,
            x_D, x_D * SP$traits[[1]]@domEff)
colnames(tmp) = c("genotype", "gv", "mu0", "x_A", "x_A_a", "x_D", "x_D_d")
head(tmp)
# ... genetic value for the individual i
id = 1 # change to inspect other individuals
tmp[id, "gv"]
# ... and its genotypic parameterization
tmp[id, "mu0"] + tmp[id, "x_A_a"] + tmp[id, "x_D_d"]

## -----------------------------------------------------------------------------
# Observed genotype frequencies
tmp = table(popQtlGenos)
(fObs = prop.table(tmp))

# Allele frequencies and Hardy-Weinberg genotype frequencies
(tmp = mean(popQtlGenos) / 2)
(fHW = c((1 - tmp)^2, 2 * tmp * (1 - tmp), tmp^2))

# Deviation from Hardy-Weinberg equilibrium frequencies
fObs - fHW

# Calculate the quantities
ans = genParam(pop)

# Alpha under
head(ans$alpha[[1]])    # Observed frequencies
head(ans$alpha_HW[[1]]) # Hardy-Weinberg frequencies

## -----------------------------------------------------------------------------
# Quantities from
# genotypic parameterization of genetic values and
# breeding value parameterization of genetic values
tmp = cbind(popQtlGenos, ans$gv[, 1],
            ans$gv_mu[1], ans$gv_a[, 1], ans$gv_d[, 1],
            ans$mu[1], ans$bv[, 1], ans$dd[, 1])
colnames(tmp) = c("genotype", "gv",
                  "mu0", "x_A_a", "x_D_d",
                  "meanGv", "bv", "dd")
head(tmp)
# ... genetic value for the individual i
id = 2 # change to inspect other individuals
tmp[id, "gv"]
# ... and its breeding value parameterization
tmp[id, "meanGv"] + tmp[id, "bv"] + tmp[id, "dd"]

# Mean of breeding values and deviations is (numerical) zero
mean(ans$bv)
mean(ans$dd)
# ... while this is not the case for genotypic parameterization
mean(ans$gv_a)
mean(ans$gv_d)

## -----------------------------------------------------------------------------
gvs = ans$gv
# also gv(pop)
meanGv = ans$mu
# also mean(gv(pop))
# also meanG(pop)

dosage = popQtlGenos
meanDosage = mean(dosage)
x = dosage - meanDosage

# alpha from observed genotypes
cov(gvs, x) / var(x)
# also ans$alpha[[1]]

bvs = ans$bv
# also bv(pop)
# also x %*% ans$alpha[[1]]

## -----------------------------------------------------------------------------
rangeY = range(c(gvs, meanGv + bvs))
# ... genetic values (black dots)
plot(gvs ~ dosage,
     xlab = "Genotype dosage", ylab = "Genetic value",
     ylim = rangeY, pch = 21, bg = "black")
# ... least-squares line
abline(a = meanGv - ans$alpha[[1]] * meanDosage, b = ans$alpha[[1]])
# ... breeding values (white dots on the line)
points(meanGv + bvs ~ dosage, pch = 21, bg = "white")
# ... deviations between genetic values and breeding values
segments(x0 = dosage, x1 = dosage,
         y0 = gvs, y1 = meanGv + bvs)
# ... mean of genetic values and dosages (plus symbol)
points(meanGv ~ mean(dosage), pch = "+")

## -----------------------------------------------------------------------------
founderPop = quickHaplo(nInd = 100, nChr = 1, segSites = 3)
SP = SimParam$new(founderPop)
SP$addTraitAD(nQtlPerChr = 3, mean = 10, meanDD = 0.75)
SP$traits
pop = newPop(founderPop)
popQtlGenos = pullQtlGeno(pop)

# Observed genotype frequencies
tmp = apply(popQtlGenos, 2, table)
(fObs = apply(tmp, 2, prop.table))

# Allele frequencies and Hardy-Weinberg genotype frequencies
(tmp = apply(popQtlGenos, 2, mean) / 2)
(fHW = rbind((1 - tmp)^2, 2 * tmp * (1 - tmp), tmp^2))

# Deviation from Hardy-Weinberg equilibrium frequencies
fObs - fHW

# Calculate the quantities
ans = genParam(pop)

# Alphas under
head(ans$alpha[[1]])    # Observed frequencies
head(ans$alpha_HW[[1]]) # Hardy-Weinberg frequencies

# Quantities from
# breeding value parameterization of genetic values and
# genotypic parameterization of genetic values
tmp = cbind(ans$gv[, 1],
            ans$gv_mu[1], ans$gv_a[, 1], ans$gv_d[, 1],
            ans$mu[1], ans$bv[, 1], ans$dd[, 1])
colnames(tmp) = c("gv",
                  "mu0", "x_A_a", "x_D_d",
                  "meanGv", "bv", "dd")
head(tmp)

## -----------------------------------------------------------------------------
founderPop = quickHaplo(nInd = 100, nChr = 2, segSites = 20)
SP = SimParam$new(founderPop)
SP$addTraitADE(nQtlPerChr = 10, mean = 10, meanDD = 0.75, relAA = 0.2)
pop = newPop(founderPop)
ans = genParam(pop)
tmp = cbind(ans$gv[, 1],
            ans$gv_mu[1], ans$gv_a[, 1], ans$gv_d[, 1], ans$gv_aa[, 1],
            ans$mu[1], ans$bv[, 1], ans$dd[, 1], ans$aa[, 1])
colnames(tmp) = c("gv",
                  "mu0", "x_A_a", "x_D_d", "x_AA_aa",
                  "meanGv", "bv", "dd", "aa")
head(tmp)

# G = mu + A + D + AA, for every individual
resid = ans$gv[, 1] - (ans$mu[1] + ans$bv[, 1] + ans$dd[, 1] + ans$aa[, 1])
max(abs(resid)) # numerical zero for every individual

## -----------------------------------------------------------------------------
tmp = cbind(ans$gv[, 1],
            ans$gv_mu[1], ans$gv_a[, 1], ans$gv_d[, 1], ans$gv_aa[, 1], ans$gv_n[, 1],
            ans$mu[1], ans$bv[, 1], ans$dd[, 1], ans$aa[, 1], ans$nd[, 1])
colnames(tmp) = c("gv",
                  "mu0", "x_A_a", "x_D_d", "x_AA_aa", "x_N_n",
                  "meanGv", "bv", "dd", "aa", "nd")
head(tmp)

# G = mu + A + D + AA = mu + A + N, for every individual
resid = (ans$dd[, 1] + ans$aa[, 1]) - ans$nd[, 1]
max(abs(resid)) # numerical zero for every individual

## -----------------------------------------------------------------------------
# Genetic variance
# = variance of genetic values
varG(pop)
k = (nInd(pop) - 1) / nInd(pop)
# also var(gv(pop)) * k
# also popVar(gv(pop))
# also ans$varG

var(gv(pop)) # sample variance, while AlphaSimR reports population variance!

# Difference between population and sample variance in a small population
nSmall = 10
# ... population variance
popVar(gv(pop[1:nSmall]))
# ... sample variance
var(gv(pop[1:nSmall]))

## -----------------------------------------------------------------------------
# Additive genetic variance
# = variance of additive genetic (breeding) values
varA(pop)
# also var(bv(pop)) * k
# also popVar(bv(pop))
# also ans$varA

# Dominance genetic variance
# = variance of dominance deviations
varD(pop)
# also var(dd(pop)) * k
# also popVar(dd(pop))
# also ans$varD

# Epistatic genetic variance
# = variance of epistatic deviations
varAA(pop)
# also var(aa(pop)) * k
# also popVar(aa(pop))
# also ans$varAA

# Non-additive genetic variance
# = variance of non-additive deviations
varN(pop)
# also var(nd(pop)) * k
# also popVar(nd(pop))
# also ans$varN

# Covariances (we describe these covariance in a separate section below)
# Because G = A + D + ... = A + N then
# varG = varA + varD + 2covAD + ... = varA + varN + 2covAN
(tmp = popVar(cbind(bv(pop), dd(pop), aa(pop), nd(pop))))
# also cov(bv(pop), dd(pop)) * k
# also cov(bv(pop), aa(pop)) * k
# also cov(dd(pop), aa(pop)) * k
# also cov(bv(pop), nd(pop)) * k
varG(pop)
sum(tmp[1:3, 1:3]) # varA + varD + 2covAD + ...
sum(tmp[c(1,4), c(1,4)]) # varA + varN + 2covAN

## -----------------------------------------------------------------------------
set.seed(113)
founderDemo = quickHaplo(nInd = 1000, nChr = 1, segSites = 1)
SPDemo = SimParam$new(founderDemo)
SPDemo$addTraitAD(nQtlPerChr = 1, meanDD = 0.5)
popDemo = newPop(founderDemo, simParam = SPDemo)

# One representative of each dosage, repeated to set exact frequencies
representatives = match(0:2, as.vector(pullQtlGeno(popDemo, simParam = SPDemo)))
stopifnot(!anyNA(representatives))
popHWE = popDemo[rep(representatives, times = c(64, 32, 4))]
popInbred = popDemo[rep(representatives, times = c(80, 0, 20))]

referenceSummary = function(population) {
  g = genParam(population, simParam = SPDemo)
  c(p = mean(pullQtlGeno(population, simParam = SPDemo)) / 2,
    alpha = g$alpha[[1]][1],
    alpha_HW = g$alpha_HW[[1]][1],
    genicVarG = unname(g$genicVarG[1]),
    HW_adjustment = unname(g$covG_HW[1]),
    observed_frequency_VarG = unname(g$genicVarG[1] + g$covG_HW[1]),
    observed_VarG = g$varG[1, 1],
    LD_adjustment = unname(g$covG_L[1]),
    genicVarD = unname(g$genicVarD[1]),
    observed_VarD = g$varD[1, 1])
}
round(rbind(HWE = referenceSummary(popHWE),
            inbred = referenceSummary(popInbred)), 4)

## -----------------------------------------------------------------------------
set.seed(1234)
founderSelection = quickHaplo(nInd = 1000, nChr = 10, segSites = 20)
SPSelection = SimParam$new(founderSelection)
SPSelection$addTraitA(nQtlPerChr = 10)
SPSelection$setVarE(h2 = 0.5)
popSelection = newPop(founderSelection, simParam = SPSelection)

selectedPop = selectInd(popSelection, nInd = 100, simParam = SPSelection)
selectionProgeny = randCross(selectedPop, nCrosses = 1000, simParam = SPSelection)

beforeSelection = genParam(popSelection, simParam = SPSelection)
afterSelection = genParam(selectionProgeny, simParam = SPSelection)

round(c(varA_before = beforeSelection$varA[1,1],
        varA_after = afterSelection$varA[1,1],
        genicA_before = unname(beforeSelection$genicVarA[1]),
        genicA_after = unname(afterSelection$genicVarA[1]),
        covA_L_before = unname(beforeSelection$covA_L[1]),
        covA_L_after = unname(afterSelection$covA_L[1])), 3)

## -----------------------------------------------------------------------------
# The full variance identity
varG = ans$varG[1,1]
varG_via_partsADAA = unname(ans$varA[1,1] + ans$varD[1,1] + ans$varAA[1,1] +
  2*(ans$covAD_L[1] + ans$covAAA_L[1] + ans$covDAA_L[1]))
varG_via_partsAN = unname(ans$varA[1,1] + ans$varN[1,1] +
  2*(ans$covAN_L[1]))

c(varG = varG, varG_via_partsADAA = varG_via_partsADAA, varG_via_partsAN = varG_via_partsAN)

