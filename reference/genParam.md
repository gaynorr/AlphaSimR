# Sumarize genetic parameters

Calculates genetic and genic additive and dominance variances for an
object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Usage

``` r
genParam(pop, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

- varA:

  an nTrait by nTrait matrix of additive genetic variances

- varD:

  an nTrait by nTrait matrix of dominance genetic variances

- varAA:

  an nTrait by nTrait matrix of additive-by-additive genetic variances

- varG:

  an nTrait by nTrait matrix of total genetic variances

- genicVarA:

  an nTrait vector of additive genic variances

- genicVarD:

  an nTrait vector of dominance genic variances

- genicVarAA:

  an nTrait vector of additive-by-additive genic variances

- genicVarG:

  an nTrait vector of total genic variances

- covA_HW:

  an nTrait vector of additive covariances due to non-random mating

- covD_HW:

  an nTrait vector of dominance covariances due to non-random mating

- covAA_HW:

  an nTrait vector of additive-by-additive covariances due to non-random
  mating

- covG_HW:

  an nTrait vector of total genic covariances due to non-random mating

- covA_L:

  an nTrait vector of additive covariances due to linkage disequilibrium

- covD_L:

  an nTrait vector of dominance covariances due to linkage
  disequilibrium

- covAA_L:

  an nTrait vector of additive-by-additive covariances due to linkage
  disequilibrium

- covAD_L:

  an nTrait vector of additive by dominance covariances due to linkage
  disequilibrium

- covAAA_L:

  an nTrait vector of additive by additive-by-additive covariances due
  to linkage disequilibrium

- covDAA_L:

  an nTrait vector of dominance by additive-by-additive covariances due
  to linkage disequilibrium

- covG_L:

  an nTrait vector of total genic covariances due to linkage
  disequilibrium

- mu:

  an nTrait vector of trait means

- mu_HW:

  an nTrait vector of expected trait means under random mating

- gv:

  a matrix of genetic values with dimensions nInd by nTraits

- bv:

  a matrix of breeding values with dimensions nInd by nTraits

- dd:

  a matrix of dominance deviations with dimensions nInd by nTraits

- aa:

  a matrix of additive-by-additive epistatic deviations with dimensions
  nInd by nTraits

- gv_mu:

  an nTrait vector of intercepts with dimensions nInd by nTraits

- gv_a:

  a matrix of additive genetic values with dimensions nInd by nTraits

- gv_d:

  a matrix of dominance genetic values with dimensions nInd by nTraits

- gv_aa:

  a matrix of additive-by-additive genetic values with dimensions nInd
  by nTraits

- alpha:

  a list of average allele subsitution effects with length nTraits

- alpha_HW:

  a list of average allele subsitution effects at Hardy-Weinberg
  equilibrium with length nTraits

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)
SP$setVarE(h2=0.5)

#Create population
pop = newPop(founderPop, simParam=SP)
ans = genParam(pop, simParam=SP)
```
