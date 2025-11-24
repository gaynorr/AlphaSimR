# Combine genomes of individuals

This function is designed to model the pairing of gametes. The male and
female individuals are treated as gametes, so the ploidy of newly
created individuals will be the sum of it parents.

## Usage

``` r
mergeGenome(females, males, crossPlan, simParam = NULL)
```

## Arguments

- females:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  for female parents.

- males:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  for male parents.

- crossPlan:

  a matrix with two column representing female and male parents. Either
  integers for the position in population or character strings for the
  IDs.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create population
pop = newPop(founderPop, simParam=SP)

#Cross individual 1 with individual 10
crossPlan = matrix(c(1,10), nrow=1, ncol=2)
pop2 = mergeGenome(pop, pop, crossPlan, simParam=SP)
```
