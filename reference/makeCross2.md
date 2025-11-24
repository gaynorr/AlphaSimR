# Make designed crosses

Makes crosses between two populations using a user supplied crossing
plan.

## Usage

``` r
makeCross2(females, males, crossPlan, nProgeny = 1, simParam = NULL)
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

- nProgeny:

  number of progeny per cross

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
pop2 = makeCross2(pop, pop, crossPlan, simParam=SP)
```
