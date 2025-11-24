# Self individuals

Creates selfed progeny from each individual in a population. Only works
when sexes is "no".

## Usage

``` r
self(pop, nProgeny = 1, parents = NULL, keepParents = TRUE, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- nProgeny:

  total number of selfed progeny per individual

- parents:

  an optional vector of indices for allowable parents

- keepParents:

  should previous parents be used for mother and father.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create population
pop = newPop(founderPop, simParam=SP)

#Self pollinate each individual
pop2 = self(pop, simParam=SP)
```
