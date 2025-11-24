# Lose individuals at random

Samples individuals at random to remove from the population. The user
supplies a probability for the individuals to be removed from the
population.

## Usage

``` r
attrition(pop, p)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- p:

  the expected proportion of individuals that will be lost to attrition.

## Value

an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=100, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create population
pop = newPop(founderPop, simParam=SP)

#Lose an expected 5% of individuals
pop = attrition(pop, p=0.05)
```
