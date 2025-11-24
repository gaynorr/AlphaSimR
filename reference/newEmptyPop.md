# Creates an empty population

Creates an empty
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
object with user defined ploidy and other parameters taken from
simParam.

## Usage

``` r
newEmptyPop(ploidy = 2L, simParam = NULL)
```

## Arguments

- ploidy:

  the ploidy of the population

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
with zero individuals

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create empty population
pop = newEmptyPop(simParam=SP)
isPop(pop)
#> [1] TRUE
```
