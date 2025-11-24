# Additive-by-additive epistatic deviations

Returns additive-by-additive epistatic deviations for all traits

## Usage

``` r
aa(pop, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

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
aa(pop, simParam=SP)
#>       Trait1
#>  [1,]      0
#>  [2,]      0
#>  [3,]      0
#>  [4,]      0
#>  [5,]      0
#>  [6,]      0
#>  [7,]      0
#>  [8,]      0
#>  [9,]      0
#> [10,]      0
```
