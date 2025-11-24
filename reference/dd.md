# Dominance deviations

Returns dominance deviations for all traits

## Usage

``` r
dd(pop, simParam = NULL)
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
dd(pop, simParam=SP)
#>            Trait1
#>  [1,]  0.27435866
#>  [2,]  0.04265155
#>  [3,] -0.09582656
#>  [4,]  0.43421754
#>  [5,] -0.04657432
#>  [6,]  0.16935884
#>  [7,] -0.40243292
#>  [8,] -0.17625453
#>  [9,]  0.26253349
#> [10,] -0.46203174
```
