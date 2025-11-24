# Estimated breeding value

A wrapper for accessing the ebv slot

## Usage

``` r
ebv(pop)
```

## Arguments

- pop:

  a
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  or similar object

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
pop@ebv = matrix(rnorm(pop@nInd), nrow=pop@nInd, ncol=1)
ebv(pop)
#>               [,1]
#>  [1,] -0.009052613
#>  [2,] -0.135856744
#>  [3,] -0.336250668
#>  [4,]  0.357796527
#>  [5,]  0.340498580
#>  [6,]  3.097230704
#>  [7,]  0.677352741
#>  [8,]  0.385994964
#>  [9,] -0.838842632
#> [10,] -0.530027590
```
