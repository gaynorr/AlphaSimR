# Genetic value

A wrapper for accessing the gv slot

## Usage

``` r
gv(pop)
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
gv(pop)
#>            Trait1
#>  [1,] -0.55017237
#>  [2,]  1.73862445
#>  [3,]  1.24166483
#>  [4,] -0.44588830
#>  [5,]  1.15386227
#>  [6,] -1.10455047
#>  [7,]  0.51983994
#>  [8,] -0.06286323
#>  [9,] -1.61759726
#> [10,] -0.87291985
```
