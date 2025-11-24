# Breeding value

Returns breeding values for all traits

## Usage

``` r
bv(pop, simParam = NULL)
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
bv(pop, simParam=SP)
#>            Trait1
#>  [1,] -0.59514020
#>  [2,]  0.01471392
#>  [3,]  0.80541083
#>  [4,]  1.35860478
#>  [5,] -1.22240126
#>  [6,]  1.84002591
#>  [7,] -1.22946991
#>  [8,]  0.21801518
#>  [9,] -0.54879548
#> [10,] -0.64096376
```
