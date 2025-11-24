# Phenotype

A wrapper for accessing the pheno slot

## Usage

``` r
pheno(pop)
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
pheno(pop)
#>           Trait1
#>  [1,] -0.4640616
#>  [2,] -0.8502050
#>  [3,] -1.5174515
#>  [4,]  1.5898559
#>  [5,]  1.7206612
#>  [6,] -2.8513408
#>  [7,] -0.1504351
#>  [8,]  0.1121931
#>  [9,]  0.7380158
#> [10,]  0.5339829
```
