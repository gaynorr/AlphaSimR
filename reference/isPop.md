# Test if object is of a Population class

Utilify function to test if object is of a Population class

## Usage

``` r
isPop(x)
```

## Arguments

- x:

  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)
isPop(pop)
#> [1] TRUE
isPop(SP)
#> [1] FALSE
```
