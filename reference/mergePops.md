# Merge list of populations

Rapidly merges a list of populations into a single population

## Usage

``` r
mergePops(popList)
```

## Arguments

- popList:

  a list containing
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  elements or a
  [`MultiPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MultiPop-class.md)

## Value

Returns a
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create a list of populations and merge list
pop = newPop(founderPop, simParam=SP)
pop@misc$tmp = rnorm(n=10)
pop@misc$tmp2 = rnorm(n=10)

popList = list(pop, pop)
pop2 = mergePops(popList)
```
