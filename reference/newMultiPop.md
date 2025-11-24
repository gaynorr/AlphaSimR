# Create new Multi Population

Creates a new
[`MultiPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MultiPop-class.md)
from one or more
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
and/or
[`MultiPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MultiPop-class.md)
objects.

## Usage

``` r
newMultiPop(...)
```

## Arguments

- ...:

  one or more
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  and/or
  [`MultiPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MultiPop-class.md)
  objects.

## Value

Returns an object of
[`MultiPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MultiPop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)
megaPop = newMultiPop(pop=pop)
isMultiPop(megaPop)
#> [1] TRUE
```
