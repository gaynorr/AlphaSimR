# Create new population

Creates an initial
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
from an object of
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
or
[`NamedMapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/NamedMapPop-class.md).
The function is intended for use with output from functions such as
[`runMacs`](https://gaynorr.github.io/AlphaSimR/reference/runMacs.md),
[`newMapPop`](https://gaynorr.github.io/AlphaSimR/reference/newMapPop.md),
or
[`quickHaplo`](https://gaynorr.github.io/AlphaSimR/reference/quickHaplo.md).

## Usage

``` r
newPop(rawPop, simParam = NULL, ...)
```

## Arguments

- rawPop:

  an object of
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
  or
  [`NamedMapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/NamedMapPop-class.md)

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

- ...:

  additional arguments used internally

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Details

Note that `newPop` takes genomes from the `rawPop` and uses them without
recombination! Hence, if you call `newPop(rawPop = founderGenomes)`
twice, you will get two sets of individuals with different id but the
same genomes. To get genetically different sets of individuals you can
subset the `rawPop` input, say first half for one set and the second
half for the other set.

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

#Misc
pop@misc$tmp1 = rnorm(n=2)
pop@misc$tmp2 = rnorm(n=2)

#MiscPop
pop@miscPop$tmp1 = sum(pop@misc$tmp1)
pop@miscPop$tmp2 = sum(pop@misc$tmp2)
```
