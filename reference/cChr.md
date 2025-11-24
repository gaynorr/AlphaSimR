# Combine MapPop chromosomes

Merges the chromosomes of multiple
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
or
[`NamedMapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/NamedMapPop-class.md)
objects. Each MapPop must have the same number of chromosomes

## Usage

``` r
cChr(...)
```

## Arguments

- ...:

  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
  or
  [`NamedMapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/NamedMapPop-class.md)
  objects to be combined

## Value

Returns an object of
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

## Examples

``` r
pop1 = quickHaplo(nInd=10, nChr=1, segSites=10)
pop2 = quickHaplo(nInd=10, nChr=1, segSites=10)

combinedPop = cChr(pop1, pop2)
```
