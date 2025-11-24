# Sample haplotypes from a MapPop

Creates a new
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
from an existing
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
by randomly sampling haplotypes.

## Usage

``` r
sampleHaplo(mapPop, nInd, inbred = FALSE, ploidy = NULL, replace = TRUE)
```

## Arguments

- mapPop:

  the
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)
  used to sample haplotypes

- nInd:

  the number of individuals to create

- inbred:

  should new individuals be fully inbred

- ploidy:

  new ploidy level for organism. If NULL, the ploidy level of the mapPop
  is used.

- replace:

  should haplotypes be sampled with replacement

## Value

an object of
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

## Examples

``` r
founderPop = quickHaplo(nInd=2,nChr=1,segSites=11,inbred=TRUE)
founderPop = sampleHaplo(mapPop=founderPop,nInd=20)
```
