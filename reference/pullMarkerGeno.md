# Pull marker genotypes

Retrieves genotype data for user specified loci

## Usage

``` r
pullMarkerGeno(pop, markers, asRaw = FALSE, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`RawPop-class`](https://gaynorr.github.io/AlphaSimR/reference/RawPop-class.md)
  or
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

- markers:

  a character vector. Indicates the names of the loci to be retrieved.

- asRaw:

  return in raw (byte) format

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md),
  not used if pop is
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

## Value

Returns a matrix of genotypes.

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=15)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$addSnpChip(5)

#Create population
pop = newPop(founderPop, simParam=SP)

#Pull genotype data for first two markers on chromosome one.
#Marker name is consistent with default naming in AlphaSimR.
pullMarkerGeno(pop, markers=c("1_1","1_2"), simParam=SP)
#>    1_1 1_2
#> 1    1   0
#> 2    1   2
#> 3    1   2
#> 4    0   1
#> 5    0   0
#> 6    0   1
#> 7    1   2
#> 8    1   1
#> 9    0   0
#> 10   1   0
```
