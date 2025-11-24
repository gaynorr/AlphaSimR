# Pull marker haplotypes

Retrieves haplotype data for user specified loci

## Usage

``` r
pullMarkerHaplo(pop, markers, haplo = "all", asRaw = FALSE, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`RawPop-class`](https://gaynorr.github.io/AlphaSimR/reference/RawPop-class.md)
  or
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

- markers:

  a character vector. Indicates the names of the loci to be retrieved

- haplo:

  either "all" for all haplotypes or an integer for a single set of
  haplotypes. Use a value of 1 for female haplotypes and a value of 2
  for male haplotypes in diploids.

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
SP$setTrackRec(TRUE)

#Create population
pop = newPop(founderPop, simParam=SP)

#Pull haplotype data for first two markers on chromosome one.
#Marker name is consistent with default naming in AlphaSimR.
pullMarkerHaplo(pop, markers=c("1_1","1_2"), simParam=SP)
#>      1_1 1_2
#> 1_1    1   0
#> 1_2    1   0
#> 2_1    1   0
#> 2_2    0   1
#> 3_1    1   1
#> 3_2    1   1
#> 4_1    0   1
#> 4_2    0   1
#> 5_1    0   1
#> 5_2    0   1
#> 6_1    1   0
#> 6_2    1   1
#> 7_1    0   0
#> 7_2    1   0
#> 8_1    1   1
#> 8_2    0   1
#> 9_1    0   1
#> 9_2    0   1
#> 10_1   0   1
#> 10_2   1   0
```
