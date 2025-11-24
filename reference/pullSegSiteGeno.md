# Pull segregating site genotypes

Retrieves genotype data for all segregating sites

## Usage

``` r
pullSegSiteGeno(pop, chr = NULL, asRaw = FALSE, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`RawPop-class`](https://gaynorr.github.io/AlphaSimR/reference/RawPop-class.md)
  or
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

- chr:

  a vector of chromosomes to retrieve. If NULL, all chromosome are
  retrieved.

- asRaw:

  return in raw (byte) format

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md),
  not used if pop is
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

## Value

Returns a matrix of genotypes

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
pullSegSiteGeno(pop, simParam=SP)
#>    1_1 1_2 1_3 1_4 1_5 1_6 1_7 1_8 1_9 1_10 1_11 1_12 1_13 1_14 1_15
#> 1    1   1   0   2   1   2   1   1   1    1    1    1    2    2    1
#> 2    1   0   0   0   0   1   1   1   0    2    0    0    2    2    1
#> 3    2   0   1   1   1   2   0   1   1    2    2    0    1    1    1
#> 4    1   2   1   2   0   2   2   1   0    0    1    0    1    0    0
#> 5    1   0   2   1   1   1   1   1   2    2    2    1    2    0    1
#> 6    1   0   1   1   1   2   1   1   2    2    2    0    2    1    1
#> 7    1   1   1   1   1   2   0   1   1    0    1    0    1    1    1
#> 8    0   1   1   1   2   2   2   2   1    1    1    2    1    1    1
#> 9    1   1   1   2   1   0   1   2   0    1    2    1    1    1    1
#> 10   1   2   0   1   0   2   2   0   1    1    0    1    0    1    1
```
