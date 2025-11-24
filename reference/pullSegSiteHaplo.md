# Pull seg site haplotypes

Retrieves haplotype data for all segregating sites

## Usage

``` r
pullSegSiteHaplo(
  pop,
  haplo = "all",
  chr = NULL,
  asRaw = FALSE,
  simParam = NULL
)
```

## Arguments

- pop:

  an object of
  [`RawPop-class`](https://gaynorr.github.io/AlphaSimR/reference/RawPop-class.md)
  or
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

- haplo:

  either "all" for all haplotypes or an integer for a single set of
  haplotypes. Use a value of 1 for female haplotypes and a value of 2
  for male haplotypes in diploids.

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

Returns a matrix of haplotypes

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
pullSegSiteHaplo(pop, simParam=SP)
#>      1_1 1_2 1_3 1_4 1_5 1_6 1_7 1_8 1_9 1_10 1_11 1_12 1_13 1_14 1_15
#> 1_1    0   0   0   0   1   1   1   0   0    1    1    0    1    0    1
#> 1_2    0   1   0   0   0   0   0   0   0    1    1    1    1    0    0
#> 2_1    1   1   1   0   0   0   0   1   0    0    1    1    0    0    1
#> 2_2    0   0   0   0   1   1   1   1   0    0    1    0    1    0    0
#> 3_1    1   0   0   1   0   1   0   1   0    0    0    1    1    1    1
#> 3_2    0   1   0   1   0   0   0   0   0    1    1    1    0    0    0
#> 4_1    0   1   0   0   1   0   0   0   1    0    1    1    0    1    0
#> 4_2    0   0   0   0   0   1   1   1   1    1    0    1    0    1    1
#> 5_1    0   0   1   1   0   0   1   0   1    1    0    0    0    0    1
#> 5_2    0   1   1   1   0   0   1   1   1    1    1    1    0    0    0
#> 6_1    1   1   0   0   0   1   0   0   0    1    0    0    1    0    0
#> 6_2    0   1   0   0   1   0   0   0   0    1    0    1    0    0    0
#> 7_1    0   1   0   1   1   0   1   0   1    1    0    1    0    0    0
#> 7_2    0   0   0   0   1   0   1   0   1    0    0    1    1    1    1
#> 8_1    0   0   1   0   1   0   1   0   0    1    0    1    0    0    0
#> 8_2    1   1   1   1   1   0   0   1   1    1    1    0    0    1    1
#> 9_1    0   1   0   0   1   0   0   1   0    0    0    0    1    1    0
#> 9_2    0   1   1   0   0   1   0   1   0    0    1    1    0    1    1
#> 10_1   1   1   0   1   1   1   1   1   1    0    0    1    1    0    1
#> 10_2   1   1   1   0   1   1   0   0   1    0    1    0    1    0    1
```
