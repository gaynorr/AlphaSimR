# Pull SNP haplotypes

Retrieves SNP haplotype data

## Usage

``` r
pullSnpHaplo(
  pop,
  snpChip = 1,
  haplo = "all",
  chr = NULL,
  asRaw = FALSE,
  simParam = NULL
)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- snpChip:

  an integer. Indicates which SNP chip's haplotypes to retrieve.

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
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns a matrix of SNP haplotypes.

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
pullSnpHaplo(pop, simParam=SP)
#>      1_5 1_6 1_9 1_14 1_15
#> 1_1    1   0   0    0    1
#> 1_2    0   1   1    0    0
#> 2_1    0   1   1    1    0
#> 2_2    1   0   0    1    1
#> 3_1    1   0   0    0    0
#> 3_2    1   1   1    0    1
#> 4_1    1   1   1    0    0
#> 4_2    1   0   0    1    1
#> 5_1    1   0   0    0    0
#> 5_2    0   0   0    0    0
#> 6_1    0   0   1    1    0
#> 6_2    0   0   1    0    1
#> 7_1    1   1   0    1    1
#> 7_2    0   0   1    0    1
#> 8_1    1   0   0    1    1
#> 8_2    1   0   1    1    0
#> 9_1    0   1   1    0    1
#> 9_2    1   0   1    0    1
#> 10_1   1   0   1    0    0
#> 10_2   0   1   0    1    1
```
