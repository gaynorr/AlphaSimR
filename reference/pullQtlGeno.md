# Pull QTL genotypes

Retrieves QTL genotype data

## Usage

``` r
pullQtlGeno(pop, trait = 1, chr = NULL, asRaw = FALSE, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- trait:

  an integer. Indicates which trait's QTL genotypes to retrieve.

- chr:

  a vector of chromosomes to retrieve. If NULL, all chromosome are
  retrieved.

- asRaw:

  return in raw (byte) format

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns a matrix of QTL genotypes.

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
pullQtlGeno(pop, simParam=SP)
#>    1_2 1_4 1_6 1_7 1_8 1_9 1_12 1_13 1_14 1_15
#> 1    2   1   2   1   1   1    0    1    2    1
#> 2    1   2   1   1   0   2    1    0    2    1
#> 3    0   2   0   1   1   2    1    1    0    1
#> 4    2   1   2   2   0   1    1    1    0    1
#> 5    1   1   0   1   0   2    1    0    0    0
#> 6    2   2   2   1   1   1    1    0    1    1
#> 7    1   1   1   1   1   2    1    1    0    2
#> 8    2   1   2   1   1   0    0    2    1    2
#> 9    1   1   1   0   1   1    1    1    1    2
#> 10   1   1   0   2   1   1    0    2    1    1
```
