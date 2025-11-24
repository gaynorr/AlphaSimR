# Pull IBD haplotypes

Retrieves IBD haplotype data

## Usage

``` r
pullIbdHaplo(pop, chr = NULL, snpChip = NULL, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- chr:

  a vector of chromosomes to retrieve. If NULL, all chromosomes are
  retrieved.

- snpChip:

  an integer indicating which SNP array loci are to be retrieved. If
  NULL, all sites are retrieved.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns a matrix of IBD haplotypes.

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
pullIbdHaplo(pop, simParam=SP)
#>      1_1 1_2 1_3 1_4 1_5 1_6 1_7 1_8 1_9 1_10 1_11 1_12 1_13 1_14 1_15
#> 1_1    1   1   1   1   1   1   1   1   1    1    1    1    1    1    1
#> 1_2    2   2   2   2   2   2   2   2   2    2    2    2    2    2    2
#> 2_1    3   3   3   3   3   3   3   3   3    3    3    3    3    3    3
#> 2_2    4   4   4   4   4   4   4   4   4    4    4    4    4    4    4
#> 3_1    5   5   5   5   5   5   5   5   5    5    5    5    5    5    5
#> 3_2    6   6   6   6   6   6   6   6   6    6    6    6    6    6    6
#> 4_1    7   7   7   7   7   7   7   7   7    7    7    7    7    7    7
#> 4_2    8   8   8   8   8   8   8   8   8    8    8    8    8    8    8
#> 5_1    9   9   9   9   9   9   9   9   9    9    9    9    9    9    9
#> 5_2   10  10  10  10  10  10  10  10  10   10   10   10   10   10   10
#> 6_1   11  11  11  11  11  11  11  11  11   11   11   11   11   11   11
#> 6_2   12  12  12  12  12  12  12  12  12   12   12   12   12   12   12
#> 7_1   13  13  13  13  13  13  13  13  13   13   13   13   13   13   13
#> 7_2   14  14  14  14  14  14  14  14  14   14   14   14   14   14   14
#> 8_1   15  15  15  15  15  15  15  15  15   15   15   15   15   15   15
#> 8_2   16  16  16  16  16  16  16  16  16   16   16   16   16   16   16
#> 9_1   17  17  17  17  17  17  17  17  17   17   17   17   17   17   17
#> 9_2   18  18  18  18  18  18  18  18  18   18   18   18   18   18   18
#> 10_1  19  19  19  19  19  19  19  19  19   19   19   19   19   19   19
#> 10_2  20  20  20  20  20  20  20  20  20   20   20   20   20   20   20
```
