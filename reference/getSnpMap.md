# Get SNP genetic map

Retrieves the genetic map for a given SNP chip.

## Usage

``` r
getSnpMap(snpChip = 1, sex = "A", simParam = NULL)
```

## Arguments

- snpChip:

  an integer. Indicates which SNP chip's map to retrieve.

- sex:

  determines which sex specific map is returned. Options are "A" for
  average map, "F" for female map, and "M" for male map. All options are
  equivalent if not using sex specific maps.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns a data.frame with:

- id:

  Unique identifier for the SNP

- chr:

  Chromosome containing the SNP

- site:

  Segregating site on the chromosome

- pos:

  Genetic map position

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addSnpChip(5)

#Pull SNP map
getSnpMap(snpChip=1, simParam=SP)
#>        id chr site       pos
#> 1_1   1_1   1    1 0.0000000
#> 1_5   1_5   1    5 0.4444444
#> 1_7   1_7   1    7 0.6666667
#> 1_9   1_9   1    9 0.8888889
#> 1_10 1_10   1   10 1.0000000
```
