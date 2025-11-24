# Get genetic map

Retrieves the genetic map for all loci.

## Usage

``` r
getGenMap(object = NULL, sex = "A")
```

## Arguments

- object:

  where to retrieve the genetic map. Can be an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)
  or
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md).
  If NULL, the function will look for a SimParam object called "SP" in
  your global environment.

- sex:

  determines which sex specific map is returned. Options are "A" for
  average map, "F" for female map, and "M" for male map. All options are
  equivalent if not using sex specific maps or using pulling from a
  MapPop.

## Value

Returns a data.frame with:

- id:

  Unique identifier for locus

- chr:

  Chromosome containing the locus

- pos:

  Genetic map position

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
getGenMap(founderPop)
#>      id chr       pos
#> 1   1_1   1 0.0000000
#> 2   1_2   1 0.1111111
#> 3   1_3   1 0.2222222
#> 4   1_4   1 0.3333333
#> 5   1_5   1 0.4444444
#> 6   1_6   1 0.5555556
#> 7   1_7   1 0.6666667
#> 8   1_8   1 0.7777778
#> 9   1_9   1 0.8888889
#> 10 1_10   1 1.0000000
```
