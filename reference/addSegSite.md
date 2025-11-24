# Add segregating site to MapPop

This function allows for adding a new segregating site with user
supplied genotypes to a MapPop. The position of the site is set using a
genetic map position.

## Usage

``` r
addSegSite(mapPop, siteName, chr, mapPos, haplo)
```

## Arguments

- mapPop:

  an object of
  [`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

- siteName:

  name to give the segregating site

- chr:

  which chromosome to add the site

- mapPos:

  genetic map position of site in Morgans

- haplo:

  haplotypes for the site

## Value

an object of
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

## Examples

``` r
# Creates a populations of 10 outbred individuals
# Their genome consists of 1 chromosome and 2 segregating sites
founderPop = quickHaplo(nInd=10,nChr=1,segSites=2)

# Add a locus a the 0.5 Morgan map position
haplo = matrix(sample(x=0:1, size=20, replace=TRUE), ncol=1)

founderPop2 = addSegSite(founderPop, siteName="x", chr=1, mapPos=0.5, haplo=haplo)

pullSegSiteHaplo(founderPop2)
#>      1_1 x 1_2
#> 1_1    1 0   1
#> 1_2    0 0   0
#> 2_1    1 0   1
#> 2_2    1 0   0
#> 3_1    0 1   0
#> 3_2    0 0   0
#> 4_1    1 1   1
#> 4_2    0 1   0
#> 5_1    0 1   1
#> 5_2    1 0   0
#> 6_1    1 1   0
#> 6_2    1 0   0
#> 7_1    0 1   0
#> 7_2    1 0   1
#> 8_1    0 1   1
#> 8_2    1 0   1
#> 9_1    0 0   1
#> 9_2    1 0   1
#> 10_1   1 1   0
#> 10_2   0 0   1
```
