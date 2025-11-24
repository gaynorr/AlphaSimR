# RRBLUP Memory Usage

Estimates the amount of RAM needed to run the
[`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md) and
its related functions for a given training population size. Note that
this function may underestimate total usage.

## Usage

``` r
RRBLUPMemUse(nInd, nMarker, model = "REG")
```

## Arguments

- nInd:

  the number of individuals in the training population

- nMarker:

  the number of markers per individual

- model:

  either "REG", "GCA", or "SCA" for
  [`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md)
  [`RRBLUP_GCA`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP_GCA.md)
  and
  [`RRBLUP_SCA`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP_SCA.md)
  respectively.

## Value

Returns an estimate for the required gigabytes of RAM

## Examples

``` r
RRBLUPMemUse(nInd=1000, nMarker=5000)
#> [1] 0.06412
```
