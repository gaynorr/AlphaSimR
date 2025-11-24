# Calculate Smith-Hazel weights

Calculates weights for Smith-Hazel index given economice weights and
phenotypic and genotypic variance-covariance matrices.

## Usage

``` r
smithHazel(econWt, varG, varP)
```

## Arguments

- econWt:

  vector of economic weights

- varG:

  the genetic variance-covariance matrix

- varP:

  the phenotypic variance-covariance matrix

## Value

a vector of weight for calculating index values

## Examples

``` r
G = 1.5*diag(2)-0.5
E = diag(2)
P = G+E
wt = c(1,1)
smithHazel(wt, G, P)
#>           [,1]
#> [1,] 0.3333333
#> [2,] 0.3333333
```
