# Population variance

Calculates the population variance matrix as opposed to the sample
variance matrix calculated by [`var`](https://rdrr.io/r/stats/cor.html).
i.e. divides by n instead of n-1

## Usage

``` r
popVar(X)
```

## Arguments

- X:

  an n by m matrix

## Value

an m by m variance-covariance matrix
