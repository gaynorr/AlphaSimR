# Solve Multikernel RR-BLUP

Solves a univariate mixed model with multiple random effects.

## Usage

``` r
solveRRBLUPMK(y, X, Mlist, maxIter = 40L)
```

## Arguments

- y:

  a matrix with n rows and 1 column

- X:

  a matrix with n rows and x columns

- Mlist:

  a list of M matrices

- maxIter:

  maximum number of iteration
