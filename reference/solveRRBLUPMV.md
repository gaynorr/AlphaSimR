# Solve Multivariate RR-BLUP

Solves a multivariate mixed model of form \\Y=X\beta+Mu+e\\

## Usage

``` r
solveRRBLUPMV(Y, X, M, maxIter = 1000L, tol = 1e-06)
```

## Arguments

- Y:

  a matrix with n rows and q columns

- X:

  a matrix with n rows and x columns

- M:

  a matrix with n rows and m columns

- maxIter:

  maximum number of iteration

- tol:

  tolerance for convergence
