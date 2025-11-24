# Solve Multivariate Model

Solves a multivariate mixed model of form \\Y=X\beta+Zu+e\\

## Usage

``` r
solveMVM(Y, X, Z, K, tol = 1e-06, maxIter = 1000L)
```

## Arguments

- Y:

  a matrix with n rows and q columns

- X:

  a matrix with n rows and x columns

- Z:

  a matrix with n rows and m columns

- K:

  a matrix with m rows and m columns

- tol:

  tolerance for convergence

- maxIter:

  maximum number of iteration
