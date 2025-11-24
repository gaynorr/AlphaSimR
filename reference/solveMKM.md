# Solve Multikernel Model

Solves a univariate mixed model with multiple random effects.

## Usage

``` r
solveMKM(y, X, Zlist, Klist, maxIter = 40L, tol = 1e-04)
```

## Arguments

- y:

  a matrix with n rows and 1 column

- X:

  a matrix with n rows and x columns

- Zlist:

  a list of Z matrices

- Klist:

  a list of K matrices

- maxIter:

  maximum number of iteration

- tol:

  tolerance for convergence
