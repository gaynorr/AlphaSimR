# Solve RR-BLUP with EM

Solves a univariate mixed model of form \\y=X\beta+Mu+e\\ using the
Expectation-Maximization algorithm.

## Usage

``` r
solveRRBLUP_EM(Y, X, M, Vu, Ve, tol, maxIter, useEM)
```

## Arguments

- Y:

  a matrix with n rows and 1 column

- X:

  a matrix with n rows and x columns

- M:

  a matrix with n rows and m columns

- Vu:

  initial guess for variance of marker effects

- Ve:

  initial guess for error variance

- tol:

  tolerance for declaring convergence

- maxIter:

  maximum iteration for attempting convergence

- useEM:

  should EM algorithm be used. If false, no estimation of variance
  components is performed. The initial values are treated as true.
