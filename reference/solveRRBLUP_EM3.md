# Solve RR-BLUP with EM and 3 random effects

Solves a univariate mixed model of form
\\y=X\beta+M_1u_1+M_2u_2+M_3u_3+e\\ using the Expectation-Maximization
algorithm.

## Usage

``` r
solveRRBLUP_EM3(Y, X, M1, M2, M3, Vu1, Vu2, Vu3, Ve, tol, maxIter, useEM)
```

## Arguments

- Y:

  a matrix with n rows and 1 column

- X:

  a matrix with n rows and x columns

- M1:

  a matrix with n rows and m1 columns

- M2:

  a matrix with n rows and m2 columns

- M3:

  a matrix with n rows and m3 columns

- Vu1:

  initial guess for variance of the first marker effects

- Vu2:

  initial guess for variance of the second marker effects

- Vu3:

  initial guess for variance of the second marker effects

- Ve:

  initial guess for error variance

- tol:

  tolerance for declaring convergence

- maxIter:

  maximum iteration for attempting convergence

- useEM:

  should EM algorithm be used. If false, no estimation of variance
  components is performed. The initial values are treated as true.
