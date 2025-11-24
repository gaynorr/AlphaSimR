# Convert a normal (Gaussian) trait to an ordered categorical (threshold) trait

Convert a normal (Gaussian) trait to an ordered categorical (threshold)
trait

## Usage

``` r
asCategorical(
  x,
  p = NULL,
  mean = 0,
  var = 1,
  threshold = c(-Inf, 0, Inf),
  include.lowest = TRUE,
  right = FALSE
)
```

## Arguments

- x:

  matrix, values for one or more traits (if not a matrix, we cast to a
  matrix)

- p:

  NULL, numeric or list, when `NULL` the `threshold` argument takes
  precedence; when numeric, provide a vector of probabilities of
  categories to convert continuous values into categories for a single
  trait (if probabilities do not sum to 1, another category is added and
  a warning is raised); when list, provide a list of numeric
  probabilities - list node with `NULL` will skip conversion for a
  specific trait (see examples); internally `p` is converted to
  `threshold` hence input `threshold` is overwritten

- mean:

  numeric, assumed mean(s) of the normal (Gaussian) trait(s); used only
  when `p` is given

- var:

  numeric, assumed variance(s) of the normal (Gaussian) trait(s); used
  only when `p` is given

- threshold:

  NULL, numeric or list, when numeric, provide a vector of threshold
  values to convert continuous values into categories for a single trait
  (the thresholds specify left-closed and right-opened intervals \[t1,
  t2), which can be changed with `include.lowest` and `right`; ensure
  you add `-Inf` and `Inf` or min and max to cover the whole range of
  values; otherwise you will get `NA` values); when list, provide a list
  of numeric thresholds - list node with `NULL` will skip conversion for
  a specific trait (see examples)

- include.lowest:

  logical, see [`cut`](https://rdrr.io/r/base/cut.html)

- right:

  logical, see [`cut`](https://rdrr.io/r/base/cut.html)

## Value

matrix of values with some traits recorded as ordered categories in the
form of 1:nC with nC being the number of categories.

## Details

If input trait is normal (Gaussian) then this function generates a
categorical trait according to the ordered probit model.

## Examples

``` r
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
SP = SimParam$new(founderPop)
trtMean = c(0, 0)
trtVarG = c(1, 2)
SP$addTraitA(nQtlPerChr = 10, mean = trtMean, var = trtVarG,
             corA = matrix(data = c(1.0, 0.6,
                                    0.6, 1.0), ncol = 2))
pop = newPop(founderPop)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found
trtVarE = c(1, 1)
trtVarP = trtVarG + trtVarE
pop = setPheno(pop, varE = trtVarE)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found
pheno(pop)
#> Error: object 'pop' not found

#Convert a single input trait
asCategorical(x = pheno(pop)[, 1])
#> Error: object 'pop' not found

#Demonstrate threshold argument (in units of pheno SD)
asCategorical(x = pheno(pop)[, 1], threshold = c(-1, 0, 1) * sqrt(trtVarP[1]))
#> Error: object 'pop' not found
asCategorical(x = pheno(pop)[, 1], threshold = c(-Inf, -1, 0, 1, Inf) * sqrt(trtVarP[1]))
#> Error: object 'pop' not found
asCategorical(x = pheno(pop)[, 1], threshold = c(-Inf, 0, Inf))
#> Error: object 'pop' not found

#Demonstrate p argument
asCategorical(x = pheno(pop)[, 1], p = 0.5, var = trtVarP[1])
#> Error: object 'pop' not found
asCategorical(x = pheno(pop)[, 1], p = c(0.5, 0.5), var = trtVarP[1])
#> Error: object 'pop' not found
asCategorical(x = pheno(pop)[, 1], p = c(0.25, 0.5, 0.25), var = trtVarP[1])
#> Error: object 'pop' not found

#Convert multiple input traits (via threshold or p argument)
try(asCategorical(x = pheno(pop)))
#> Error in eval(expr, envir) : object 'pop' not found
asCategorical(x = pheno(pop),
              threshold = list(c(-Inf, 0, Inf),
                               NULL))
#> Error: object 'pop' not found
try(asCategorical(x = pheno(pop), p = c(0.5, 0.5)))
#> Error in eval(expr, envir) : object 'pop' not found
asCategorical(x = pheno(pop),
              p = list(c(0.5, 0.5),
                       NULL),
              mean = trtMean, var = trtVarP)
#> Error: object 'pop' not found

asCategorical(x = pheno(pop),
              threshold = list(c(-Inf, 0, Inf),
                               c(-Inf, -2, -1, 0, 1, 2, Inf) * sqrt(trtVarP[2])))
#> Error: object 'pop' not found
q = c(-2, -1, 0, 1, 2)
p = pnorm(q)
p = c(p[1], p[2]-p[1], p[3]-p[2], p[4]-p[3], p[5]-p[4], 1-p[5])
asCategorical(x = pheno(pop),
              p = list(c(0.5, 0.5),
                       p),
              mean = trtMean, var = trtVarP)
#> Error: object 'pop' not found
```
