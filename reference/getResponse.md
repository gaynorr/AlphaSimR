# Returns a vector response from a population

Returns a vector response from a population

## Usage

``` r
getResponse(pop, trait, use, simParam = NULL, ...)
```

## Arguments

- pop:

  an object of class Pop or HybirdPop

- trait:

  a vector or custom function

- use:

  a character ("rand", "gv", "ebv", "pheno", or "bv"; note that "bv"
  doesn't work on class HybridPop)

- simParam:

  simulation parameters are only used when use="bv"

- ...:

  are additional arguments passed to trait when trait is a function
