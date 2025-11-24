# Usefulness criterion

Calculates the usefulness criterion

## Usage

``` r
usefulness(
  pop,
  trait = 1,
  use = "gv",
  p = 0.1,
  selectTop = TRUE,
  simParam = NULL,
  ...
)
```

## Arguments

- pop:

  and object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  or
  [`HybridPop-class`](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.md)

- trait:

  the trait for selection. Either a number indicating a single trait or
  a function returning a vector of length nInd.

- use:

  select on genetic values (`gv`, default), estimated breeding values
  (`ebv`), breeding values (`bv`), or phenotypes (`pheno`)

- p:

  the proportion of individuals selected

- selectTop:

  selects highest values if true. Selects lowest values if false.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

- ...:

  additional arguments if using a function for trait

## Value

Returns a numeric value

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)

#Determine usefulness of population
usefulness(pop, simParam=SP)
#> [1] 1

#Should be equivalent to GV of best individual
max(gv(pop))
#> [1] 1
```
