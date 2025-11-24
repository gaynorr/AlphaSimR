# Select open pollinating plants

This function models selection in an open pollinating plant population.
It allows for varying the percentage of selfing. The function also
provides an option for modeling selection as occuring before or after
pollination.

## Usage

``` r
selectOP(
  pop,
  nInd,
  nSeeds,
  probSelf = 0,
  pollenControl = FALSE,
  trait = 1,
  use = "pheno",
  selectTop = TRUE,
  candidates = NULL,
  simParam = NULL,
  ...
)
```

## Arguments

- pop:

  and object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  or
  [`MultiPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MultiPop-class.md)

- nInd:

  the number of plants to select

- nSeeds:

  number of seeds per plant

- probSelf:

  percentage of seeds expected from selfing. Value ranges from 0 to 1.

- pollenControl:

  are plants selected before pollination

- trait:

  the trait for selection. Either a number indicating a single trait or
  a function returning a vector of length nInd. The function must work
  on a vector or matrix of `use` values as `trait(pop@use, ...)` -
  depending on what `use` is. See the examples and
  [`selIndex`](https://gaynorr.github.io/AlphaSimR/reference/selIndex.md).

- use:

  the selection criterion. Either a character (genetic values "gv",
  estimated breeding values "ebv", breeding values "bv", phenotypes
  "pheno", or randomly "rand") or a function returning a vector of
  length nInd. The function must work on `pop` as `use(pop, trait, ...)`
  or as `trait(pop@use, ...)` depending on what `trait` is. See the
  examples.

- selectTop:

  selects highest values if true. Selects lowest values if false.

- candidates:

  an optional vector of eligible selection candidates.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

- ...:

  additional arguments if using a function for `trait` and `use`

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
or
[`MultiPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MultiPop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)

#Create population
pop = newPop(founderPop, simParam=SP)

#Create new population by selecting the best 3 plant
#Assuming 50% selfing in plants and 10 seeds per plant
pop2 = selectOP(pop, nInd=3, nSeeds=10, probSelf=0.5, simParam=SP)
```
