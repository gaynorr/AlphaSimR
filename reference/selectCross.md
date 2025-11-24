# Select and randomly cross

This is a wrapper that combines the functionalities of
[`randCross`](https://gaynorr.github.io/AlphaSimR/reference/randCross.md)
and
[`selectInd`](https://gaynorr.github.io/AlphaSimR/reference/selectInd.md).
The purpose of this wrapper is to combine both selection and crossing in
one function call that minimized the amount of intermediate populations
created. This reduces RAM usage and simplifies code writing. Note that
this wrapper does not provide the full functionality of either function.

## Usage

``` r
selectCross(
  pop,
  nInd = NULL,
  nFemale = NULL,
  nMale = NULL,
  nCrosses,
  nProgeny = 1,
  trait = 1,
  use = "pheno",
  selectTop = TRUE,
  simParam = NULL,
  ...,
  balance = TRUE
)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- nInd:

  the number of individuals to select. These individuals are selected
  without regards to sex and it supercedes values for nFemale and nMale.
  Thus if the simulation uses sexes, it is likely better to leave this
  value as NULL and use nFemale and nMale instead.

- nFemale:

  the number of females to select. This value is ignored if nInd is set.

- nMale:

  the number of males to select. This value is ignored if nInd is set.

- nCrosses:

  total number of crosses to make

- nProgeny:

  number of progeny per cross

- trait:

  the trait for selection. Either a number indicating a single trait or
  a function returning a vector of length nInd.

- use:

  select on genetic values "gv", estimated breeding values "ebv",
  breeding values "bv", phenotypes "pheno", or randomly "rand"

- selectTop:

  selects highest values if true. Selects lowest values if false.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

- ...:

  additional arguments if using a function for trait

- balance:

  if using sexes, this option will balance the number of progeny per
  parent. This argument occurs after ..., so the argument name must be
  matched exactly.

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

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

#Select 4 individuals and make 8 crosses
pop2 = selectCross(pop, nInd=4, nCrosses=8, simParam=SP)
```
