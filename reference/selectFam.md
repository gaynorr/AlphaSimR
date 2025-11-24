# Select families

Selects a subset of full-sib families from a population.

## Usage

``` r
selectFam(
  pop,
  nFam,
  trait = 1,
  use = "pheno",
  sex = "B",
  famType = "B",
  selectTop = TRUE,
  returnPop = TRUE,
  candidates = NULL,
  simParam = NULL,
  ...
)
```

## Arguments

- pop:

  and object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md),
  [`HybridPop-class`](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.md)
  or
  [`MultiPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MultiPop-class.md)

- nFam:

  the number of families to select

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

- sex:

  which sex to select. Use "B" for both, "F" for females and "M" for
  males. If the simulation is not using sexes, the argument is ignored.

- famType:

  which type of family to select. Use "B" for full-sib families, "F" for
  half-sib families on female side and "M" for half-sib families on the
  male side.

- selectTop:

  selects highest values if true. Selects lowest values if false.

- returnPop:

  should results be returned as a
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md).
  If FALSE, only the index of selected individuals is returned.

- candidates:

  an optional vector of eligible selection candidates.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

- ...:

  additional arguments if using a function for `trait` and `use`

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md),
[`HybridPop-class`](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.md)
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

#Create 3 biparental families with 10 progeny
pop2 = randCross(pop, nCrosses=3, nProgeny=10, simParam=SP)

#Select best 2 families
pop3 = selectFam(pop2, 2, simParam=SP)
```
