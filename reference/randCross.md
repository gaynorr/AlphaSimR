# Make random crosses

A wrapper for
[`makeCross`](https://gaynorr.github.io/AlphaSimR/reference/makeCross.md)
that randomly selects parental combinations for all possible
combinantions.

## Usage

``` r
randCross(
  pop,
  nCrosses,
  nProgeny = 1,
  balance = TRUE,
  parents = NULL,
  ignoreSexes = FALSE,
  simParam = NULL
)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- nCrosses:

  total number of crosses to make

- nProgeny:

  number of progeny per cross

- balance:

  if using sexes, this option will balance the number of progeny per
  parent

- parents:

  an optional vector of indices for allowable parents

- ignoreSexes:

  should sexes be ignored

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create population
pop = newPop(founderPop, simParam=SP)

#Make 10 crosses
pop2 = randCross(pop, 10, simParam=SP)
```
