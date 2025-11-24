# Make random crosses

A wrapper for
[`makeCross2`](https://gaynorr.github.io/AlphaSimR/reference/makeCross2.md)
that randomly selects parental combinations for all possible
combinantions between two populations.

## Usage

``` r
randCross2(
  females,
  males,
  nCrosses,
  nProgeny = 1,
  balance = TRUE,
  femaleParents = NULL,
  maleParents = NULL,
  ignoreSexes = FALSE,
  simParam = NULL
)
```

## Arguments

- females:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  for female parents.

- males:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  for male parents.

- nCrosses:

  total number of crosses to make

- nProgeny:

  number of progeny per cross

- balance:

  this option will balance the number of progeny per parent

- femaleParents:

  an optional vector of indices for allowable female parents

- maleParents:

  an optional vector of indices for allowable male parents

- ignoreSexes:

  should sex be ignored

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
pop2 = randCross2(pop, pop, 10, simParam=SP)
```
