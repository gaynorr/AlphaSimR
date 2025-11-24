# Pedigree cross

Creates a
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
from a generic pedigree and a set of founder individuals.

The way in which the user supplied pedigree is used depends on the value
of matchID. If matchID is TRUE, the IDs in the user supplied pedigree
are matched against founderNames. If matchID is FALSE, founder
individuals in the user supplied pedigree are randomly sampled from
founderPop.

## Usage

``` r
pedigreeCross(
  founderPop,
  id,
  mother,
  father,
  matchID = FALSE,
  maxCycle = 100,
  DH = NULL,
  nSelf = NULL,
  useFemale = TRUE,
  simParam = NULL
)
```

## Arguments

- founderPop:

  a
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- id:

  a vector of unique identifiers for individuals in the pedigree. The
  values of these IDs are seperate from the IDs in the founderPop if
  matchID=FALSE.

- mother:

  a vector of identifiers for the mothers of individuals in the
  pedigree. Must match one of the elements in the id vector or they will
  be treated as unknown.

- father:

  a vector of identifiers for the fathers of individuals in the
  pedigree. Must match one of the elements in the id vector or they will
  be treated as unknown.

- matchID:

  indicates if the IDs in founderPop should be matched to the id
  argument. See details.

- maxCycle:

  the maximum number of loops to make over the pedigree to sort it.

- DH:

  an optional vector indicating if an individual should be made a
  doubled haploid.

- nSelf:

  an optional vector indicating how many generations an individual
  should be selfed.

- useFemale:

  If creating DH lines, should female recombination rates be used. This
  parameter has no effect if, recombRatio=1.

- simParam:

  an object of 'SimParam' class

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create population
pop = newPop(founderPop, simParam=SP)

#Pedigree for a biparental cross with 7 generations of selfing
id = 1:10
mother = c(0,0,1,3:9)
father = c(0,0,2,3:9)
pop2 = pedigreeCross(pop, id, mother, father, simParam=SP)
```
