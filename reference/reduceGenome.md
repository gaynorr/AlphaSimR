# Create individuals with reduced ploidy

Creates new individuals from gametes. This function was created to model
the creation of diploid potatoes from tetraploid potatoes. It can be
used on any population with an even ploidy level. The newly created
individuals will have half the ploidy level of the originals. The
reduction can occur with or without genetic recombination.

## Usage

``` r
reduceGenome(
  pop,
  nProgeny = 1,
  useFemale = TRUE,
  keepParents = TRUE,
  simRecomb = TRUE,
  simParam = NULL
)
```

## Arguments

- pop:

  an object of 'Pop' superclass

- nProgeny:

  total number of progeny per individual

- useFemale:

  should female recombination rates be used.

- keepParents:

  should previous parents be used for mother and father.

- simRecomb:

  should genetic recombination be modeled.

- simParam:

  an object of 'SimParam' class

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create population
pop = newPop(founderPop, simParam=SP)

#Create individuals with reduced ploidy
pop2 = reduceGenome(pop, simParam=SP)
```
