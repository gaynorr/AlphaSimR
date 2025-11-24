# Generates DH lines

Creates DH lines from each individual in a population. Only works with
diploid individuals. For polyploids, use
[`reduceGenome`](https://gaynorr.github.io/AlphaSimR/reference/reduceGenome.md)
and
[`doubleGenome`](https://gaynorr.github.io/AlphaSimR/reference/doubleGenome.md).

## Usage

``` r
makeDH(pop, nDH = 1, useFemale = TRUE, keepParents = TRUE, simParam = NULL)
```

## Arguments

- pop:

  an object of 'Pop' superclass

- nDH:

  total number of DH lines per individual

- useFemale:

  should female recombination rates be used.

- keepParents:

  should previous parents be used for mother and father.

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

#Create 1 DH for each individual
pop2 = makeDH(pop, simParam=SP)
```
