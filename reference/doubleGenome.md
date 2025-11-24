# Double the ploidy of individuals

Creates new individuals with twice the ploidy. This function was created
to model the formation of tetraploid potatoes from diploid potatoes.
This function will work on any population.

## Usage

``` r
doubleGenome(pop, keepParents = TRUE, simParam = NULL)
```

## Arguments

- pop:

  an object of 'Pop' superclass

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

#Create individuals with doubled ploidy
pop2 = doubleGenome(pop, simParam=SP)
```
