# Add Random Mutations

Adds random mutations to individuals in a population. Note that any
existing phenotypes or EBVs are kept. Thus, the user will need to run
[`setPheno`](https://gaynorr.github.io/AlphaSimR/reference/setPheno.md)
and/or
[`setEBV`](https://gaynorr.github.io/AlphaSimR/reference/setEBV.md) to
generate new phenotypes or EBVs that reflect changes introduced by the
new mutations.

## Usage

``` r
mutate(pop, mutRate = 2.5e-08, returnPos = FALSE, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- mutRate:

  rate of new mutations

- returnPos:

  should the positions of mutations be returned

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
if returnPos=FALSE or a list containing a
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
and a data.frame containing the postions of mutations if returnPos=TRUE

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)

#Introduce mutations
pop = mutate(pop, simParam=SP)
```
