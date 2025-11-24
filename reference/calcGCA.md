# Calculate GCA

Calculate general combining ability of test crosses. Intended for output
from hybridCross using the "testcross" option, but will work for any
population.

## Usage

``` r
calcGCA(pop, use = "pheno")
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  or
  [`HybridPop-class`](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.md)

- use:

  tabulate either genetic values "gv", estimated breeding values "ebv",
  or phenotypes "pheno"

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10, inbred=TRUE)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)

#Make crosses for full diallele
pop2 = hybridCross(pop, pop, simParam=SP)
GCA = calcGCA(pop2, use="gv")
```
