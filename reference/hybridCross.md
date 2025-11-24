# Hybrid crossing

A convenient function for hybrid plant breeding simulations. Allows for
easy specification of a test cross scheme and/or creation of an object
of
[`HybridPop-class`](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.md).
Note that the
[`HybridPop-class`](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.md)
should only be used if the parents were created using the
[`makeDH`](https://gaynorr.github.io/AlphaSimR/reference/makeDH.md)
function or
[`newPop`](https://gaynorr.github.io/AlphaSimR/reference/newPop.md)
using inbred founders. The id for new individuals is
\[mother_id\]\_\[father_id\]

## Usage

``` r
hybridCross(
  females,
  males,
  crossPlan = "testcross",
  returnHybridPop = FALSE,
  simParam = NULL
)
```

## Arguments

- females:

  female population, an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- males:

  male population, an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- crossPlan:

  either "testcross" for all possible combinations or a matrix with two
  columns for designed crosses

- returnHybridPop:

  should results be returned as
  [`HybridPop-class`](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.md).
  If false returns results as
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md).
  Population must be fully inbred if TRUE.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)

#Create population
pop = newPop(founderPop, simParam=SP)

#Make crosses for full diallele
pop2 = hybridCross(pop, pop, simParam=SP)
```
