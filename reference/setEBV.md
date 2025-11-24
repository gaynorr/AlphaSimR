# Set estimated breeding values (EBV)

Adds genomic estimated values to a populations's EBV slot using output
from a genomic selection functions. The genomic estimated values can be
either estimated breeding values, estimated genetic values, or estimated
general combining values.

## Usage

``` r
setEBV(
  pop,
  solution,
  value = "gv",
  targetPop = NULL,
  append = FALSE,
  simParam = NULL
)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- solution:

  an object of
  [`RRsol-class`](https://gaynorr.github.io/AlphaSimR/reference/RRsol-class.md)

- value:

  the genomic value to be estimated. Can be either "gv", "bv", "female",
  or "male".

- targetPop:

  an optional target population that can be used when value is "bv",
  "female", or "male". When supplied, the allele frequency in the
  targetPop is used to set these values.

- append:

  should estimated values be appended to existing data in the EBV slot.
  If TRUE, a new column is added. If FALSE, existing data is replaced
  with the new estimates.

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)
SP$addSnpChip(10)

#Create population
pop = newPop(founderPop, simParam=SP)

#Run GS model and set EBV
ans = RRBLUP(pop, simParam=SP)
pop = setEBV(pop, ans, simParam=SP)

#Evaluate accuracy
cor(gv(pop), ebv(pop))
#>        est_GV_Trait1
#> Trait1     0.7102514
```
