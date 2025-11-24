# Calculate Mendelian sampling

Calculate Mendelian sampling

## Usage

``` r
mendelianSampling(
  pop,
  parents = NULL,
  mothers = NULL,
  fathers = NULL,
  use = "gv",
  simParam = NULL
)
```

## Arguments

- pop:

  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  with individuals whose parent average will be calculated

- parents:

  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  with mothers and fathers of individuals in `pop`; if `NULL` must
  provide `mothers` and `fathers`

- mothers:

  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  with mothers of individuals in `pop`; if `NULL` must provide `parents`

- fathers:

  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  with fathers of individuals in `pop`; if `NULL` must provide `parents`

- use:

  character, calculate using
  `"`[`gv`](https://gaynorr.github.io/AlphaSimR/reference/gv.md)`"`,
  `"`[`bv`](https://gaynorr.github.io/AlphaSimR/reference/bv.md)`"`,
  `"`[`ebv`](https://gaynorr.github.io/AlphaSimR/reference/ebv.md)`"`,
  or
  `"`[`pheno`](https://gaynorr.github.io/AlphaSimR/reference/pheno.md)`"`

- simParam:

  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)
  object

## Value

a matrix of Mendelian samplings with dimensions nInd by nTraits

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitAD(10, meanDD=0.5)
SP$setVarE(h2=0.5)

#Create population
pop = newPop(founderPop, simParam=SP)
pop2 = randCross(pop, nCrosses=10, nProgeny=2)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found
mendelianSampling(pop2, parents = pop)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found
mendelianSampling(pop2, mothers = pop, fathers = pop)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found
```
