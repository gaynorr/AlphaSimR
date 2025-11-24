# Set GCA as phenotype

Calculates general combining ability from a set of testers and returns
these values as phenotypes for a population.

## Usage

``` r
setPhenoGCA(
  pop,
  testers,
  use = "pheno",
  h2 = NULL,
  H2 = NULL,
  varE = NULL,
  corE = NULL,
  reps = 1,
  fixEff = 1L,
  p = NULL,
  inbred = FALSE,
  onlyPheno = FALSE,
  simParam = NULL
)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- testers:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- use:

  true genetic value (`gv`) or phenotypes (`pheno`, default)

- h2:

  a vector of desired narrow-sense heritabilities for each trait. See
  details in
  [`setPheno`](https://gaynorr.github.io/AlphaSimR/reference/setPheno.md).

- H2:

  a vector of desired broad-sense heritabilities for each trait. See
  details in
  [`setPheno`](https://gaynorr.github.io/AlphaSimR/reference/setPheno.md).

- varE:

  error (co)variances for traits. See details in
  [`setPheno`](https://gaynorr.github.io/AlphaSimR/reference/setPheno.md).

- corE:

  an optional matrix for correlations between errors. See details in
  [`setPheno`](https://gaynorr.github.io/AlphaSimR/reference/setPheno.md).

- reps:

  number of replications for phenotype. See details in
  [`setPheno`](https://gaynorr.github.io/AlphaSimR/reference/setPheno.md).

- fixEff:

  fixed effect to assign to the population. Used by genomic selection
  models only.

- p:

  the p-value for the environmental covariate used by GxE traits. If
  NULL, a value is sampled at random.

- inbred:

  are both pop and testers fully inbred. They are only fully inbred if
  created by
  [`newPop`](https://gaynorr.github.io/AlphaSimR/reference/newPop.md)
  using inbred founders or by the
  [`makeDH`](https://gaynorr.github.io/AlphaSimR/reference/makeDH.md)
  function

- onlyPheno:

  should only the phenotype be returned, see return

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
or a matrix if onlyPheno=TRUE

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10, inbred=TRUE)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)

#Set phenotype to average per
pop2 = setPhenoGCA(pop, pop, use="gv", inbred=TRUE, simParam=SP)
```
