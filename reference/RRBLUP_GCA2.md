# RR-BLUP GCA Model 2

Fits an RR-BLUP model that estimates seperate marker effects for females
and males. This implementation is meant for situations where
[`RRBLUP_GCA`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP_GCA.md)
is too slow. Note that RRBLUP_GCA2 is only faster in certain situations.
Most users should use
[`RRBLUP_GCA`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP_GCA.md).

## Usage

``` r
RRBLUP_GCA2(
  pop,
  traits = 1,
  use = "pheno",
  snpChip = 1,
  useQtl = FALSE,
  maxIter = 10,
  VuF = NULL,
  VuM = NULL,
  Ve = NULL,
  useEM = TRUE,
  tol = 1e-06,
  simParam = NULL,
  ...
)
```

## Arguments

- pop:

  a
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  to serve as the training population

- traits:

  an integer indicating the trait to model, a trait name, or a function
  of the traits returning a single value.

- use:

  train model using phenotypes "pheno", genetic values "gv", estimated
  breeding values "ebv", breeding values "bv", or randomly "rand"

- snpChip:

  an integer indicating which SNP chip genotype to use

- useQtl:

  should QTL genotypes be used instead of a SNP chip. If TRUE, snpChip
  specifies which trait's QTL to use, and thus these QTL may not match
  the QTL underlying the phenotype supplied in traits.

- maxIter:

  maximum number of iterations for convergence.

- VuF:

  marker effect variance for females. If value is NULL, a reasonable
  starting point is chosen automatically.

- VuM:

  marker effect variance for males. If value is NULL, a reasonable
  starting point is chosen automatically.

- Ve:

  error variance. If value is NULL, a reasonable starting point is
  chosen automatically.

- useEM:

  use EM to solve variance components. If false, the initial values are
  considered true.

- tol:

  tolerance for EM algorithm convergence

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

- ...:

  additional arguments if using a function for traits

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
ans = RRBLUP_GCA2(pop, simParam=SP)
pop = setEBV(pop, ans, simParam=SP)

#Evaluate accuracy
cor(gv(pop), ebv(pop))
#>        est_GV_Trait1
#> Trait1     0.7304646
```
