# Fast RR-BLUP

Solves an RR-BLUP model for genomic predictions given known variance
components. This implementation is meant as a fast and low memory
alternative to
[`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md) or
[`RRBLUP2`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP2.md).
Unlike the those functions, the fastRRBLUP does not fit fixed effects
(other than the intercept) or account for unequal replication.

## Usage

``` r
fastRRBLUP(
  pop,
  traits = 1,
  use = "pheno",
  snpChip = 1,
  useQtl = FALSE,
  maxIter = 1000,
  Vu = NULL,
  Ve = NULL,
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
  of the traits returning a single value. Only univariate models are
  supported.

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

  maximum number of iterations.

- Vu:

  marker effect variance. If value is NULL, a reasonable value is chosen
  automatically.

- Ve:

  error variance. If value is NULL, a reasonable value is chosen
  automatically.

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
ans = fastRRBLUP(pop, simParam=SP)
pop = setEBV(pop, ans, simParam=SP)

#Evaluate accuracy
cor(gv(pop), ebv(pop))
#>        est_GV_Trait1
#> Trait1     0.7636568
```
