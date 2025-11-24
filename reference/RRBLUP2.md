# RR-BLUP Model 2

Fits an RR-BLUP model for genomic predictions. This implementation is
meant for situations where
[`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md) is
too slow. Note that RRBLUP2 is only faster in certain situations, see
details below. Most users should use
[`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md).

## Usage

``` r
RRBLUP2(
  pop,
  traits = 1,
  use = "pheno",
  snpChip = 1,
  useQtl = FALSE,
  maxIter = 10,
  Vu = NULL,
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
  of the traits returning a single value. Unlike
  [`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md),
  only univariate models are supported.

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

  marker effect variance. If value is NULL, a reasonable starting point
  is chosen automatically.

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

## Details

The RRBLUP2 function works best when the number of markers is not too
large. This is because it solves the RR-BLUP problem by setting up and
solving Henderson's mixed model equations. Solving these equations
involves a square matrix with dimensions equal to the number of fixed
effects plus the number of random effects (markers). Whereas the
[`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md)
function solves the RR-BLUP problem using the EMMA approach. This
approach involves a square matrix with dimensions equal to the number of
phenotypic records. This means that the RRBLUP2 function uses less
memory than RRBLUP when the number of markers is approximately equal to
or smaller than the number of phenotypic records.

The RRBLUP2 function is not recommend for cases where the variance
components are unknown. This is uses the EM algorithm to solve for
unknown variance components, which is generally considerably slower than
the EMMA approach of
[`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md). The
number of iterations for the EM algorithm is set by maxIter. The default
value is typically too small for convergence. When the algorithm fails
to converge a warning is displayed, but results are given for the last
iteration. These results may be "good enough". However we make no claim
to this effect, because we can not generalize to all possible use cases.

The RRBLUP2 function can quickly solve the mixed model equations without
estimating variance components. The variance components are set by
defining Vu and Ve. Estimation of components is suppressed by setting
useEM to false. This may be useful if the model is being retrained
multiple times during the simulation. You could run
[`RRBLUP`](https://gaynorr.github.io/AlphaSimR/reference/RRBLUP.md)
function the first time the model is trained, and then use the variance
components from this output for all future runs with the RRBLUP2
functions. Again, we can make no claim to the general robustness of this
approach.

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
ans = RRBLUP2(pop, simParam=SP)
pop = setEBV(pop, ans, simParam=SP)

#Evaluate accuracy
cor(gv(pop), ebv(pop))
#>        est_GV_Trait1
#> Trait1     0.3905707
```
