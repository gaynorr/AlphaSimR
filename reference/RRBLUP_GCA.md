# RR-BLUP GCA Model

Fits an RR-BLUP model that estimates seperate marker effects for females
and males. Useful for predicting GCA of parents in single cross hybrids.
Can also predict performance of specific single cross hybrids.

## Usage

``` r
RRBLUP_GCA(
  pop,
  traits = 1,
  use = "pheno",
  snpChip = 1,
  useQtl = FALSE,
  maxIter = 40L,
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
ans = RRBLUP_GCA(pop, simParam=SP)
pop = setEBV(pop, ans, simParam=SP)

#Evaluate accuracy
cor(gv(pop), ebv(pop))
#>        est_GV_Trait1
#> Trait1     0.8112236
```
