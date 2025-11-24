# Edit genome - the top QTL

Edits the top QTL (with the largest additive effect) to a homozygous
state for the allele increasing. Only nonfixed QTL are edited The gv
slot is recalculated to reflect the any changes due to editing, but
other slots remain the same.

## Usage

``` r
editGenomeTopQtl(pop, ind, nQtl, trait = 1, increase = TRUE, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- ind:

  a vector of individuals to edit

- nQtl:

  number of QTL to edit

- trait:

  which trait effects should guide selection of the top QTL

- increase:

  should the trait value be increased or decreased

- simParam:

  an object of
  [`SimParam`](https://gaynorr.github.io/AlphaSimR/reference/SimParam.md)

## Value

Returns an object of
[`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)

#Create population
pop = newPop(founderPop, simParam=SP)

#Change up to 10 loci for individual 1
pop2 = editGenomeTopQtl(pop, ind=1, nQtl=10, simParam=SP)
```
