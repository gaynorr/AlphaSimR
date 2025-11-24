# Edit genome

Edits selected loci of selected individuals to a homozygous state for
either the 1 or 0 allele. The gv slot is recalculated to reflect the any
changes due to editing, but other slots remain the same.

## Usage

``` r
editGenome(pop, ind, chr, segSites, allele, simParam = NULL)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)

- ind:

  a vector of individuals to edit

- chr:

  a vector of chromosomes to edit. Length must match length of segSites.

- segSites:

  a vector of segregating sites to edit. Length must match length of
  chr.

- allele:

  either 0 or 1 for desired allele

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

#Change individual 1 to homozygous for the 1 allele
#at locus 1, chromosome 1
pop2 = editGenome(pop, ind=1, chr=1, segSites=1,
                  allele=1, simParam=SP)
```
