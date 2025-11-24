# Quick founder haplotype simulation

Rapidly simulates founder haplotypes by randomly sampling 0s and 1s.
This is equivalent to having all loci with allele frequency 0.5 and
being in linkage equilibrium.

## Usage

``` r
quickHaplo(nInd, nChr, segSites, genLen = 1, ploidy = 2L, inbred = FALSE)
```

## Arguments

- nInd:

  number of individuals to simulate

- nChr:

  number of chromosomes to simulate

- segSites:

  number of segregating sites per chromosome

- genLen:

  genetic length of chromosomes

- ploidy:

  ploidy level of organism

- inbred:

  should founder individuals be inbred

## Value

an object of
[`MapPop-class`](https://gaynorr.github.io/AlphaSimR/reference/MapPop-class.md)

## Examples

``` r
# Creates a populations of 10 outbred individuals
# Their genome consists of 1 chromosome and 100 segregating sites
founderPop = quickHaplo(nInd=10,nChr=1,segSites=100)
```
