# Phenotypic variance

Returns phenotypic variance for all traits

## Usage

``` r
varP(pop)
```

## Arguments

- pop:

  an object of
  [`Pop-class`](https://gaynorr.github.io/AlphaSimR/reference/Pop-class.md)
  or
  [`HybridPop-class`](https://gaynorr.github.io/AlphaSimR/reference/HybridPop-class.md)

## Examples

``` r
#Create founder haplotypes
founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)

#Set simulation parameters
SP = SimParam$new(founderPop)
SP$addTraitA(10)
SP$setVarE(h2=0.5)

#Create population
pop = newPop(founderPop, simParam=SP)
varP(pop)
#>          Trait1
#> Trait1 2.492171
```
