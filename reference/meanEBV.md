# Mean estimated breeding values

Returns the mean estimated breeding values for all traits

## Usage

``` r
meanEBV(pop)
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
trtH2 = 0.5
SP$setVarE(h2=trtH2)

#Create population
pop = newPop(founderPop, simParam=SP)
pop@ebv = trtH2 * (pop@pheno - meanP(pop)) #ind performance based EBV
meanEBV(pop)
#>        Trait1 
#> -8.326673e-18 
```
