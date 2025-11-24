# Variance of estimated breeding values

Returns variance of estimated breeding values for all traits

## Usage

``` r
varEBV(pop)
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
varA(pop)
#> Error in get("SP", envir = .GlobalEnv): object 'SP' not found
varEBV(pop)
#>           Trait1
#> Trait1 0.2847924
```
