# Sample epistatic effects

Samples epistatic effects from a normal distribution or gamma
distribution with a variance relative to the variance of the additive
effects

## Usage

``` r
sampEpiEff(qtlLoci, nTraits, corr, gamma, shape, relVar)
```

## Arguments

- qtlLoci:

  total number of loci

- nTraits:

  number of traits

- corr:

  correlation between epistatic effects

- gamma:

  indicator of trait should use a gamma distribution

- shape:

  gamma distribution shape parameter

- relVar:

  desired variance for epistatic effects

## Value

a matrix with dimensions qtlLoci by nTraits
