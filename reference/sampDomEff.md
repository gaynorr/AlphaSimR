# Sample dominance effects

Samples dominance deviation effects from a normal distribution and uses
previously sampled additive effects to form dominance effects

## Usage

``` r
sampDomEff(qtlLoci, nTraits, addEff, corDD, meanDD, varDD)
```

## Arguments

- qtlLoci:

  total number of loci

- nTraits:

  number of traits

- addEff:

  previously sampled additive effects

- corDD:

  correlation between dominance degrees

- meanDD:

  mean value of dominance degrees

- varDD:

  variance of dominance degrees

## Value

a matrix with dimensions qtlLoci by nTraits
