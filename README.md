# AlphaSimR

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/AlphaSimR)](https://cran.r-project.org/package=AlphaSimR)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/AlphaSimR)](https://cran.r-project.org/package=AlphaSimR)
[![R-CMD-check](https://github.com/gaynorr/AlphaSimR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gaynorr/AlphaSimR/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/gaynorr/AlphaSimR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/gaynorr/AlphaSimR)
<!-- badges: end -->

**Stochastic simulation of plant and animal breeding programs in R.**

AlphaSimR simulates breeding programs down to the DNA sequence of every
individual. You write the breeding program as an R script, one step at a
time, so a simulation can be as simple or as detailed as the question
needs.

- **Founder genomes** simulated from a demographic history with the
  bundled coalescent simulator MaCS, generated quickly at random, or
  imported from real haplotypes and genetic maps.
- **Traits** with additive, dominance, epistatic and
  genotype-by-environment effects, any number of them, genetically
  correlated, in diploid and polyploid species.
- **Breeding operations**: crossing designs, selfing, doubled haploids,
  phenotyping, and selection on phenotypes, true values, estimated
  breeding values or selection indices.
- **Genomic and pedigree prediction** with built-in RR-BLUP models,
  including dominance and combining ability models for hybrid breeding,
  and pedigree BLUP.
- **Bookkeeping** of pedigrees, recombination events and true genetic
  parameters, so results can be checked against the truth.

## Installation

Install the released version from CRAN:

``` r
install.packages("AlphaSimR")
```

Or the development version, which may be unstable, from the `devel`
branch on GitHub:

``` r
# install.packages("remotes")
remotes::install_github("gaynorr/AlphaSimR@devel")
```

## Example

Twenty generations of truncation selection on a trait controlled by
10,000 QTL:

``` r
library(AlphaSimR)

# Founder genomes and simulation parameters
founderPop = quickHaplo(nInd = 1000, nChr = 10, segSites = 1000)
SP = SimParam$new(founderPop)
SP$addTraitA(nQtlPerChr = 1000)
SP$setSexes("yes_sys")

# The breeding program
pop = newPop(founderPop)
genMean = meanG(pop)
for (generation in 1:20) {
  pop = selectCross(pop, nFemale = 500, nMale = 25, use = "gv",
                    nCrosses = 1000)
  genMean = c(genMean, meanG(pop))
}

plot(0:20, genMean, type = "l", xlab = "Generation",
     ylab = "Mean genetic value")
```

[Introduction to AlphaSimR](https://gaynorr.github.io/AlphaSimR/articles/AlphaSimR.html)
explains each step.

## Learn more

Start with the introduction, then read whichever topics your simulation
needs.

| Article | What it covers |
|:--|:--|
| [Introduction to AlphaSimR](https://gaynorr.github.io/AlphaSimR/articles/AlphaSimR.html) | A first simulation, step by step |
| [The Genome](https://gaynorr.github.io/AlphaSimR/articles/Genome.html) | Founder genomes, genetic maps and SNP chips |
| [Traits](https://gaynorr.github.io/AlphaSimR/articles/traits.html) | Specifying genetic architecture and phenotypes |
| [Simulating a Plant Breeding Program](https://gaynorr.github.io/AlphaSimR/articles/SimulatePlantBreedingProgram.html) | A complete multi-stage clonal breeding program |
| [Selection](https://gaynorr.github.io/AlphaSimR/articles/Selection.html) | Selection criteria, indices and response to selection |
| [Genomic Selection](https://gaynorr.github.io/AlphaSimR/articles/GenomicSelection.html) | Training and using genomic prediction models |
| [Quantitative Genetics](https://gaynorr.github.io/AlphaSimR/articles/QuanGen.html) | The theory behind the genetic values and variances reported |

The [full list of articles](https://gaynorr.github.io/AlphaSimR/articles/)
also covers importing data, pedigrees, recombination, combining
ability, optimal contribution selection and parallel computing. A free
course on breeding program simulation with AlphaSimR is available on
[edX](https://www.edx.org/learn/animal-breeding/the-university-of-edinburgh-breeding-programme-modelling-with-alphasimr).

## Getting help

- Questions about using AlphaSimR: [GitHub Discussions](https://github.com/gaynorr/AlphaSimR/discussions)
- Bug reports and feature requests: [GitHub Issues](https://github.com/gaynorr/AlphaSimR/issues)

Contributions are welcome; see the
[contributing guide](https://github.com/gaynorr/AlphaSimR/blob/devel/.github/CONTRIBUTING.md).
Please note that AlphaSimR is released with a
[Code of Conduct](https://github.com/gaynorr/AlphaSimR/blob/devel/.github/CODE_OF_CONDUCT.md).
By contributing to this project you agree to abide by its terms.

## Citation

If you use AlphaSimR, please cite:

Gaynor, R. C., Gorjanc, G., and Hickey, J. M. (2021). AlphaSimR: an R
package for breeding program simulations. *G3 Genes|Genomes|Genetics*,
11(2), jkaa017. <https://doi.org/10.1093/g3journal/jkaa017>

For plant breeding programs, see also: Bančič, J., Greenspoon, P.,
Gaynor, R. C., and Gorjanc, G. (2025). Plant breeding simulations with
AlphaSimR. *Crop Science*, 65(1), e21312.
<https://doi.org/10.1002/csc2.21312>

`citation("AlphaSimR")` gives both in BibTeX form. AlphaSimR is the
successor to AlphaSim
([Faux et al., 2016](https://doi.org/10.3835/plantgenome2016.02.0013)).
