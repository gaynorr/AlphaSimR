# README #

[![R-CMD-check](https://github.com/gaynorr/AlphaSimR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gaynorr/AlphaSimR/actions/workflows/R-CMD-check.yaml)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/AlphaSimR)](https://cran.r-project.org/package=AlphaSimR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/AlphaSimR)](https://cran.r-project.org/package=AlphaSimR)
[![](http://cranlogs.r-pkg.org/badges/AlphaSimR)](https://cran.r-project.org/package=AlphaSimR)

The successor to the 'AlphaSim' software for breeding program simulation (Faux et al., 2016; https://doi.org/10.3835/plantgenome2016.02.0013). Used for stochastic simulations of breeding programs to the level of DNA sequence for every individual. Contained is a wide range of functions for modeling common tasks in a breeding program, such as selection and crossing. These functions allow for constructing simulations of highly complex plant and animal breeding programs via scripting in the R software environment. Such simulations can be used to evaluate overall breeding program performance and conduct research into breeding program design, such as implementation of genomic selection. Included is the 'Markovian Coalescent Simulator' ('MaCS') for fast simulation of biallelic sequences according to a population demographic history (Chen et al., 2009; https://doi.org/10.1101/gr.083634.108).

## Publication

Gaynor, R. Chris, Gregor Gorjanc, and John M. Hickey. 2021. AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

## Download

[AlphaSimR](https://cran.r-project.org/package=AlphaSimR) is available on CRAN.

To install use:

    install.packages('AlphaSimR')

The development version of AlphaSimR (potentially unstable) can be accessed from the devel branch on GitHub.

To install use:

    devtools::install_github(repo="gaynorr/AlphaSimR@devel")

To install with vignettes use:

    devtools::install_github(repo="gaynorr/AlphaSimR@devel", build_vignettes=TRUE)

