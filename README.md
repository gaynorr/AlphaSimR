# README #

CRAN stats: [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/AlphaSimR)](https://cran.r-project.org/package=AlphaSimR)
[![](http://cranlogs.r-pkg.org/badges/grand-total/AlphaSimR)](https://cran.r-project.org/package=AlphaSimR)
[![](http://cranlogs.r-pkg.org/badges/AlphaSimR)](https://cran.r-project.org/package=AlphaSimR)

R CMD checks: [![CRAN](https://cranchecks.info/badges/summary/AlphaSimR?label=CRAN)](https://cran.r-project.org/web/checks/check_results_AlphaSimR.html)
[![R universe](https://gaynorr.r-universe.dev/AlphaSimR/badges/checks?label=R-universe)](https://gaynorr.r-universe.dev/AlphaSimR)
[![GitHub](https://img.shields.io/github/actions/workflow/status/gaynorr/AlphaSimR/R-CMD-check.yaml?label=GitHub)](https://github.com/gaynorr/AlphaSimR/actions/workflows/R-CMD-check.yaml)

The successor to the 'AlphaSim' software for breeding program simulation (Faux et al., 2016; https://doi.org/10.3835/plantgenome2016.02.0013). Used for stochastic simulations of breeding programs to the level of DNA sequence for every individual. Contained is a wide range of functions for modeling common tasks in a breeding program, such as selection and crossing. These functions allow for constructing simulations of highly complex plant and animal breeding programs via scripting in the R software environment. Such simulations can be used to evaluate overall breeding program performance and conduct research into breeding program design, such as implementation of genomic selection. Included is the 'Markovian Coalescent Simulator' ('MaCS') for fast simulation of biallelic sequences according to a population demographic history (Chen et al., 2009; https://doi.org/10.1101/gr.083634.108).

## Publications

Gaynor, R Chris and Gorjanc, Gregor and Hickey, John M (2021) AlphaSimR: an R package for breeding program simulations. G3 Gene|Genomes|Genetics, 11(2):jkaa017. https://doi.org/10.1093/g3journal/jkaa017.

Bančič, Jon and Greenspoon, Philip and Gaynor, R Chris and Gorjanc, Gregor (2024) Plant breeding simulations with AlphaSimR. Crop Science, 65(1):e21312. https://doi.org/10.1002/csc2.21312.

## Online course

Free course on breeding program simulations and how to use AlphaSimR is available on [edX](https://www.edx.org/learn/animal-breeding/the-university-of-edinburgh-breeding-programme-modelling-with-alphasimr).

## Download

[AlphaSimR](https://cran.r-project.org/package=AlphaSimR) is available on CRAN.

To install use:

    install.packages('AlphaSimR')

The development version of AlphaSimR (potentially unstable) can be accessed from the devel branch on GitHub.

To install use:

    devtools::install_github(repo="gaynorr/AlphaSimR@devel")

To install with vignettes use:

    devtools::install_github(repo="gaynorr/AlphaSimR@devel", build_vignettes=TRUE)
