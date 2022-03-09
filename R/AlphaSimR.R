#' @useDynLib AlphaSimR, .registration = TRUE
#' @import Rcpp
#' @importFrom methods new validObject is .hasSlot 
#' @importFrom methods show classLabel
#' @importFrom stats aggregate rnorm qnorm var
#' @importFrom stats coef dnorm lm pnorm qgamma na.omit
#' @importFrom stats model.matrix rbinom runif cov2cor
#' @importFrom utils combn read.table write.table
#' @importFrom R6 R6Class

#' @description
#' The successor to the 'AlphaSim' software for breeding program 
#' simulation [Faux et al. (2016) <doi:10.3835/plantgenome2016.02.0013>]. 
#' Used for stochastic simulations of breeding programs to the level of DNA 
#' sequence for every individual. Contained is a wide range of functions for 
#' modeling common tasks in a breeding program, such as selection and crossing. 
#' These functions allow for constructing simulations of highly complex plant and 
#' animal breeding programs via scripting in the R software environment. Such 
#' simulations can be used to evaluate overall breeding program performance and 
#' conduct research into breeding program design, such as implementation of 
#' genomic selection. Included is the 'Markovian Coalescent Simulator' ('MaCS') 
#' for fast simulation of biallelic sequences according to a population 
#' demographic history [Chen et al. (2009) <doi:10.1101/gr.083634.108>].
#' 
#' Please see the introductory vignette for instructions for using this package. 
#' The vignette can be viewed using the following command: 
#' \code{vignette("intro",package="AlphaSimR")}
#' @keywords internal
"_PACKAGE"