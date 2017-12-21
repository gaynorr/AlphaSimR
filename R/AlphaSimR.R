#' @useDynLib AlphaSimR, .registration = TRUE
#' @import Rcpp RcppArmadillo
#' @importFrom methods new validObject
#' @importFrom stats aggregate rnorm qnorm var
#' @importFrom stats coef dnorm lm pnorm qgamma
#' @importFrom stats model.matrix
#' @importFrom utils combn read.table write.table
#' @importFrom R6 R6Class

#' @title AlphaSimR: Breeding Program Simulations
#'
#' @description
#' This package contains classes and functions for 
#' simulating plant and animal breeding programs.
#'
#' @docType package
#' @name AlphaSimR-package
NULL