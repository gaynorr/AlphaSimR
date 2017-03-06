#' @useDynLib AlphaSimR
#' @import Rcpp RcppArmadillo

#' @title Simulation parameters
#' 
#' @description Container for global simulation parameters
#' 
#' @slot ploidy ploidy level of species
#' @slot nChr number of chromosomes
#' @slot nTraits number of traits
#' @slot nSnpChips number of SNP chips
#' @slot segSites segregating sites per chromosome
#' @slot genMaps list of chromsome genetic maps
#' @slot traits list of trait
#' @slot snpChips list of SNP chips
#'
#' @export
setClass("SimParams",
         slots=c(ploidy="integer",
                 nChr="integer",
                 nTraits="integer",
                 nSnpChips="integer",
                 segSites="integer",
                 genMaps="list",
                 traits="list",
                 snpChips="list"))


setValidity("SimParams",function(object){
  errors = character()
  if(object@nChr!=length(object@segSites)){
    errors = c(errors,"nChr!=length(segSites)")
  }
  if(object@nChr!=length(object@genMap)){
    errors = c(errors,"nChr!=length(genMaps)")
  }
  if(object@nTraits!=length(object@traits)){
    errors = c(errors,"nTraits!=length(traits)")
  }
  if(object@nSnpChips!=length(object@snpChips)){
    errors = c(errors,"nSnpChips!=length(snpChips)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})
