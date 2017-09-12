#' @title Simulation parameters
#' 
#' @description 
#' Container for global simulation parameters. Saving this object 
#' as SIMPARAM will allow it to be accessed by function defaults.
#' 
#' @slot ploidy ploidy level of species
#' @slot nChr number of chromosomes
#' @slot nTraits number of traits
#' @slot nSnpChips number of SNP chips
#' @slot segSites segregating sites per chromosome
#' @slot gender is gender used for mating
#' @slot genMaps "matrix" of chromsome genetic maps
#' @slot recombRatio ratio of genetic recombination in 
#' females relative to male
#' @slot traits list of trait
#' @slot snpChips list of SNP chips
#' @slot potQtl list of potential QTL segregating sites
#' @slot potSnp list of potential SNP segregating sites
#' @slot lastId last ID number assigned
#'
#' @export
setClass("SimParam",
         slots=c(ploidy="integer",
                 nChr="integer",
                 nTraits="integer",
                 nSnpChips="integer",
                 segSites="integer",
                 gender="character",
                 genMaps="matrix",
                 recombRatio="numeric",
                 traits="list",
                 snpChips="list",
                 potQtl="list",
                 potSnp="list",
                 lastId="numeric"))


setValidity("SimParam",function(object){
  errors = character()
  if(object@nChr!=length(object@segSites)){
    errors = c(errors,"nChr!=length(segSites)")
  }
  if(object@nChr!=length(object@genMaps)){
    errors = c(errors,"nChr!=length(genMaps)")
  }
  if(object@nTraits!=length(object@traits)){
    errors = c(errors,"nTraits!=length(traits)")
  }
  if(object@nSnpChips!=length(object@snpChips)){
    errors = c(errors,"nSnpChips!=length(snpChips)")
  }
  if(object@nChr!=length(object@potQtl)){
    errors = c(errors,"nChr!=length(potQtl)")
  }
  if(object@nChr!=length(object@potSnp)){
    errors = c(errors,"nChr!=length(potSnp)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})
