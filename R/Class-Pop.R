#Pop----
#' @title Population
#'
#' @slot nInd number of individuals
#' @slot geno list of chromosome genotypes
#'
#' @return
#' @export
#'
#' @examples
setClass("Pop",
         slots=c(nInd="integer",
                 geno="list"))

#TraitPop----
#' @title Traited Population
#' @description Extends \code{\link{Pop-class}}
#' @slot gv matrix of genetic values
#'
#' @return
#' @export
#'
#' @examples
setClass("TraitPop",
         slots=c(gv="matrix"),
         contains="Pop")

setValidity("TraitPop",function(object){
  errors = character()
  if(object@nInd!=nrow(object@gv)){
    errors = c(errors,"nInd!=nrow(gv)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#PhenoPop----
#' @title Phenotyped Population
#' @description Extends \code{\link{TraitPop-class}}
#' @slot pheno matrix of phenotypes
#'
#' @return
#' @export
#'
#' @examples
setClass("PhenoPop",
         slots=c(pheno="matrix"),
         contains="TraitPop")

setValidity("PhenoPop",function(object){
  errors = character()
  if(object@nInd!=nrow(object@pheno)){
    errors = c(errors,"nInd!=nrow(pheno)")
  }
  if(ncol(object@gv)!=ncol(object@pheno)){
    errors = c(errors,"ncol(gv)!=ncol(pheno)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#PedPop----
#' @title Pedigreed Population
#' @description Extends \code{\link{PhenoPop-class}}
#' @slot id individual's identifier
#' @slot par1 individual's female parent
#' @slot par2 individual's male parent
#'
#' @return
#' @export
#'
#' @examples
setClass("PedPop",
         slots=c(id="character",
                 par1="character",
                 par2="character"),
         contains="TraitPop")

setValidity("PedPop",function(object){
  errors = character()
  if(object@nInd!=length(object@id)){
    errors = c(errors,"nInd!=length(id)")
  }
  if(object@nInd!=length(object@par1)){
    errors = c(errors,"nInd!=length(par1)")
  }
  if(object@nInd!=length(object@par2)){
    errors = c(errors,"nInd!=length(par2)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})
