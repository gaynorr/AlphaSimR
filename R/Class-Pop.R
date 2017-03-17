# The Pop superclass contains genotypes
# and summary data for multiple individuals


#Pop----
#' @title Population
#' 
#' @description Population class
#' 
#' @slot nInd number of individuals
#' @slot gender gender of individuals
#' @slot geno list containing chromosome genotypes
#' 
#' @export
setClass("Pop",
         slots=c(nInd="integer",
                 gender="character",
                 geno="list"))
setValidity("Pop",function(object){
  errors = character()
  if(object@nInd!=length(object@gender)){
    errors = c(errors,"nInd!=length(gender)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#TraitPop----
#' @title Population with traits
#' 
#' @description Extends \code{\link{Pop-class}}
#' 
#' @slot gv matrix of genetic values
#' @slot pheno matrix of phenotypic values
#'
#' @export
setClass("TraitPop",
         slots=c(gv="matrix",
                 pheno="matrix"),
         contains="Pop")

setValidity("TraitPop",function(object){
  errors = character()
  if(object@nInd!=nrow(object@gv)){
    errors = c(errors,"nInd!=nrow(gv)")
  }
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
#' @title Population with pedigree
#' 
#' @description Extends \code{\link{TraitPop-class}}
#' 
#' @slot id individual's identifier
#' @slot par1 individual's female parent
#' @slot par2 individual's male parent
#'
#' @export
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
