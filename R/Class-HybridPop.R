#HybridPop----
#' @title Population
#' 
#' @description Population class
#' 
#' @slot nInd number of individuals
#' @slot nTraits number of traits
#' @slot gv matrix of genotypic values
#' @slot pheno matrix of phenotypic values
#' @slot id hybrid's id
#' @slot par1 hybrid's first parent/mother
#' @slot par2 hybrid's second parent/father
#' 
#' @export
setClass("HybridPop",
         slots=c(nInd="integer",
                 nTraits="integer",
                 gv="matrix",
                 pheno="matrix",
                 id="character",
                 par1="character",
                 par2="character"))

setValidity("HybridPop",function(object){
  errors = character()
  if(object@nInd!=nrow(object@gv)){
    errors = c(errors,"nInd!=nrow(gv)")
  }
  if(object@nInd!=nrow(object@pheno)){
    errors = c(errors,"nInd!=nrow(pheno)")
  }
  if(ncol(object@gv)!=object@nTraits){
    errors = c(errors,"ncol(gv)!=nTraits")
  }
  if(ncol(object@pheno)!=object@nTraits){
    errors = c(errors,"ncol(pheno)!=nTraits")
  }
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

setMethod("[",
          signature(x = "HybridPop"),
          function(x, i, j=NULL, ..., drop = TRUE){
            x@id = x@id[i]
            x@par1 = x@par1[i]
            x@par2 = x@par2[i]
            x@gv = matrix(x@gv[i,],ncol=x@nTraits)
            x@pheno = matrix(x@pheno[i,],ncol=x@nTraits)
            x@nInd = length(x@id)
            validObject(x)
            return(x)
          }
)
