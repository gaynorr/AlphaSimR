#HybridPop----
#' @title Hybrid population
#' 
#' @description
#' A lightweight version of \code{\link{Pop-class}} for hybrid lines.
#' Memory is saved by not storing genotypic data.
#' 
#' @param x a 'HybridPop'
#' @param i index of individuals
#' @param ... additional 'HybridPop' objects 
#' 
#' @slot nInd number of individuals
#' @slot id an individual's identifier
#' @slot mother the identifier of the individual's mother
#' @slot father the identifier of the individual's father
#' @slot nTraits number of traits
#' @slot gv matrix of genetic values. When using GxE traits,
#' gv reflects gv when p=0.5. Dimensions are nInd by nTraits.
#' @slot pheno matrix of phenotypic values. Dimensions are
#' nInd by nTraits.
#' @slot gxe list containing GxE slopes for GxE traits
#' 
#' @export
setClass("HybridPop",
         slots=c(nInd="integer",
                 id="character",
                 mother="character",
                 father="character",
                 nTraits="integer",
                 gv="matrix",
                 pheno="matrix",
                 gxe="list"))

setValidity("HybridPop",function(object){
  errors = character()
  if(any(grepl(" ",object@id,fixed=TRUE))){
    errors = c(errors,"id can not contain spaces")
  }
  if(any(grepl(" ",object@mother,fixed=TRUE))){
    errors = c(errors,"mother can not contain spaces")
  }
  if(any(grepl(" ",object@father,fixed=TRUE))){
    errors = c(errors,"father can not contain spaces")
  }
  if(object@nInd!=length(object@id)){
    errors = c(errors,"nInd!=length(id)")
  }
  if(object@nInd!=length(object@mother)){
    errors = c(errors,"nInd!=length(mother)")
  }
  if(object@nInd!=length(object@father)){
    errors = c(errors,"nInd!=length(father)")
  }
  if(object@nInd!=nrow(object@gv)){
    errors = c(errors,"nInd!=nrow(gv)")
  }
  if(object@nInd!=nrow(object@pheno)){
    errors = c(errors,"nInd!=nrow(pheno)")
  }
  if(object@nTraits!=ncol(object@gv)){
    errors = c(errors,"nTraits!=ncol(gv)")
  }
  if(object@nTraits!=ncol(object@pheno)){
    errors = c(errors,"nTraits!=ncol(pheno)")
  }
  if(object@nTraits!=length(object@gxe)){
    errors = c(errors,"nTraits!=length(gxe)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#' @describeIn HybridPop Extract HybridPop using index or id
setMethod("[",
          signature(x = "HybridPop"),
          function(x, i){
            if(is.character(i)){
              i = x@id%in%i
            }
            x@id = x@id[i]
            x@mother = x@mother[i]
            x@father = x@father[i]
            x@gv = x@gv[i,,drop=FALSE]
            x@pheno = x@pheno[i,,drop=FALSE]
            x@nInd = length(x@id)
            if(x@nTraits>=1){
              for(trait in 1:x@nTraits){
                if(!is.null(x@gxe[[trait]])){
                  x@gxe[[trait]] = x@gxe[[trait]][i]
                }
              }
            }
            validObject(x)
            return(x)
          }
)

#' @describeIn HybridPop Combine multiple HybridPops
setMethod("c",
          signature(x = "HybridPop"),
          function (x, ...){
            for(y in list(...)){
              stopifnot(class(y)=="HybridPop")
              x@nInd = x@nInd+y@nInd
              x@id = c(x@id,y@id)
              x@mother = c(x@mother,y@mother)
              x@father = c(x@father,y@father)
              x@gv = rbind(x@gv,y@gv)
              x@pheno = rbind(x@pheno,y@pheno)
              if(x@nTraits>=1){
                for(trait in 1:x@nTraits){
                  if(!is.null(x@gxe[[trait]])){
                    x@gxe[[trait]] = c(x@gxe[[trait]],y@gxe[[trait]])
                  }
                }
              }
            }
            validObject(x)
            return(x)
          }
)
