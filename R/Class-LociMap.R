# The LociMap superclass contains SNP/QTL locations
# Trait classes add QTL effects

#LociMap----
#' @title Loci metadata
#' 
#' @description used for both SNPs and QTLs
#' 
#' @slot nLoci total number of loci
#' @slot lociPerChr number of loci per chromosome
#' @slot lociLoc physical position of loci
#'
#' @export
setClass("LociMap",
         slots=c(nLoci="integer",
                 lociPerChr="integer",
                 lociLoc="integer"))

setValidity("LociMap",function(object){
  errors = character()
  if(object@nLoci!=sum(object@lociPerChr)){
    errors = c(errors,"nLoci!=sum(lociPerChr)")
  }
  if(object@nLoci!=length(object@lociLoc)){
    errors = c(errors,"nLoci!=length(lociLoc)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#TraitA----
#' @title Additive trait
#' 
#' @description Extends \code{\link{LociMap-class}} 
#' to model additive traits
#' 
#' @slot addEff additive effects
#' @slot intercept adjustment factor for gc
#'
#' @export
setClass("TraitA",
         slots=c(addEff="numeric",
                 intercept="numeric"),
         contains="LociMap")

setValidity("TraitA",function(object){
  errors = character()
  if(object@nLoci!=length(object@addEff)){
    errors = c(errors,"nLoci!=length(addEff)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#TraitAD----
#' @title Additive and dominance trait
#' 
#' @description Extends \code{\link{TraitA-class}} 
#' to add dominance
#' 
#' @slot domEff dominance effects
#'
#' @export
setClass("TraitAD",
         slots=c(domEff="numeric"),
         contains="TraitA")

setValidity("TraitAD",function(object){
  errors = character()
  if(object@nLoci!=length(object@domEff)){
    errors = c(errors,"nLoci!=length(domEff)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#TraitAG----
#' @title Additive and GxE trait
#' 
#' @description Extends \code{\link{TraitA-class}} 
#' to add GxE effects
#' 
#' @slot gxeEff GxE effects
#' @slot varGxeLoci GxE variance for loci
#'
#' @export
setClass("TraitAG",
         slots=c(gxeEff="numeric",
                 varGxeLoci="numeric"),
         contains="TraitA")

setValidity("TraitAG",function(object){
  errors = character()
  if(object@nLoci!=length(object@gxeEff)){
    errors = c(errors,"nLoci!=length(gxeEff)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#TraitADG----
#' @title Additive, dominance and GxE trait
#' 
#' @description Extends \code{\link{TraitAD-class}} 
#' to add GxE effects
#' 
#' @slot gxeEff GxE effects
#' @slot varGxeLoci GxE variance for loci
#'
#' @export
setClass("TraitADG",
         slots=c(gxeEff="numeric",
                 varGxeLoci="numeric"),
         contains="TraitAD")

setValidity("TraitADG",function(object){
  errors = character()
  if(object@nLoci!=length(object@gxeEff)){
    errors = c(errors,"nLoci!=length(gxeEff)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#RRsol----
#' @title RR-BLUP Solution
#' 
#' @description Extends \code{\link{LociMap-class}} 
#' to contain estimated effects from \code{\link{RRBLUP}}
#' 
#' @slot markerEff GEBVs for markers
#' @slot fixEff Estimates for fixed effects
#'
#' @export
setClass("RRsol",
         slots=c(markerEff="matrix",
                 fixEff="matrix"),
         contains="LociMap")

setValidity("RRsol",function(object){
  errors = character()
  if(!is.numeric(object@markerEff)){
    errors = c(errors,"!is.numeric(markerEff)")
  }
  if(!is.numeric(object@fixEff)){
    errors = c(errors,"!is.numeric(fixEff)")
  }
  if(object@nLoci!=nrow(object@markerEff)){
    errors = c(errors,"nLoci!=nrow(markerEff)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#GCAsol----
#' @title RR-BLUP GCA Solution
#' 
#' @description Extends \code{\link{LociMap-class}} 
#' to contain estimated effects from \code{\link{RRBLUP_GCA}}
#' 
#' @slot femaleEff marker GCA for "female" pool
#' @slot maleEff marker GCA for "male" pool
#' @slot fixEff Estimates for fixed effects
#'
#' @export
setClass("GCAsol",
         slots=c(femaleEff="matrix",
                 maleEff="matrix",
                 fixEff="matrix"),
         contains="LociMap")

setValidity("GCAsol",function(object){
  errors = character()
  if(!is.numeric(object@femaleEff)){
    errors = c(errors,"!is.numeric(femaleEff)")
  }
  if(!is.numeric(object@maleEff)){
    errors = c(errors,"!is.numeric(maleEff)")
  }
  if(!is.numeric(object@fixEff)){
    errors = c(errors,"!is.numeric(fixEff)")
  }
  if(object@nLoci!=nrow(object@femaleEff)){
    errors = c(errors,"nLoci!=nrow(femaleEff)")
  }
  if(object@nLoci!=nrow(object@maleEff)){
    errors = c(errors,"nLoci!=nrow(maleEff)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#SCAsol----
#' @title RR-BLUP SCA Solution
#' 
#' @description Extends \code{\link{GCAsol-class}} 
#' to contain estimated effects from \code{\link{RRBLUP_SCA}}
#' 
#' @slot scaEff marker SCA effects
#'
#' @export
setClass("SCAsol",
         slots=c(scaEff="matrix"),
         contains="GCAsol")

setValidity("SCAsol",function(object){
  errors = character()
  if(!is.numeric(object@scaEff)){
    errors = c(errors,"!is.numeric(scaEff)")
  }
  if(object@nLoci!=nrow(object@scaEff)){
    errors = c(errors,"nLoci!=nrow(scaEff)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

#getGv----
setGeneric("getGv",function(object,...){
  standardGeneric("getGv")
})
setMethod("getGv",signature("TraitA"),
          function(object,pop,...){
            getGvA(object,pop)
          })
setMethod("getGv",signature("TraitAD"),
          function(object,pop,...){
            getGvAD(object,pop)
          })
setMethod("getGv",signature("TraitAG"),
          function(object,pop,w){
            z = qnorm(w,sd=sqrt(object@varGxeLoci))
            getGvAG(object,pop,z)
          })
setMethod("getGv",signature("TraitADG"),
          function(object,pop,w){
            z = qnorm(w,sd=sqrt(object@varGxeLoci))
            getGvADG(object,pop,z)
          })
#getHybridGv----
setGeneric("getHybridGv",function(object,...){
  standardGeneric("getHybridGv")
})
setMethod("getHybridGv",signature("TraitA"),
          function(object,fPop,fPar,mPop,mPar,...){
            getHybridGvA(object,fPop,fPar,mPop,mPar)
          })
setMethod("getHybridGv",signature("TraitAD"),
          function(object,fPop,fPar,mPop,mPar,...){
            getHybridGvAD(object,fPop,fPar,mPop,mPar)
          })
setMethod("getHybridGv",signature("TraitAG"),
          function(object,fPop,fPar,mPop,mPar,w){
            z = qnorm(w,sd=sqrt(object@varGxeLoci))
            getHybridGvAG(object,fPop,fPar,mPop,mPar,z)
          })
setMethod("getHybridGv",signature("TraitADG"),
          function(object,fPop,fPar,mPop,mPar,w){
            z = qnorm(w,sd=sqrt(object@varGxeLoci))
            getHybridGvADG(object,fPop,fPar,mPop,mPar,z)
          })
