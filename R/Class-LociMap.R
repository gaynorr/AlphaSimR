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

#TraitAE----
#' @title Additive and epistatic trait
#' 
#' @description Extends \code{\link{TraitA-class}} 
#' to add epistasis
#' 
#' @slot epiEff epistatic effects
#'
#' @export
setClass("TraitAE",
         slots=c(epiEff="matrix"),
         contains="TraitA")

setValidity("TraitAE",function(object){
  errors = character()
  if(object@nLoci!=(2*nrow(object@epiEff))){
    errors = c(errors,"nLoci!=2*nrow(epiEff)")
  }
  if(ncol(object@epiEff)!=3){
    errors = c(errors,"ncol(epiEff)!=3")
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

#TraitADE----
#' @title Additive, dominance, and epistatic trait
#' 
#' @description Extends \code{\link{TraitAD-class}} 
#' to add epistasis
#' 
#' @slot epiEff epistatic effects
#'
#' @export
setClass("TraitADE",
         slots=c(epiEff="matrix"),
         contains="TraitAD")

setValidity("TraitADE",function(object){
  errors = character()
  if(object@nLoci!=(2*nrow(object@epiEff))){
    errors = c(errors,"nLoci!=2*nrow(epiEff)")
  }
  if(ncol(object@epiEff)!=3){
    errors = c(errors,"ncol(epiEff)!=3")
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
#' @slot gxeInt GxE intercept
#' @slot envVar Environmental variance
#'
#' @export
setClass("TraitAG",
         slots=c(gxeEff="numeric",
                 gxeInt="numeric",
                 envVar="numeric"),
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

#TraitAEG----
#' @title Additive, epistasis and GxE trait
#' 
#' @description Extends \code{\link{TraitAE-class}} 
#' to add GxE effects
#' 
#' @slot gxeEff GxE effects
#' @slot gxeInt GxE intercept
#' @slot envVar Environmental variance
#'
#' @export
setClass("TraitAEG",
         slots=c(gxeEff="numeric",
                 gxeInt="numeric",
                 envVar="numeric"),
         contains="TraitAE")

setValidity("TraitAEG",function(object){
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
#' @slot gxeInt GxE intercept
#' @slot envVar Environmental variance
#'
#' @export
setClass("TraitADG",
         slots=c(gxeEff="numeric",
                 gxeInt="numeric",
                 envVar="numeric"),
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

#TraitADEG----
#' @title Additive, dominance, epistasis, and GxE trait
#' 
#' @description Extends \code{\link{TraitADE-class}} 
#' to add GxE effects
#' 
#' @slot gxeEff GxE effects
#' @slot gxeInt GxE intercept
#' @slot envVar Environmental variance
#'
#' @export
setClass("TraitADEG",
         slots=c(gxeEff="numeric",
                 gxeInt="numeric",
                 envVar="numeric"),
         contains="TraitADE")

setValidity("TraitADEG",function(object){
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
#' @slot Vu Estimated marker variance
#' @slot Ve Estimated error variance
#' @slot LL Log-likelihood
#' @slot iter Number of iterations for convergence
#'
#' @export
setClass("RRsol",
         slots=c(markerEff="matrix",
                 fixEff="matrix",
                 Vu="matrix",
                 Ve="matrix",
                 LL="numeric",
                 iter="numeric"),
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



#RRDsol----
#' @title RR-BLUP Solution with Dominance
#' 
#' @description Extends \code{\link{LociMap-class}} 
#' to contain estimated effects from \code{\link{RRBLUP}}
#' 
#' @slot markerEff GEBVs for markers
#' @slot addEff additive effects
#' @slot domEff dominance effects
#' @slot fixEff Estimates for fixed effects
#' @slot Vu Estimated marker variance
#' @slot Ve Estimated error variance
#' @slot LL Log-likelihood
#' @slot iter Number of iterations for convergence
#'
#' @export
setClass("RRDsol",
         slots=c(markerEff="matrix",
                 addEff="matrix",
                 domEff="matrix",
                 fixEff="matrix",
                 Vu="matrix",
                 Ve="matrix",
                 LL="numeric",
                 iter="numeric"),
         contains="LociMap")

setValidity("RRDsol",function(object){
  errors = character()
  if(!is.numeric(object@markerEff)){
    errors = c(errors,"!is.numeric(markerEff)")
  }
  if(!is.numeric(object@addEff)){
    errors = c(errors,"!is.numeric(addEff)")
  }
  if(!is.numeric(object@domEff)){
    errors = c(errors,"!is.numeric(domEff)")
  }
  if(!is.numeric(object@fixEff)){
    errors = c(errors,"!is.numeric(fixEff)")
  }
  if(object@nLoci!=nrow(object@markerEff)){
    errors = c(errors,"nLoci!=nrow(markerEff)")
  }
  if(object@nLoci!=nrow(object@domEff)){
    errors = c(errors,"nLoci!=nrow(domEff)")
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
#' @slot Vu Estimated marker variances in order: 
#' female effects, male effects and SCA effects 
#' (\code{\link{SCAsol-class}} only)
#' @slot Ve Estimated error variance
#' @slot LL Log-likelihood
#' @slot iter Number of iterations for convergence
#'
#' @export
setClass("GCAsol",
         slots=c(femaleEff="matrix",
                 maleEff="matrix",
                 fixEff="matrix",
                 Vu="matrix",
                 Ve="matrix",
                 LL="numeric",
                 iter="numeric"),
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
#' @slot d dominance effect
#'
#' @export
setClass("SCAsol",
         slots=c(d="matrix"),
         contains="GCAsol")

setValidity("SCAsol",function(object){
  errors = character()
  if(!is.numeric(object@d)){
    errors = c(errors,"!is.numeric(d)")
  }
  if(object@nLoci!=nrow(object@d)){
    errors = c(errors,"nLoci!=nrow(d)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})
