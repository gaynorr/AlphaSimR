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
#' @slot name optional name for LociMap object
#'
#' @export
setClass("LociMap",
         slots=c(nLoci="integer",
                 lociPerChr="integer",
                 lociLoc="integer",
                 name="character"),
         prototype=list(name=NA_character_))

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

# Test if object is of LociMap class
isLociMap = function(x) {
  ret = is(x, class2 = "LociMap")
  return(ret)
}

#TraitA----
#' @title Additive trait
#'
#' @description Extends \code{\link{LociMap-class}}
#' to model additive traits
#'
#' @slot addEff additive effects
#' @slot intercept adjustment factor for gv
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

# Test if object is of a TraitA class
isTraitA = function(x) {
  ret = is(x, class2 = "TraitA")
  return(ret)
}

#TraitA2----
#' @title Sex specific additive trait
#'
#' @description Extends \code{\link{TraitA-class}}
#' to model seperate additive effects for parent of
#' origin. Used exclusively for genomic selection.
#'
#' @slot addEffMale additive effects
#'
#' @export
setClass("TraitA2",
         slots=c(addEffMale="numeric"),
         contains="TraitA")

setValidity("TraitA2",function(object){
  errors = character()
  if(object@nLoci!=length(object@addEffMale)){
    errors = c(errors,"nLoci!=length(addEffMale)")
  }
  if(length(errors)==0){
    return(TRUE)
  }else{
    return(errors)
  }
})

# Test if object is of a TraitA2 class
isTraitA2 = function(x) {
  ret = is(x, class2 = "TraitA2")
  return(ret)
}

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

# Test if object is of a TraitAE class
isTraitAE = function(x) {
  ret = is(x, class2 = "TraitAE")
  return(ret)
}

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

# Test if object is of a TraitAD class
isTraitAD = function(x) {
  ret = is(x, class2 = "TraitAD")
  return(ret)
}

#TraitA2D----
#' @title Sex specific additive and dominance trait
#'
#' @description Extends \code{\link{TraitA2-class}}
#' to add dominance
#'
#' @slot domEff dominance effects
#'
#' @export
setClass("TraitA2D",
         slots=c(domEff="numeric"),
         contains="TraitA2")

setValidity("TraitA2D",function(object){
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

# Test if object is of a TraitA2D class
isTraitA2D = function(x) {
  ret = is(x, class2 = "TraitA2D")
  return(ret)
}

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

# Test if object is of a TraitADE class
isTraitADE = function(x) {
  ret = is(x, class2 = "TraitADE")
  return(ret)
}

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

# Test if object is of a TraitAG class
isTraitAG = function(x) {
  ret = is(x, class2 = "TraitAG")
  return(ret)
}

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

# Test if object is of a TraitAEG class
isTraitAEG = function(x) {
  ret = is(x, class2 = "TraitAEG")
  return(ret)
}

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

# Test if object is of a TraitADG class
isTraitADG = function(x) {
  ret = is(x, class2 = "TraitADG")
  return(ret)
}

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

# Test if object is of a TraitADEG class
isTraitADEG = function(x) {
  ret = is(x, class2 = "TraitADEG")
  return(ret)
}
