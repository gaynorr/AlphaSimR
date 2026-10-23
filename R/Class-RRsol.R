#RRsol----
#' @title RR-BLUP Solution
#' 
#' @description Contains output from AlphaSimR's genomic 
#' selection functions.
#' 
#' @slot gv Trait(s) for estimating genetic values
#' @slot bv Trait(s) for estimating breeding values
#' @slot female Trait(s) for estimating GCA in the female pool
#' @slot male Trait(s) for estimating GCA in the male pool
#' @slot Vu Estimated marker variance(s). For a model fitted to several
#' traits, the covariance matrix of the marker effects among traits. For a
#' multiple trait model with several random terms, such as
#' \code{\link{RRBLUP_D}}, \code{\link{RRBLUP_GCA}} or
#' \code{\link{RRBLUP_SCA}}, a block diagonal matrix holding one covariance
#' matrix per term, with rows and columns named "term:trait".
#' @slot Ve Estimated error variance, or for several traits the covariance
#' matrix of the errors among traits
#' @slot model The model that was fitted, named by the function that fits it
#' without a trailing 2, so that \code{\link{RRBLUP}} and
#' \code{\link{RRBLUP2}} both give "RRBLUP". It lets a fit be passed back as
#' varComp to the same model only. Empty when no model is recorded, as for an
#' object created directly.
#'
#' @export
setClass("RRsol",
         slots=c(gv="list",
                 bv="list",
                 female="list",
                 male="list",
                 Vu="matrix",
                 Ve="matrix",
                 model="character"),
         prototype=prototype(model=""))

# Test if object is of a RRsol class
isRRsol = function(x) {
  ret = is(x, class2 = "RRsol")
  return(ret)
}
