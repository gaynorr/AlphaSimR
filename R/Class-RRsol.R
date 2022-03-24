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
#' @slot Vu Estimated marker variance(s)
#' @slot Ve Estimated error variance
#'
#' @export
setClass("RRsol",
         slots=c(gv="list",
                 bv="list",
                 female="list",
                 male="list",
                 Vu="matrix",
                 Ve="matrix"))

#' @describeIn RRsol Test if object is of a RRsol class
#' @export
isRRsol = function(x) {
  ret = is(x, class2 = "RRsol")
  return(ret)
}
