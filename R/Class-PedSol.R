#PedSol----
#' @title Pedigree BLUP solution
#'
#' @description Contains output from \code{\link{pedBLUP}}, which is used
#' with \code{\link{setEBV}} to set estimated breeding values or genetic
#' values.
#'
#' @slot iid the internal IDs of the individuals in the fitted pedigree,
#' which are the phenotyped individuals and all of their ancestors
#' @slot ebv the estimated breeding values of those individuals, as
#' deviations from the base population
#' @slot fixEff the estimated fixed effects, one for each level of the
#' fixEff slot of the training population, named by level
#' @slot center the mean estimated breeding value of the records in the
#' training population, which setEBV removes when value is "bv"
#' @slot trait the name of the trait
#' @slot Vu the estimated additive genetic variance, as a 1 by 1 matrix
#' @slot Ve the estimated error variance, as a 1 by 1 matrix
#' @slot model the model that was fitted, "pedBLUP"
#'
#' @export
setClass("PedSol",
         slots=c(iid="integer",
                 ebv="numeric",
                 fixEff="numeric",
                 center="numeric",
                 trait="character",
                 Vu="matrix",
                 Ve="matrix",
                 model="character"),
         prototype=prototype(model="pedBLUP"))
