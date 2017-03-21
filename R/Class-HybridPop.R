#HybridPop----
#' @title Population
#' 
#' @description Population class
#' 
#' @slot nInd number of individuals
#' @slot gv matrix of genotypic values
#' @slot pheno matrix of phenotypic values
#' @slot id hybrid's id
#' @slot par1 hybrid's first parent/mother
#' @slot par2 hybrid's second parent/father
#' 
#' @export
setClass("HybridPop",
         slots=c(nInd="integer",
                 gv="matrix",
                 pheno="matrix",
                 id="character",
                 par1="character",
                 par2="character"))


