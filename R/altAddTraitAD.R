#' @title Alternative method for adding an AD trait
#' 
#' @description An alternative method for adding a trait with additive 
#' and dominance effects to an AlphaSimR simulation. The function attempts 
#' to create a trait matching user defined values for number of QTL, inbreeding 
#' depression, additive genetic variance and dominance genetic variance.
#' 
#' @param nQtlPerChr number of QTLs per chromosome. Can be a single value or nChr values.
#' @param mean desired mean of the trait
#' @param varA desired additive variance
#' @param varD desired dominance variance
#' @param inbrDepr desired inbreeding depression
#' @param limMeanDD limits for meanDD, see details
#' @param limVarDD limits for varDD, see details
#' @param name optional name for trait
#' @param simParam an object of 'SimParam' class
#' 
#' @details 
#' TO DO
#' 
#' @return NULL, the trait is added direction to simParam
#' 
#' @export
altAddTraitAD = function(nQtlPerChr, mean, varA, varD, inbrDepr, 
                          limMeanDD = c(0, 2), 
                          limVarDD = c(0, 0.3),
                          name = NULL, simParam = NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  
  # a_std, d_mu, d_std
}
