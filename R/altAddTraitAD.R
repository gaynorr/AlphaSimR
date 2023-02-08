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
#' @param silent should summary details be printed to the console
#' @param simParam an object of 'SimParam' class
#' 
#' @details 
#' TO DO
#' 
#' @return NULL, the trait is added directly to simParam
#' 
#' @export
altAddTraitAD = function(nQtlPerChr, mean, varA, varD, inbrDepr, 
                         limMeanDD = c(0, 1.5), 
                         limVarDD = c(0, 0.3),
                         name = NULL, silent=FALSE, 
                         simParam = NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  # Check that simParam still has it's founderPop
  stopifnot("'simParam$founderPop' must have a 'RawPop'" = 
              is(simParam$founderPop, "RawPop"))
  
  # Add a placeholder trait
  simParam$addTraitAD(nQtlPerChr = nQtlPerChr, 
                      name = name)
  
  # Initialize the tuner object using objects from simParam
  tuner = TuneAD$new(LociMap = simParam$traits[[simParam$nTraits]],
                     Pop = simParam$founderPop,
                     mean_ = mean,
                     varA_ = varA,
                     varD_ = varD,
                     inbrDepr_ = inbrDepr,
                     nThreads_ = simParam$nThreads)
  
  # Run optim to optimize meanDD and varDD
  optOut = optim(par = c(mean(limMeanDD), mean(sqrt(limVarDD))),
                 fn = tuner$objective, 
                 gr = NULL,
                 method = "L-BFGS-B",
                 lower = c(limMeanDD[1], sqrt(limVarDD[1])),
                 upper = c(limMeanDD[2], sqrt(limVarDD[2])))
  
  # Sumarize parameters
  meanDD = optOut$par[1]
  varDD = optOut$par[2]^2
  output = tuner$finalize(meanDD_ = meanDD, 
                          stdDevDD_ = sqrt(varDD))
  
  # Set new trait
  nTraits = simParam$nTraits
  trait = simParam$traits[[nTraits]]
  trait@addEff = c(output$a)
  trait@domEff = c(output$d)
  trait@intercept = c(output$intercept)
  simParam$switchTrait(nTraits, trait)
  
  # Report trait details
  if(!silent){
    cat("New trait called", simParam$traitNames[nTraits], "was added \n")
    cat("Dominance variance is", output$varD, "\n")
    cat("Inbreeding depression is", output$inbrDepr, "\n")
    cat("Used meanDD equals", meanDD, "\n")
    cat("Used varDD equals", varDD, "\n")
  }
  
  rm(tuner)
  
  return(invisible())
}
