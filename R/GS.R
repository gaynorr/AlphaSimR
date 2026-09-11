#' @title Convert traits to a vector of names
#' 
#' @description This function processes the traits arguments from GS models. It usually 
#' receives a number indicating the trait, but may also be the trait name 
#' itself or a custom function.
#'
#' @param traits the traits argument from a GS model
#' @param simParam simulation parameters. If \code{NULL}, the function uses
#' the object named \code{SP} from the global environment.
#'
#' @returns a vector of names for traits
#' 
#' @keywords internal
convertTraitsToNames = function(traits, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.character(traits)){
    # Suspect trait is a name
    take = match(traits, simParam$traitNames)
    if(any(is.na(take))){
      stop("'",traits[is.na(take)],"' did not match any trait names")
    }
    traits = take
  }else if(is.function(traits)){
    traits = "Custom Function"
  }else{
    traits = simParam$traitNames[traits]
  }
  return(traits)
}

#' @title Fast RR-BLUP
#'
#' @description
#' Solves an RR-BLUP model for genomic predictions. This implementation is
#' meant as a fast and low memory alternative to \code{\link{RRBLUP}} or
#' \code{\link{RRBLUP2}}. Fixed effects are fit in the same way as
#' \code{\link{RRBLUP}}, using the levels of the population's fixEff slot.
#'
#' The mixed model equations are solved by preconditioned conjugate
#' gradient, iterating over the genotypes rather than over a stored
#' coefficient matrix \insertCite{stranden_1999}{AlphaSimR}. Neither the
#' coefficient matrix nor a numeric copy of the genotypes is formed, so
#' memory use stays close to the size of the genotypes, and the iteration
#' is shared across nThreads.
#'
#' Variance components are estimated by REML when they are not supplied.
#' REML needs a matrix decomposition whose cost is the square of the number
#' of records in memory and their cube in time, which is what this function
#' is built to avoid, so the estimate is taken from a random subset of the
#' records set by nSubsample. Variance components are nuisance parameters
#' for genomic prediction and are estimated far more precisely than a
#' breeding program needs, so a few thousand records give values good
#' enough to shrink with. The marker effects are always solved for using
#' every record.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, a trait name,
#' or a function of the traits returning a single value. Only univariate models
#' are supported.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations.
#' @param Vu marker effect variance. If value is NULL, this variance
#' and Ve are estimated by REML. See details.
#' @param Ve error variance. If value is NULL, this variance and Vu
#' are estimated by REML. See details.
#' @param nSubsample the number of records used to estimate the variance
#' components. Ignored when Vu and Ve are supplied. A value of zero uses
#' every record, which is only advisable for small training populations,
#' because the memory needed grows with the square of this value.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = fastRRBLUP(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
fastRRBLUP = function(pop, traits=1, use="pheno", snpChip=1,
                      useQtl=FALSE, maxIter=1000, Vu=NULL, Ve=NULL,
                      nSubsample=5000L, simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  # Sort out Vu and Ve. Both are estimated unless both are supplied.
  estVarComp = is.null(Vu) | is.null(Ve)
  if(estVarComp){
    # Values are not used, but something has to be passed
    Vu = 1
    Ve = 1
    nSubsample = as.integer(nSubsample)
    if((nSubsample>0L) & (nSubsample<nrow(y))){
      # Sampling here rather than in C++ keeps it under set.seed()
      subset = sort(sample.int(nrow(y), nSubsample))
    }else{
      subset = seq_len(nrow(y))
    }
  }else{
    subset = 1L
  }

  #Fit model
  ans = callFastRRBLUP(y,fixEff,pop@geno,lociPerChr,
                       lociLoc,Vu,Ve,maxIter,
                       estVarComp,as.integer(subset),
                       nThreads)

  bv = new("TraitA",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$alpha),
           intercept=c(ans$beta),
           name=paste0("est_BV_",traits))

  gv = new("TraitA",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$alpha),
           intercept=c(ans$mu),
           name=paste0("est_GV_",traits))

  output = new("RRsol",
               bv = list(bv),
               gv = list(gv),
               female = as.list(NULL),
               male = as.list(NULL),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}



#' @title RR-BLUP Model
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait or traits to model, a vector of trait names,
#' or a function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. Only used
#' when number of traits is greater than 1.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
RRBLUP = function(pop, traits=1, use="pheno", snpChip=1,
                  useQtl=FALSE, maxIter=1000L,
                  simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  #Fit model
  if(ncol(y)>1){
    ans = callRRBLUP_MV(y, fixEff, pop@geno, lociPerChr,
                        lociLoc, maxIter, nThreads)
  }else{
    ans = callRRBLUP(y, fixEff, pop@geno, lociPerChr, lociLoc,
                     nThreads)
  }

  markerEff=ans$u

  bv = gv = vector("list",ncol(y))
  
  for(i in seq_len(ncol(y))){
    bv[[i]] = new("TraitA",
                  nLoci=nLoci,
                  lociPerChr=lociPerChr,
                  lociLoc=lociLoc,
                  addEff=ans$alpha[,i],
                  intercept=ans$beta[i],
                  name=paste0("est_BV_",traits[i]))

    gv[[i]] = new("TraitA",
                  nLoci=nLoci,
                  lociPerChr=lociPerChr,
                  lociLoc=lociLoc,
                  addEff=ans$alpha[,i],
                  intercept=ans$mu[i],
                  name=paste0("est_GV_",traits[i]))
  }

  output = new("RRsol",
               bv = bv,
               gv = gv,
               female = as.list(NULL),
               male = as.list(NULL),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}

#' @title RR-BLUP Model 2
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions. This implementation is
#' meant for situations where \code{\link{RRBLUP}} is too slow. Note that
#' RRBLUP2 is only faster in certain situations, see details below. Most
#' users should use \code{\link{RRBLUP}}.
#'
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, a trait name, or a
#' function of the traits returning a single value. Unlike \code{\link{RRBLUP}},
#' only univariate models are supported.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations.
#' @param Vu marker effect variance. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param Ve error variance. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param useEM use EM to solve variance components. If false,
#' the initial values are considered true.
#' @param tol tolerance for EM algorithm convergence
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @details
#' The RRBLUP2 function works best when the number of markers is not
#' too large. This is because it solves the RR-BLUP problem by setting
#' up and solving Henderson's mixed model equations. Solving these equations
#' involves a square matrix with dimensions equal to the number of fixed
#' effects plus the number of random effects (markers). Whereas the \code{\link{RRBLUP}}
#' function solves the RR-BLUP problem using the EMMA approach. This approach involves
#' a square matrix with dimensions equal to the number of phenotypic records. This means
#' that the RRBLUP2 function uses less memory than RRBLUP when the number of markers
#' is approximately equal to or smaller than the number of phenotypic records.
#'
#' The RRBLUP2 function is not recommend for cases where the variance components are
#' unknown. This is uses the EM algorithm to solve for unknown variance components,
#' which is generally considerably slower than the EMMA approach of \code{\link{RRBLUP}}.
#' The number of iterations for the EM algorithm is set by maxIter. The default value
#' is typically too small for convergence. When the algorithm fails to converge a
#' warning is displayed, but results are given for the last iteration. These results may
#' be "good enough". However we make no claim to this effect, because we can not generalize
#' to all possible use cases.
#'
#' The RRBLUP2 function can quickly solve the mixed model equations without estimating variance
#' components. The variance components are set by defining Vu and Ve. Estimation of components
#' is suppressed by setting useEM to false. This may be useful if the model is being retrained
#' multiple times during the simulation. You could run \code{\link{RRBLUP}} function the first
#' time the model is trained, and then use the variance components from this output for all
#' future runs with the RRBLUP2 functions. Again, we can make no claim to the general robustness
#' of this approach.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP2(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
RRBLUP2 = function(pop, traits=1, use="pheno", snpChip=1,
                   useQtl=FALSE, maxIter=10, Vu=NULL, Ve=NULL,
                   useEM=TRUE, tol=1e-6, simParam=NULL,
                   nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  # Sort out Vu and Ve
  if(is.function(traits)){
    if(is.null(Vu)){
      Vu = var(y)/nLoci
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot(length(traits)==1)
    if(is.null(Vu)){
      Vu = 2*simParam$varA[traits]/nLoci
      if(is.na(Vu)){
        Vu = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traits]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }

  #Fit model
  ans = callRRBLUP2(y, fixEff, pop@geno, lociPerChr,
                    lociLoc, Vu, Ve, tol, maxIter, useEM,
                    nThreads)

  bv = new("TraitA",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$alpha),
           intercept=c(ans$beta),
           name=paste0("est_BV_",traits))

  gv = new("TraitA",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$alpha),
           intercept=c(ans$mu),
           name=paste0("est_GV_",traits))

  output = new("RRsol",
               bv = list(bv),
               gv = list(gv),
               female = as.list(NULL),
               male = as.list(NULL),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}

#' @title RR-BLUP Model with Dominance
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions that includes
#' dominance effects.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, a trait name, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. Only used
#' when number of traits is greater than 1.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP_D(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
RRBLUP_D = function(pop, traits=1, use="pheno", snpChip=1,
                    useQtl=FALSE, maxIter=40L,
                    simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_D(y, fixEff, pop@geno, lociPerChr,
                     lociLoc, maxIter, nThreads)

  bv = new("TraitA",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$alpha),
           intercept=c(ans$beta),
           name=paste0("est_BV_",traits))

  gv = new("TraitAD",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$a),
           domEff=c(ans$d),
           intercept=c(ans$mu),
           name=paste0("est_GV_",traits))

  output = new("RRsol",
               bv = list(bv),
               gv = list(gv),
               female = as.list(NULL),
               male = as.list(NULL),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}


#' @title RR-BLUP with Dominance Model 2
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions that includes
#' dominance effects. This implementation is meant for situations where
#' \code{\link{RRBLUP_D}} is too slow. Note that RRBLUP_D2
#' is only faster in certain situations. Most users should use
#' \code{\link{RRBLUP_D}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, a trait name, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. Only used
#' when number of traits is greater than 1.
#' @param Va marker effect variance for additive effects. If value is NULL,
#' a reasonable starting point is chosen automatically.
#' @param Vd marker effect variance for dominance effects. If value is NULL,
#' a reasonable starting point is chosen automatically.
#' @param Ve error variance. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param useEM use EM to solve variance components. If false,
#' the initial values are considered true.
#' @param tol tolerance for EM algorithm convergence
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP_D2(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
RRBLUP_D2 = function(pop, traits=1, use="pheno", snpChip=1,
                     useQtl=FALSE, maxIter=10, Va=NULL, Vd=NULL,
                     Ve=NULL, useEM=TRUE, tol=1e-6,
                     simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  # Sort out Va, Vd and Ve
  if(is.function(traits)){
    if(is.null(Va)){
      Va = var(y)/nLoci
    }
    if(is.null(Vd)){
      Vd = var(y)/nLoci
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot(length(traits)==1)
    if(is.null(Va)){
      Va = 2*simParam$varA[traits]/nLoci
      if(is.na(Va)){
        Va = var(y)/nLoci
      }
    }
    if(is.null(Vd)){
      Vd = 2*simParam$varA[traits]/nLoci
      if(is.na(Vd)){
        Vd = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traits]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }

  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_D2(y, fixEff, pop@geno, lociPerChr,
                      lociLoc, maxIter, Va, Vd, Ve, tol, useEM,
                      nThreads)

  bv = new("TraitA",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$alpha),
           intercept=c(ans$beta),
           name=paste0("est_BV_",traits))

  gv = new("TraitAD",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$a),
           domEff=c(ans$d),
           intercept=c(ans$mu),
           name=paste0("est_GV_",traits))

  output = new("RRsol",
               bv = list(bv),
               gv = list(gv),
               female = as.list(NULL),
               male = as.list(NULL),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}

#' @title RR-BLUP GCA Model
#'
#' @description
#' Fits an RR-BLUP model that estimates seperate marker effects for
#' females and males. Useful for predicting GCA of parents
#' in single cross hybrids. Can also predict performance of specific
#' single cross hybrids.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, a trait name, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations for convergence.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP_GCA(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
RRBLUP_GCA = function(pop, traits=1, use="pheno", snpChip=1,
                      useQtl=FALSE, maxIter=40L,
                      simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_GCA(y, fixEff, pop@geno,
                       lociPerChr, lociLoc, maxIter,
                       nThreads)

  gv = new("TraitA2",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$alpha1),
           addEffMale=c(ans$alpha2),
           intercept=c(ans$mu),
           name=paste0("est_GV_",traits))

  female = new("TraitA",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               addEff=c(ans$alpha1),
               intercept=c(ans$beta1),
               name=paste0("est_female_",traits))

  male = new("TraitA",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               addEff=c(ans$alpha2),
               intercept=c(ans$beta2),
             name=paste0("est_male_",traits))

  output = new("RRsol",
               gv = list(gv),
               bv = as.list(NULL),
               female = list(female),
               male = list(male),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}

#' @title RR-BLUP GCA Model 2
#'
#' @description
#' Fits an RR-BLUP model that estimates seperate marker effects for
#' females and males. This implementation is meant for situations where
#' \code{\link{RRBLUP_GCA}} is too slow. Note that RRBLUP_GCA2
#' is only faster in certain situations. Most users should use
#' \code{\link{RRBLUP_GCA}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, a trait name, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations for convergence.
#' @param VuF marker effect variance for females. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param VuM marker effect variance for males. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param Ve error variance. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param useEM use EM to solve variance components. If false,
#' the initial values are considered true.
#' @param tol tolerance for EM algorithm convergence
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP_GCA2(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
RRBLUP_GCA2 = function(pop, traits=1, use="pheno", snpChip=1,
                       useQtl=FALSE, maxIter=10, VuF=NULL, VuM=NULL,
                       Ve=NULL, useEM=TRUE, tol=1e-6,
                       simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  # Sort out VuF, VuM and Ve
  if(is.function(traits)){
    if(is.null(VuF)){
      VuF = var(y)/nLoci
    }
    if(is.null(VuM)){
      VuM = var(y)/nLoci
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot(length(traits)==1)
    if(is.null(VuF)){
      VuF = 2*simParam$varA[traits]/nLoci
      if(is.na(VuF)){
        VuF = var(y)/nLoci
      }
    }
    if(is.null(VuM)){
      VuM = 2*simParam$varA[traits]/nLoci
      if(is.na(VuM)){
        VuM = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traits]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }

  #Fit model
  stopifnot(ncol(y)==1)

  ans = callRRBLUP_GCA2(y, fixEff, pop@geno,
                        lociPerChr, lociLoc, maxIter,
                        VuF, VuM, Ve, tol, useEM,
                        nThreads)

  gv = new("TraitA2",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$alpha1),
           addEffMale=c(ans$alpha2),
           intercept=c(ans$mu),
           name=paste0("est_GV_",traits))

  female = new("TraitA",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               addEff=c(ans$alpha1),
               intercept=c(ans$beta1),
               name=paste0("est_female_",traits))

  male = new("TraitA",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               addEff=c(ans$alpha2),
               intercept=c(ans$beta2),
             name=paste0("est_male_",traits))

  output = new("RRsol",
               gv = list(gv),
               bv = as.list(NULL),
               female = list(female),
               male = list(male),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}

#' @title RR-BLUP SCA Model
#'
#' @description
#' An extention of \code{\link{RRBLUP_GCA}} that adds dominance effects.
#' Note that we have not seen any consistent benefit of this model over
#' \code{\link{RRBLUP_GCA}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, a trait name, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations for convergence.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP_SCA(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
RRBLUP_SCA = function(pop, traits=1, use="pheno", snpChip=1,
                      useQtl=FALSE, maxIter=40L,
                      simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_SCA(y, fixEff, pop@geno,
                       lociPerChr, lociLoc, maxIter,
                       nThreads)

  gv = new("TraitA2D",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$a1),
           addEffMale=c(ans$a2),
           domEff=c(ans$d),
           intercept=c(ans$mu),
           name=paste0("est_GV_",traits))

  female = new("TraitA",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               addEff=c(ans$alpha1),
               intercept=c(ans$beta1),
               name=paste0("est_female_",traits))

  male = new("TraitA",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               addEff=c(ans$alpha2),
               intercept=c(ans$beta2),
             name=paste0("est_male_",traits))

  output = new("RRsol",
               gv = list(gv),
               bv = as.list(NULL),
               female = list(female),
               male = list(male),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}

#' @title RR-BLUP SCA Model 2
#'
#' @description
#' Fits an RR-BLUP model that estimates seperate additive effects for
#' females and males and a dominance effect. This implementation is meant
#' for situations where \code{\link{RRBLUP_SCA}} is too slow. Note that
#' RRBLUP_SCA2 is only faster in certain situations. Most users should use
#' \code{\link{RRBLUP_SCA}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, a trait name, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations for convergence.
#' @param VuF marker effect variance for females. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param VuM marker effect variance for males. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param VuD marker effect variance for dominance. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param Ve error variance. If value is NULL, a
#' reasonable starting point is chosen automatically.
#' @param useEM use EM to solve variance components. If false,
#' the initial values are considered true.
#' @param tol tolerance for EM algorithm convergence
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP_SCA2(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
RRBLUP_SCA2 = function(pop, traits=1, use="pheno", snpChip=1,
                       useQtl=FALSE, maxIter=10, VuF=NULL, VuM=NULL,
                       VuD=NULL, Ve=NULL, useEM=TRUE, tol=1e-6,
                       simParam=NULL, nThreads=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  traits = convertTraitsToNames(traits, simParam)

  fixEff = as.integer(factor(pop@fixEff))

  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }

  # Sort out VuF, VuM, VuD and Ve
  if(is.function(traits)){
    if(is.null(VuF)){
      VuF = var(y)/nLoci
    }
    if(is.null(VuM)){
      VuM = var(y)/nLoci
    }
    if(is.null(VuD)){
      VuD = var(y)/nLoci/2
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot(length(traits)==1)
    if(is.null(VuF)){
      VuF = 2*simParam$varA[traits]/nLoci
      if(is.na(VuF)){
        VuF = var(y)/nLoci
      }
    }
    if(is.null(VuM)){
      VuM = 2*simParam$varA[traits]/nLoci
      if(is.na(VuM)){
        VuM = var(y)/nLoci
      }
    }
    if(is.null(VuD)){
      VuD = simParam$varA[traits]/nLoci
      if(is.na(VuD)){
        VuD = var(y)/nLoci/2
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traits]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }

  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_SCA2(y, fixEff, pop@geno,
                        lociPerChr, lociLoc, maxIter,
                        VuF, VuM, VuD, Ve, tol, useEM,
                        nThreads)

  gv = new("TraitA2D",
           nLoci=nLoci,
           lociPerChr=lociPerChr,
           lociLoc=lociLoc,
           addEff=c(ans$a1),
           addEffMale=c(ans$a2),
           domEff=c(ans$d),
           intercept=c(ans$mu),
           name=paste0("est_GV_",traits))

  female = new("TraitA",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               addEff=c(ans$alpha1),
               intercept=c(ans$beta1),
               name=paste0("est_female_",traits))

  male = new("TraitA",
             nLoci=nLoci,
             lociPerChr=lociPerChr,
             lociLoc=lociLoc,
             addEff=c(ans$alpha2),
             intercept=c(ans$beta2),
             name=paste0("est_male_",traits))

  output = new("RRsol",
               gv = list(gv),
               bv = as.list(NULL),
               female = list(female),
               male = list(male),
               Vu = as.matrix(ans$Vu),
               Ve = as.matrix(ans$Ve))

  return(output)
}

#' @title Set estimated breeding values (EBV)
#'
#' @description
#' Adds genomic estimated values to a populations's EBV
#' slot using output from a genomic selection functions.
#' The genomic estimated values can be either estimated
#' breeding values, estimated genetic values, or
#' estimated general combining values.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param solution an object of \code{\link{RRsol-class}}
#' @param value the genomic value to be estimated. Can be
#' either "gv", "bv", "female", or "male".
#' @param targetPop an optional target population that can
#' be used when value is "bv", "female", or "male". When
#' supplied, the allele frequency in the targetPop is used
#' to set these values.
#' @param append should estimated values be appended to
#' existing data in the EBV slot. If TRUE, a new column is
#' added. If FALSE, existing data is replaced with the
#' new estimates.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=20)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model and set EBV
#' ans = RRBLUP(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' @export
setEBV = function(pop, solution, value="gv", targetPop=NULL,
                  append=FALSE, simParam=NULL, nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  nTraits = length(solution@gv)

  ebv = matrix(NA_real_,
               nrow=pop@nInd,
               ncol=nTraits)

  # Placeholder names
  colnames(ebv) = as.character(1:nTraits)

  value = tolower(value)

  if(value=="gv"){
    for(i in seq_len(nTraits)){
      tmp = getGv(solution@gv[[i]],pop,nThreads)
      ebv[,i] = tmp[[1]]
      colnames(ebv)[i] = solution@gv[[i]]@name
    }

  }else if(value=="bv"){

    if(is.null(targetPop)){

      if(length(solution@bv)==0){
        stop("This genomic selection model does not produce breeding value estimates.")
      }
      
      for(i in seq_len(nTraits)){
        tmp = getGv(solution@bv[[i]],pop,nThreads)
        ebv[,i] = tmp[[1]]
        colnames(ebv)[i] = solution@bv[[i]]@name
      }

    }else{
      
      for(i in seq_len(nTraits)){
        trait = solution@gv[[i]]
        if(.hasSlot(trait,"addEffMale")){
          stop("This genomic selection model does not produce breeding value estimates. Try value='male' or value='female' instead.")
        }

        p = calcGenoFreq(targetPop@geno,
                         trait@lociPerChr,
                         trait@lociLoc,
                         nThreads)
        p = c(p)
        q = 1-p

        a = trait@addEff
        if(.hasSlot(trait,"domEff")){
          d = trait@domEff
        }else{
          d = rep(0, length(a))
        }

        alpha = a+d*(q-p)
        intercept = -sum((p-q)*alpha)
        trait = new("TraitA",
                    nLoci=trait@nLoci,
                    lociPerChr=trait@lociPerChr,
                    lociLoc=trait@lociLoc,
                    addEff=alpha,
                    intercept=intercept)

        tmp = getGv(trait, pop, nThreads)
        ebv[,i] = tmp[[1]]

        # changing original name from "est_GV_..." to "est_BV_..."
        tmp = solution@gv[[i]]@name
        tmp = strsplit(tmp, "_")[[1]]
        tmp[2] = "BV"
        colnames(ebv)[i] = paste(tmp,collapse="_")
      }

    }

  }else if(value=="female"){

    if(is.null(targetPop)){

      if(length(solution@female)==0){
        stop("This genomic selection model does not produce GCA estimates for females.")
      }

      for(i in seq_len(nTraits)){
        tmp = getGv(solution@female[[i]],pop,nThreads)
        ebv[,i] = tmp[[1]]
        colnames(ebv)[i] = solution@female[[i]]@name
      }

    }else{

      for(i in seq_len(nTraits)){
        trait = solution@gv[[i]]
        p = calcGenoFreq(targetPop@geno,
                         trait@lociPerChr,
                         trait@lociLoc,
                         nThreads)
        p = c(p)
        q = 1-p

        a = trait@addEff
        if(.hasSlot(trait,"domEff")){
          d = trait@domEff
        }else{
          d = rep(0, length(a))
        }

        alpha = (a+d*(q-p))/2
        intercept = -sum((p-q)*alpha)
        trait = new("TraitA",
                    nLoci=trait@nLoci,
                    lociPerChr=trait@lociPerChr,
                    lociLoc=trait@lociLoc,
                    addEff=alpha,
                    intercept=intercept)

        tmp = getGv(trait, pop, nThreads)
        ebv[,i] = tmp[[1]]

        # changing original name from "est_GV_..." to "est_female_..."
        tmp = solution@gv[[i]]@name
        tmp = strsplit(tmp, "_")[[1]]
        tmp[2] = "female"
        colnames(ebv)[i] = paste(tmp,collapse="_")
      }

    }

  }else if(value=="male"){

    if(is.null(targetPop)){

      if(length(solution@male)==0){
        stop("This genomic selection model does not produce GCA estimates for males.")
      }
      
      for(i in seq_len(nTraits)){
        tmp = getGv(solution@male[[i]],pop,nThreads)
        ebv[,i] = tmp[[1]]
        colnames(ebv)[i] = solution@male[[i]]@name
      }

    }else{

      for(i in seq_len(nTraits)){
        trait = solution@gv[[i]]
        p = calcGenoFreq(targetPop@geno,
                         trait@lociPerChr,
                         trait@lociLoc,
                         nThreads)
        p = c(p)
        q = 1-p

        if(.hasSlot(trait,"addEffMale")){
          a = trait@addEffMale
        }else{
          a = trait@addEff
        }
        if(.hasSlot(trait,"domEff")){
          d = trait@domEff
        }else{
          d = rep(0, length(a))
        }

        alpha = (a+d*(q-p))/2
        intercept = -sum((p-q)*alpha)
        trait = new("TraitA",
                    nLoci=trait@nLoci,
                    lociPerChr=trait@lociPerChr,
                    lociLoc=trait@lociLoc,
                    addEff=alpha,
                    intercept=intercept)

        tmp = getGv(trait, pop, nThreads)
        ebv[,i] = tmp[[1]]

        # changing original name from "est_BV_..." to "est_male_..."
        tmp = solution@gv[[i]]@name
        tmp = strsplit(tmp, "_")[[1]]
        tmp[2] = "male"
        colnames(ebv)[i] = paste(tmp,collapse="_")
      }

    }

  }else{
    stop(paste0("value=",value," is not a valid option"))
  }

  if(append){
    pop@ebv = cbind(pop@ebv,ebv)
  }else{
    pop@ebv = ebv
  }

  return(pop)
}

#' @title RRBLUP Memory Usage
#'
#' @description
#' Estimates the amount of RAM needed to fit one of AlphaSimR's genomic
#' selection models to a training population of a given size. The estimate
#' covers the matrices the solvers hold at their peak, which is what decides
#' whether a model can be fitted at all. It does not cover the population
#' object itself or anything else in the R session, so it is a lower bound
#' on what the whole simulation needs.
#'
#' @param nInd the number of individuals in the training population
#' @param nMarker the number of markers per individual
#' @param model the model being fitted, given as the name of the function
#' that fits it. One of "fastRRBLUP", "RRBLUP", "RRBLUP2", "RRBLUP_D",
#' "RRBLUP_D2", "RRBLUP_GCA", "RRBLUP_GCA2", "RRBLUP_SCA" or "RRBLUP_SCA2".
#' The older values "REG", "GCA" and "SCA" are still accepted and are read as
#' "RRBLUP", "RRBLUP_GCA" and "RRBLUP_SCA".
#' @param nTraits the number of traits fitted at once. Only
#' \code{\link{RRBLUP}} fits more than one.
#' @param nFixEff the number of fixed effect levels, which is the number of
#' distinct values in the population's fixEff slot.
#' @param nSubsample the number of records \code{\link{fastRRBLUP}} uses to
#' estimate variance components. Ignored by every other model, and by
#' fastRRBLUP itself when Vu and Ve are supplied.
#'
#' @details
#' The models differ in what they have to hold in memory, and the differences
#' are large enough to decide which one is usable.
#'
#' \code{\link{fastRRBLUP}} keeps the genotypes as one byte per locus and
#' iterates over them, so it holds no square matrix at all. Its estimate is
#' dominated by the genotypes themselves and by the subset of records used
#' for the variance components.
#'
#' \code{\link{RRBLUP}} decomposes a square matrix whose dimensions are the
#' smaller of the number of records and the number of markers, because a
#' matrix and its transpose share their nonzero eigenvalues.
#'
#' The numbered models set up Henderson's mixed model equations, whose
#' coefficient matrix is square with dimensions equal to the number of fixed
#' effects plus the number of random effects, so they grow with the number of
#' markers rather than the number of records.
#'
#' The GCA, SCA and dominance models fit more than one random effect and hold
#' a square matrix of the number of records for each, which makes them the
#' most demanding of the set.
#'
#' @return Returns an estimate for the required gigabytes of RAM
#'
#' @examples
#' RRBLUPMemUse(nInd=1000, nMarker=5000)
#'
#' # The same training data, fitted three ways
#' RRBLUPMemUse(nInd=5000, nMarker=2000, model="fastRRBLUP")
#' RRBLUPMemUse(nInd=5000, nMarker=2000, model="RRBLUP")
#' RRBLUPMemUse(nInd=5000, nMarker=2000, model="RRBLUP_SCA")
#'
#' @export
RRBLUPMemUse = function(nInd, nMarker, model="RRBLUP", nTraits=1L,
                        nFixEff=1L, nSubsample=5000L){
  n = as.double(nInd)
  m = as.double(nMarker)
  q = as.double(nFixEff)
  nT = as.double(nTraits)

  # Names used before the models were named after their functions
  model = switch(toupper(model),
                 "REG" = "RRBLUP",
                 "GCA" = "RRBLUP_GCA",
                 "SCA" = "RRBLUP_SCA",
                 model)

  # Genotypes are read as one byte per locus before being turned into
  # dosages, so both are held while the conversion runs
  genoBytes = n*m

  if(model=="fastRRBLUP"){
    # No dosages are stored and no square matrix is formed. The solver holds
    # six vectors of markers plus the column means, and a few of records.
    doubles = 7*m+3*n+q*n
    # Variance components come from a subset of the records. Whichever cross
    # product is smaller is the one decomposed.
    nSub = as.double(nSubsample)
    if((nSub<=0) | (nSub>n)){
      nSub = n
    }
    if(m<nSub){
      doubles = doubles+nSub*m+5*m^2
    }else{
      doubles = doubles+3*nSub^2
    }
  }else if(model=="RRBLUP"){
    if(nT>1){
      # Rotation by the eigenvectors of the record cross product, plus one
      # small inverse per record
      doubles = n*m+3*n^2+nT^2*n
    }else if(m<n){
      # Decomposition over markers, with the fixed effects projected out of
      # a second copy of the dosages
      doubles = 2*n*m+3*m^2
    }else{
      # Decomposition over records
      doubles = n*m+4*n^2
    }
  }else if(model=="RRBLUP2"){
    doubles = n*m+2*(q+m)^2
  }else if(model %in% c("RRBLUP_D","RRBLUP_GCA","RRBLUP_SCA")){
    nKernel = if(model=="RRBLUP_SCA") 3 else 2
    if(m<n){
      # The average information matrix is reached through the dosages
      doubles = 2*nKernel*n*m+(nKernel+3)*n^2
    }else{
      # and through square matrices of the records
      doubles = nKernel*n*m+(3*nKernel+3)*n^2
    }
  }else if(model %in% c("RRBLUP_D2","RRBLUP_GCA2","RRBLUP_SCA2")){
    nKernel = if(model=="RRBLUP_SCA2") 3 else 2
    doubles = nKernel*n*m+2*(q+nKernel*m)^2
  }else{
    stop(paste0("model=",model," not recognized"))
  }

  # Using base 10 rather than base 2, which is the more conservative of the
  # two and leaves a little room for what is not counted here
  return((8*doubles+genoBytes)/10^9) #GB
}
