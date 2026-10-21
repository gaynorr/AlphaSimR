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
  simParam = getSimParam(simParam)
  if(is.character(traits)){
    # Suspect trait is a name
    take = match(traits, simParam$traitNames)
    if(any(is.na(take))){
      stop("'",traits[is.na(take)],"' did not match any trait names")
    }
    # The function returns names, so a validated name is returned as it
    # stands rather than as the index it happened to match.
    traits = simParam$traitNames[take]
  }else if(is.function(traits)){
    traits = "Custom Function"
  }else{
    traits = simParam$traitNames[traits]
  }
  return(traits)
}

#' @title Convert traits to indices
#'
#' @description Returns the position of each trait within the SimParam, or
#' NULL when the response is a custom function and so has no trait behind
#' it. Vectors such as varA and varE are unnamed and indexed by position,
#' so they need this rather than \code{convertTraitsToNames}.
#'
#' @param traits an integer, a character or a function
#' @param simParam an object of class \code{\link{SimParam}}
#'
#' @keywords internal
convertTraitsToIndex = function(traits, simParam=NULL){
  simParam = getSimParam(simParam)
  if(is.function(traits)){
    return(NULL)
  }
  if(is.character(traits)){
    take = match(traits, simParam$traitNames)
    if(any(is.na(take))){
      stop("'",traits[is.na(take)],"' did not match any trait names")
    }
    return(take)
  }
  return(as.integer(traits))
}

#' @title Check that every trait is recorded
#'
#' @description The multiple trait models need every trait recorded on every
#' individual. Their solvers would stop on a missing value too, but not with
#' a message that says why, so the wrappers call this first. Every multiple
#' trait wrapper goes through it, so supporting missing records later starts
#' here.
#'
#' @param y the matrix of responses
#'
#' @keywords internal
checkCompleteRecords = function(y){
  if(anyNA(y)){
    stop("multiple trait models need every trait recorded on every individual, ",
         "but the responses contain missing values")
  }
  invisible(NULL)
}

#' @title Marker variances of a multiple trait multikernel fit
#'
#' @description Stores the covariance matrices among traits of each random
#' term of a multiple trait fit as one block diagonal matrix, so they fit the
#' Vu slot of \code{\link{RRsol-class}}. The blocks are in the order of the
#' terms, and the rows and columns are named "term:trait".
#'
#' @param Vu a list of covariance matrices, one per term
#' @param terms the names of the terms
#' @param traits the names of the traits
#'
#' @keywords internal
blockVu = function(Vu, terms, traits){
  nT = length(traits)
  out = matrix(0, length(Vu)*nT, length(Vu)*nT)
  for(i in seq_along(Vu)){
    take = (i-1)*nT + seq_len(nT)
    out[take, take] = as.matrix(Vu[[i]])
  }
  labels = paste0(rep(terms, each=nT), ":", rep(traits, length(Vu)))
  dimnames(out) = list(labels, labels)
  return(out)
}

#' @title Error variance of a model fit
#'
#' @description Ve as a matrix, with the traits' names on its rows and
#' columns when there is more than one trait
#'
#' @param Ve the error variance or covariance matrix
#' @param traits the names of the traits
#'
#' @keywords internal
traitMatrix = function(Ve, traits){
  Ve = as.matrix(Ve)
  if( (length(traits)>1) & (nrow(Ve)==length(traits)) ){
    dimnames(Ve) = list(traits, traits)
  }
  return(Ve)
}

#' @title Names of the responses
#'
#' @description The names given to a fit's traits. They are the names of the
#' traits asked for when there is one per column of the responses. A
#' response that is not one column per trait, such as random values or a
#' function of the traits returning several columns, is named "Response1",
#' "Response2" and so on, unless it is a single column, which keeps the name
#' it has always had.
#'
#' @param traits the names of the traits asked for
#' @param nResponse the number of columns of the responses
#'
#' @keywords internal
responseNames = function(traits, nResponse){
  if( (nResponse==1) | (length(traits)==nResponse) ){
    return(traits)
  }
  return(paste0("Response", seq_len(nResponse)))
}

#' @title Fast RR-BLUP fit
#'
#' @description The fast mode of \code{\link{RRBLUP}}, which solves the
#' mixed model equations by preconditioned conjugate gradient. See
#' \code{\link{RRBLUP}} for what it fits and for its arguments. simParam and
#' nThreads have already been resolved by the caller.
#'
#' @param pop the training population
#' @param y the response, as getResponse returns it
#' @param traits the trait or traits to model, as given to RRBLUP
#' @param snpChip the SNP chip, or trait whose QTL, to use
#' @param useQtl whether QTL genotypes are used instead of a SNP chip
#' @param maxIter maximum number of conjugate gradient iterations
#' @param Vu marker effect variance or covariance matrix, or NULL
#' @param Ve error variance or covariance matrix, or NULL
#' @param nSubsample the number of records used to estimate the variance
#' components
#' @param simParam an object of class \code{\link{SimParam}}
#' @param nThreads number of threads
#'
#' @keywords internal
fitFastRRBLUP = function(pop, y, traits, snpChip, useQtl, maxIter,
                         Vu, Ve, nSubsample, simParam, nThreads){
  # A function of the traits can return a plain vector, which has no columns
  # to count
  y = as.matrix(y)

  traits = convertTraitsToNames(traits, simParam)
  nTraits = ncol(y)
  traits = responseNames(traits, nTraits)
  if(nTraits>1){
    checkCompleteRecords(y)
  }

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
    if(nTraits>1){
      Vu = diag(nTraits)
      Ve = diag(nTraits)
    }else{
      Vu = 1
      Ve = 1
    }
    subset = fastSubset(nSubsample, nrow(y))
  }else{
    subset = 1L
  }

  #Fit model
  if(nTraits>1){
    Vu = as.matrix(Vu)
    Ve = as.matrix(Ve)
    if(!all(dim(Vu)==nTraits) | !all(dim(Ve)==nTraits)){
      stop("Vu and Ve must be square matrices with one row per trait")
    }
    ans = callFastRRBLUP_MV(y,fixEff,pop,lociPerChr,
                            lociLoc,Vu,Ve,maxIter,
                            estVarComp,as.integer(subset),
                            nThreads)
  }else{
    ans = callFastRRBLUP(y,fixEff,pop,lociPerChr,
                         lociLoc,Vu,Ve,maxIter,
                         estVarComp,as.integer(subset),
                         nThreads)
  }
  alpha = as.matrix(ans$alpha)

  bv = gv = vector("list",nTraits)

  for(i in seq_len(nTraits)){
    bv[[i]] = new("TraitA",
                  nLoci=nLoci,
                  lociPerChr=lociPerChr,
                  lociLoc=lociLoc,
                  addEff=alpha[,i],
                  intercept=c(ans$beta)[i],
                  name=paste0("est_BV_",traits[i]))

    gv[[i]] = new("TraitA",
                  nLoci=nLoci,
                  lociPerChr=lociPerChr,
                  lociLoc=lociLoc,
                  addEff=alpha[,i],
                  intercept=c(ans$mu)[i],
                  name=paste0("est_GV_",traits[i]))
  }

  Vu = as.matrix(ans$Vu)
  Ve = as.matrix(ans$Ve)
  if(nTraits>1){
    dimnames(Vu) = dimnames(Ve) = list(traits, traits)
  }

  output = new("RRsol",
               bv = bv,
               gv = gv,
               female = as.list(NULL),
               male = as.list(NULL),
               Vu = Vu,
               Ve = Ve,
               model = "RRBLUP")

  return(output)
}



#' @title Records for fast mode variance components
#'
#' @description The records the fast mode of the RR-BLUP models estimates
#' variance components on: a random subset of nSubsample records, or every
#' record when nSubsample is zero or at least the number of records.
#'
#' @param nSubsample the number of records to use
#' @param nRecords the number of records
#'
#' @keywords internal
fastSubset = function(nSubsample, nRecords){
  nSubsample = as.integer(nSubsample)
  if((nSubsample>0L) & (nSubsample<nRecords)){
    # Sampling here rather than in C++ keeps it under set.seed()
    return(sort(sample.int(nRecords, nSubsample)))
  }
  return(seq_len(nRecords))
}

#' @title Variance components from an earlier fit
#'
#' @description Takes Vu and Ve from a fit passed to a model as varComp,
#' after checking that it was fitted with the same model, to the same number
#' of traits and on the same markers. Vu is the variance of a marker effect,
#' so it carries over only to the markers it was estimated on. A fit that
#' records no model, as one created directly or saved before the model slot
#' existed, cannot be checked and is refused.
#'
#' @param varComp the earlier fit, an object of \code{\link{RRsol-class}}
#' @param Vu the Vu argument of the model, which must be NULL
#' @param Ve the Ve argument of the model, which must be NULL
#' @param model the name of the model being fitted
#' @param nTraits the number of traits being fitted
#' @param lociPerChr the number of markers on each chromosome
#' @param lociLoc the positions of the markers on their chromosomes
#'
#' @keywords internal
varCompFromFit = function(varComp, Vu, Ve, model, nTraits, lociPerChr,
                          lociLoc){
  if(!is.null(Vu) | !is.null(Ve)){
    stop("supply either varComp or Vu and Ve, not both")
  }
  if(!isRRsol(varComp)){
    stop("varComp must be the result of a genomic selection model, ",
         "an object of class RRsol")
  }
  fitModel = if(.hasSlot(varComp, "model")) varComp@model else ""
  if((length(fitModel)!=1L) || is.na(fitModel) || (fitModel=="")){
    stop("varComp does not record the model it was fitted with, so its ",
         "variance components cannot be checked; supply Vu and Ve instead")
  }
  if(fitModel!=model){
    stop("varComp was fitted with ", fitModel, ", not ", model)
  }
  if(nrow(varComp@Ve)!=nTraits){
    stop("varComp was fitted to ", nrow(varComp@Ve), " trait(s), not ",
         nTraits)
  }
  fitTrait = varComp@gv[[1]]
  if(!identical(as.integer(fitTrait@lociPerChr), as.integer(lociPerChr)) ||
     !identical(as.integer(fitTrait@lociLoc), as.integer(lociLoc))){
    stop("varComp was fitted on different markers, and its Vu is the ",
         "variance of an effect of those markers")
  }
  return(list(Vu=varComp@Vu, Ve=varComp@Ve))
}

#' @title Fixed variance components outside fast mode
#'
#' @description Takes the variance components supplied to a model whose
#' fastMode is FALSE, as Vu and Ve or as an earlier fit in varComp, and puts
#' them in the form the C++ solvers hold fixed: one variance per random term
#' followed by the error variance. Returns an empty vector when none are
#' supplied, in which case they are estimated by REML as usual. Holding them
#' fixed is for one trait, so supplied values are refused for several, and
#' both Vu and Ve must be given, because estimating one of them alone would
#' ignore the other without notice.
#'
#' @param Vu the Vu argument of the model
#' @param Ve the Ve argument of the model
#' @param varComp the varComp argument of the model
#' @param model the name of the model being fitted
#' @param nTerm the number of random terms of the model
#' @param nTraits the number of traits being fitted
#' @param lociPerChr the number of markers on each chromosome
#' @param lociLoc the positions of the markers on their chromosomes
#'
#' @keywords internal
directVarComp = function(Vu, Ve, varComp, model, nTerm, nTraits, lociPerChr,
                         lociLoc){
  if(is.null(Vu) & is.null(Ve) & is.null(varComp)){
    return(numeric(0))
  }
  if(nTraits>1){
    stop("Vu, Ve and varComp can be supplied for several traits only when ",
         "fastMode is TRUE")
  }
  if(!is.null(varComp)){
    vc = varCompFromFit(varComp, Vu, Ve, model, nTraits, lociPerChr,
                        lociLoc)
    Vu = vc$Vu
    Ve = vc$Ve
  }
  if(is.null(Vu) | is.null(Ve)){
    stop("supply both Vu and Ve, or neither to estimate them")
  }
  Vu = unlist(fastMKVarComp(Vu, Ve, nTerm, 1L)$Vu)
  Ve = as.numeric(Ve)
  if(length(Ve)!=1L){
    stop("Ve must be a single value for one trait")
  }
  # A term with no variance has no effects, as REML often estimates for a
  # dominance or SCA term, and the solver leaves it out
  if(any(!is.finite(c(Vu, Ve))) | any(Vu<0) | !(Ve>0)){
    stop("Vu must not be negative and Ve must be positive")
  }
  return(c(Vu, Ve))
}

#' @title Deprecation message of a numbered model
#'
#' @description The message given when one of the numbered genomic selection
#' models is called, naming the unnumbered model that replaces it and how its
#' variance components are passed there.
#'
#' @param old the name of the numbered model
#' @param new the name of the model that replaces it
#' @param vu how the numbered model's marker variances are written as Vu
#'
#' @keywords internal
numberedDeprecation = function(old, new, vu){
  return(paste0(
    old, " is deprecated. To hold the variance components fixed, pass an ",
    "earlier fit to ", new, " as varComp, or call ", new, " with Vu=", vu,
    " and Ve=Ve. To estimate them, call ", new, " without them. See help(\"",
    old, "\")."))
}

#' @title Variance components for a fast multikernel model
#'
#' @description Puts Vu and Ve in the form the fast mode of
#' \code{\link{RRBLUP_D}}, \code{\link{RRBLUP_GCA}} and
#' \code{\link{RRBLUP_SCA}} passes to C++: a list with one nTraits by nTraits
#' matrix per random term, and an nTraits by nTraits matrix. For one trait,
#' Vu may be a vector with one value per term, and for several traits a list
#' of covariance matrices or the block diagonal matrix held in the Vu slot of
#' a fit, so a fit's own values can be passed back. When either is NULL, both
#' are estimated, and placeholders are returned.
#'
#' @param Vu the marker effect variances, or NULL
#' @param Ve the error variance, or NULL
#' @param nTerm the number of random terms of the model
#' @param nTraits the number of traits fitted
#'
#' @keywords internal
fastMKVarComp = function(Vu, Ve, nTerm, nTraits){
  if(is.null(Vu) | is.null(Ve)){
    # Values are not used, but something has to be passed
    return(list(est=TRUE, Vu=rep(list(diag(nTraits)), nTerm),
                Ve=diag(nTraits)))
  }
  if(is.list(Vu)){
    if(length(Vu)!=nTerm){
      stop("Vu must have one entry per random term, ", nTerm,
           " for this model")
    }
    Vu = lapply(Vu, as.matrix)
  }else if(nTraits==1){
    Vu = as.numeric(Vu)
    if(length(Vu)!=nTerm){
      stop("Vu must have one value per random term, ", nTerm,
           " for this model")
    }
    Vu = lapply(Vu, as.matrix)
  }else{
    Vu = as.matrix(Vu)
    n = nTerm*nTraits
    if(!all(dim(Vu)==n)){
      stop("Vu must be a list of ", nTerm, " covariance matrices, or a ",
           "block diagonal matrix with ", n, " rows, as in the Vu slot of ",
           "a fit")
    }
    blocks = lapply(seq_len(nTerm),
                    function(i) (i-1L)*nTraits+seq_len(nTraits))
    offBlock = Vu
    for(b in blocks){
      offBlock[b,b] = 0
    }
    if(any(offBlock!=0)){
      stop("Vu given as a matrix must be block diagonal, as in the Vu slot ",
           "of a fit")
    }
    Vu = lapply(blocks, function(b) Vu[b,b,drop=FALSE])
  }
  Vu = lapply(Vu, unname)
  return(list(est=FALSE, Vu=Vu, Ve=unname(as.matrix(Ve))))
}

#' @title Costs of the two routes of a genomic selection model
#'
#' @description Estimates of the work done by the direct route of a model
#' and by its fast route, in floating point operations, for deciding which is
#' faster. Only the terms that grow fastest are counted, so the estimates are
#' rough, but they only have to rank the two routes.
#'
#' The direct route of the multikernel models is a REML fit on every record,
#' costed as useCoefficientSide in src/mme.h costs it, on whichever of the
#' variance matrix of the records and Henderson's equations is cheaper. The
#' fast route is the same fit on nSubsample records, given each term's codes
#' or, with more markers than records in the subset, a factor of their cross
#' product, plus the conjugate gradient solve on every record. The direct
#' route of RRBLUP forms and decomposes the smaller of the two cross products
#' of the markers, and its fast route does the same on the subset.
#'
#' The cost of the solve was measured against the REML fits it is compared
#' with, at 20,000 records and 2,000 markers: about 2,500 operations per
#' record, marker and term, growing with about the square root of the number
#' of traits, because the traits share the reads of the genotypes.
#'
#' @param model the name of the model's function
#' @param nInd the number of records
#' @param nMarker the number of markers
#' @param nTraits the number of traits
#' @param nFixEff the number of fixed effect levels
#' @param nSubsample the number of records in the subset
#'
#' @keywords internal
fastModeCosts = function(model, nInd, nMarker, nTraits, nFixEff, nSubsample){
  n = as.double(nInd)
  m = as.double(nMarker)
  nT = as.double(nTraits)
  q = as.double(nFixEff)
  nSub = as.double(nSubsample)
  remlCost = function(nRec, nCoef, nPar){
    nEval = 10
    recordSide = nEval*(nRec^3+nRec^2*nCoef)
    coefficientSide = nRec*nCoef^2+nEval*(nCoef^3+nRec*nCoef*(nPar+2))
    return(min(recordSide, coefficientSide))
  }
  cgCost = function(nKernel){
    return(2500*n*m*nKernel*sqrt(nT))
  }
  if(model=="RRBLUP"){
    s = min(n, m)
    sSub = min(nSub, m)
    # Forming the cross product, and its eigendecomposition with vectors
    direct = n*m*s+10*s^3
    fast = nSub*m*sSub+10*sSub^3+cgCost(1)
  }else{
    nKernel = if(model=="RRBLUP_SCA") 3 else 2
    if(model!="RRBLUP_GCA"){
      # The summed dominance codes
      q = q+1
    }
    nPar = (nKernel+1)*nT*(nT+1)/2
    direct = remlCost(n*nT, nT*(q+nKernel*m), nPar)
    mSub = min(m, nSub)
    fast = remlCost(nSub*nT, nT*(q+nKernel*mSub), nPar)+cgCost(nKernel)
    if(m>=nSub){
      # Each term's cross product on the subset, which the REML fit is given
      # in place of the codes
      fast = fast+nKernel*nSub^2*m
    }
  }
  return(list(direct=direct, fast=fast))
}

#' @title Whether a model takes its fast route
#'
#' @description fastMode asks for the fast route of a model, which is taken
#' when it is expected to be faster than the direct route. Supplied Vu and Ve
#' always take it, because the direct route cannot use them. A subset of
#' every record never does, because the fast route would then repeat the
#' direct route's REML fit and add the solve to it. Otherwise the fast route
#' is taken when its cost from \code{\link{fastModeCosts}} is less than half
#' that of the direct route. The margin allows for the REML fit on the subset,
#' whose number of iterations varies, taking longer than its estimate, as it
#' did by up to half again in the measurements. Setting the
#' option AlphaSimR.fastModeRoute to "fast" always takes the fast route,
#' which the tests use to exercise it on small populations.
#'
#' @param model the name of the model's function
#' @param Vu the supplied marker effect variances, or NULL
#' @param Ve the supplied error variance, or NULL
#' @param nInd the number of records
#' @param nMarker the number of markers
#' @param nTraits the number of traits
#' @param nFixEff the number of fixed effect levels
#' @param nSubsample the number of records in the subset
#'
#' @keywords internal
useFastRoute = function(model, Vu, Ve, nInd, nMarker, nTraits, nFixEff,
                        nSubsample){
  if(!is.null(Vu) & !is.null(Ve)){
    return(TRUE)
  }
  if(identical(getOption("AlphaSimR.fastModeRoute"), "fast")){
    return(TRUE)
  }
  nSubsample = as.integer(nSubsample)
  if((nSubsample<=0L) | (nSubsample>=nInd)){
    return(FALSE)
  }
  costs = fastModeCosts(model, nInd, nMarker, nTraits, nFixEff, nSubsample)
  return(costs$fast<0.5*costs$direct)
}

#' @title RR-BLUP Model
#'
#' @description
#' Fits an additive ridge regression BLUP (RR-BLUP) model for genomic
#' prediction. One trait is fitted with a single trait model and several
#' traits with a multiple trait model; see details for both.
#'
#' @details
#' \subsection{The model}{
#' The model is the ridge regression model of
#' \insertCite{whittaker2000Marker;textual}{AlphaSimR}, applied to
#' genome-wide markers by
#' \insertCite{meuwissen2001Prediction;textual}{AlphaSimR},
#' \deqn{y = X\beta + Mu + e,}
#' where y holds the records, X is the design matrix of the fixed effects,
#' formed from the levels of the population's fixEff slot, M holds the
#' marker dosages, u the marker effects and e the errors. The marker effects
#' are the BLUPs \insertCite{henderson1975Best}{AlphaSimR} given variance
#' components estimated by REML
#' \insertCite{patterson1971Recovery}{AlphaSimR}. The estimated breeding
#' values and genetic values use the same marker effects and differ only in
#' their intercepts.
#' }
#'
#' \subsection{One trait}{
#' When traits gives a single trait, or is a function returning a single
#' value, \eqn{u \sim N(0, I\sigma^2_u)} and \eqn{e \sim N(0, I\sigma^2_e)}.
#' REML is reduced to a one dimensional search over the ratio of the two
#' variances using the spectral approach of FaST-LMM
#' \insertCite{lippert2011FaST}{AlphaSimR}. One eigendecomposition, of
#' whichever of the marker and record cross products is smaller, supplies
#' the likelihood, the fixed effects and the marker effects, so the model is
#' fitted without iteration and maxIter is not used.
#' }
#'
#' \subsection{Several traits}{
#' When traits gives more than one trait, the traits are fitted together as a
#' multiple trait model. The effects of each marker on the traits have an
#' unstructured covariance matrix \eqn{V_u}, the errors of each record have an
#' unstructured covariance matrix \eqn{V_e}, and markers and records are
#' independent of each other. These matrices are estimated by REML using the
#' algorithm of \insertCite{zhou2014Efficient;textual}{AlphaSimR}, which
#' rotates the records by an eigendecomposition of the markers and the traits
#' by a transformation that diagonalizes \eqn{V_u} and \eqn{V_e} together, so
#' that each evaluation of the likelihood costs time linear in the number of
#' records. Starting values come from single trait fits, and the estimates
#' are refined by EM steps followed by Newton steps using the average
#' information matrix \insertCite{gilmour1995Average}{AlphaSimR}, for at most
#' maxIter iterations. The Vu and Ve slots of the result hold the two
#' covariance matrices.
#' }
#'
#' \subsection{Fast mode}{
#' Setting fastMode to TRUE gives a fast and low memory alternative for
#' training populations too large to fit directly. It fits the same model and
#' is used only where it is expected to be faster, see fastMode. The mixed
#' model equations are then solved by preconditioned conjugate gradient,
#' iterating over the genotypes rather than over a stored coefficient matrix
#' \insertCite{stranden1999Solving}{AlphaSimR}. Neither the coefficient
#' matrix nor a numeric copy of the genotypes is formed, so memory use stays
#' close to the size of the genotypes, and the iteration is shared across
#' nThreads.
#'
#' For several traits, every trait must be recorded on every individual. The
#' traits are separated by a canonical transformation, which is exact when
#' every trait is recorded and they share their fixed effects, and the
#' transformed traits are solved together by conjugate gradient, each with
#' its own recurrence and test of convergence. Each iteration costs about as
#' much as one iteration for each trait on its own, so a fit to several
#' traits takes about as long as fitting them one at a time, while the memory
#' used stays close to the size of the genotypes.
#'
#' Variance components are estimated by REML when they are not supplied,
#' with the same method as outside fast mode. REML needs a matrix
#' decomposition whose cost is the square of the number of records in memory
#' and their cube in time, which is what fast mode is built to avoid, so the
#' estimate is taken from a random subset of the records set by nSubsample.
#' Variance components are nuisance parameters for genomic prediction and are
#' estimated far more precisely than a breeding program needs, so a few
#' thousand records give values good enough to shrink with. The marker
#' effects are always solved for using every record.
#' }
#'
#' \subsection{Known variance components}{
#' Retraining a model whose variance components are already known, for
#' example from an earlier fit in the same simulation, needs only the solve
#' for the marker effects. Pass the earlier fit as varComp, or give Vu and Ve.
#' With fastMode TRUE the solve is by fast mode, for any number of traits.
#' With fastMode FALSE the mixed model equations of one trait are solved
#' directly, written on whichever of the markers and the records is the
#' smaller, so the result is exact rather than to the tolerance of an
#' iterative solver. This takes the place of \code{\link{RRBLUP2}} with
#' useEM FALSE, which is deprecated.
#' }
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait or traits to model, a vector of trait names,
#' or a function of the traits returning a single value. More than one trait
#' gives a multiple trait model. See details.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. In fast mode, it limits the
#' iterations of the conjugate gradient solver. Otherwise it is only used
#' when the number of traits is greater than 1, and limits the iterations of
#' the REML fit of the multiple trait model.
#' @param fastMode should the model be solved by preconditioned conjugate
#' gradient where that is expected to be faster, as it is for large training
#' populations. The cost of both ways of fitting the model is estimated from
#' the numbers of records, markers, traits and nSubsample, and fast mode is
#' used only when its cost is less than half that of the default, and never
#' when nSubsample covers every record. Otherwise the model is fitted as it
#' is when fastMode is FALSE, with the same result. Supplying Vu and Ve
#' with fastMode TRUE always uses fast mode. See details.
#' @param Vu marker effect variance, or for several traits the covariance
#' matrix of the marker effects among traits. If value is NULL, this and Ve
#' are estimated by REML. When supplied with Ve, the marker effects are
#' solved for with the variance components held fixed: by fast mode when
#' fastMode is TRUE, and directly otherwise, which is for one trait only. See
#' details.
#' @param Ve error variance, or for several traits the covariance matrix of
#' the errors among traits, which must be positive definite. If value is
#' NULL, this and Vu are estimated by REML. See Vu and details.
#' @param nSubsample the number of records used to estimate the variance
#' components when fastMode is TRUE. Ignored when Vu and Ve are supplied, and
#' when fastMode is FALSE. A value of zero uses every record, which is only
#' advisable for small training populations, because the memory needed grows
#' with the square of this value.
#' @param varComp an earlier fit of this model, from which Vu and Ve are
#' taken in place of supplying them, so that only the marker effects are
#' solved for. Not used together with Vu and Ve. The fit must be of RRBLUP
#' or RRBLUP2, to the same number of traits, on the same markers. See Vu and
#' the genomic selection vignette.
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
#' ans = RRBLUP(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' #Run GS model in fast mode
#' ans = RRBLUP(pop, fastMode=TRUE, simParam=SP)
#'
#' #Refit in fast mode, taking the variance components from that fit
#' ans = RRBLUP(pop, fastMode=TRUE, varComp=ans, simParam=SP)
#'
#' #Refit directly with the same variance components held fixed
#' ans = RRBLUP(pop, varComp=ans, simParam=SP)
#'
#' @export
RRBLUP = function(pop, traits=1, use="pheno", snpChip=1,
                  useQtl=FALSE, maxIter=1000L,
                  simParam=NULL, nThreads=NULL, ...,
                  fastMode=FALSE, Vu=NULL, Ve=NULL, nSubsample=3000L,
                  varComp=NULL){
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  # Variance components supplied outside fast mode are held fixed by the
  # direct solver below. In fast mode they are passed to the fast route,
  # which supplying them always takes.
  fixedVc = numeric(0)
  if(fastMode){
    if(useQtl){
      chip = simParam$traits[[snpChip]]
    }else{
      chip = simParam$snpChips[[snpChip]]
    }
    nLoci = chip@nLoci
    if(!is.null(varComp)){
      vc = varCompFromFit(varComp, Vu, Ve, "RRBLUP", NCOL(y),
                          chip@lociPerChr, chip@lociLoc)
      Vu = vc$Vu
      Ve = vc$Ve
    }
    if(useFastRoute("RRBLUP", Vu, Ve, NROW(y), nLoci, NCOL(y),
                    length(unique(pop@fixEff)), nSubsample)){
      return(fitFastRRBLUP(pop=pop, y=y, traits=traits, snpChip=snpChip,
                           useQtl=useQtl, maxIter=maxIter, Vu=Vu, Ve=Ve,
                           nSubsample=nSubsample, simParam=simParam,
                           nThreads=nThreads))
    }
  }
  # A function of the traits can return a plain vector, which has no columns
  # to count
  y = as.matrix(y)

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
  if(!fastMode){
    fixedVc = directVarComp(Vu, Ve, varComp, "RRBLUP", 1L, ncol(y),
                            lociPerChr, lociLoc)
  }

  #Fit model
  if(ncol(y)>1){
    ans = callRRBLUP_MV(y, fixEff, pop, lociPerChr,
                        lociLoc, maxIter, nThreads)
  }else{
    ans = callRRBLUP(y, fixEff, pop, lociPerChr, lociLoc,
                     nThreads, fixedVc)
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
               Ve = as.matrix(ans$Ve),
               model = "RRBLUP")

  return(output)
}

#' @title RR-BLUP Model 2
#'
#' @description
#' This function is deprecated. Use \code{\link{RRBLUP}}, which takes known
#' variance components as Vu and Ve or as an earlier fit in varComp.
#'
#' Fits the single trait model of \code{\link{RRBLUP}} with variance
#' components estimated by the EM algorithm, or taken as given. See details.
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
#' \subsection{The model}{
#' The model is the single trait ridge regression model of
#' \insertCite{whittaker2000Marker;textual}{AlphaSimR}, applied to
#' genome-wide markers by
#' \insertCite{meuwissen2001Prediction;textual}{AlphaSimR},
#' \deqn{y = X\beta + Mu + e,}
#' where y holds the records, X is the design matrix of the fixed effects,
#' formed from the levels of the population's fixEff slot, M holds the
#' marker dosages, \eqn{u \sim N(0, I\sigma^2_u)} the marker effects and
#' \eqn{e \sim N(0, I\sigma^2_e)} the errors. The marker effects are the BLUPs
#' \insertCite{henderson1975Best}{AlphaSimR} given the variance components.
#' When useEM is TRUE, the variance components are estimated by EM-REML
#' \insertCite{dempster1977Maximum,patterson1971Recovery}{AlphaSimR},
#' starting from Vu and Ve and iterating until the change is below tol or
#' maxIter is reached. When useEM is FALSE, Vu and Ve are taken as the true
#' values. Only one trait can be modeled.
#' }
#'
#' \subsection{Replacing it}{
#' RRBLUP2 was for retraining a model cheaply when the variance components
#' are already known. \code{\link{RRBLUP}} now does the same: called with Vu
#' and Ve, or with an earlier fit as varComp, and fastMode FALSE, it solves
#' the mixed model equations directly and gives the same result, to
#' rounding, as RRBLUP2 with useEM FALSE. With fastMode TRUE it does so for
#' several traits as well.
#'
#' When the variance components are unknown, use \code{\link{RRBLUP}}
#' without them. RRBLUP2 estimates them with the EM algorithm, which is
#' generally far slower than the spectral REML approach of
#' \code{\link{RRBLUP}}, and its default maxIter is typically too small for
#' convergence. When the algorithm fails to converge a warning is displayed
#' and the last iteration is returned.
#' }
#'
#' \subsection{How the model is solved}{
#' There are two ways to solve this model with the variance components
#' fixed, and RRBLUP2, like \code{\link{RRBLUP}} given Vu and Ve, uses
#' whichever is cheaper.
#'
#' Henderson's mixed model equations work with a square matrix whose
#' dimensions are the number of fixed effects plus the number of markers.
#' Written on the records instead, the model works with a square matrix whose
#' dimensions are the number of records, which is the same shape the EMMA
#' approach of \code{\link{RRBLUP}} works with. The two give the same answer.
#' The cost of each is the cost of the other with the number of records and
#' the number of markers exchanged, so the smaller of those two numbers names
#' the cheaper method, and the crossover sits where they are equal. Memory
#' crosses over at the same place.
#'
#' Estimating variance components needs a quantity that only Henderson's
#' equations supply, so useEM = TRUE always uses the marker side, and there
#' the old advice holds: it works best when the markers do not greatly
#' outnumber the records. With useEM = FALSE the choice is free and the
#' cheaper method is taken, so more markers than records is no longer a
#' problem. \code{\link{RRBLUPMemUse}} takes a useEM argument for this reason
#' and will tell you what either choice costs.
#' }
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
#' #Run GS model, then refit with its variance components held fixed,
#' #which replaces RRBLUP2 with useEM=FALSE
#' ans = RRBLUP(pop, simParam=SP)
#' ans = RRBLUP(pop, Vu=ans@Vu, Ve=ans@Ve, simParam=SP)
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
  .Deprecated(msg=numberedDeprecation("RRBLUP2", "RRBLUP", "Vu"))
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  # varA and varE are unnamed and indexed by position, so the numeric
  # index is taken before traits is replaced by its name.
  traitIndex = convertTraitsToIndex(traits, simParam)
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
  if(is.null(traitIndex)){
    if(is.null(Vu)){
      Vu = var(y)/nLoci
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot("this model fits a single trait, so traits must give one trait"=
                length(traits)==1)
    if(is.null(Vu)){
      Vu = 2*simParam$varA[traitIndex]/nLoci
      if(is.na(Vu)){
        Vu = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traitIndex]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }

  #Fit model
  ans = callRRBLUP2(y, fixEff, pop, lociPerChr,
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
               Ve = as.matrix(ans$Ve),
               model = "RRBLUP")

  return(output)
}

#' @title RR-BLUP Model with Dominance
#'
#' @description
#' Fits an RR-BLUP model for genomic prediction that includes dominance
#' effects. One trait is fitted with a single trait model and several traits
#' with a multiple trait model; see details for both.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait or traits to model, a vector
#' of trait names, or a function of the traits returning a single value.
#' Several traits are fitted together as a multiple trait model, which needs
#' every trait recorded on every individual. See details.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. In fast mode, it limits the
#' iterations of the conjugate gradient solver, and NULL gives 1000.
#' Otherwise it limits the iterations of the REML fit, and NULL gives 100.
#' @param fastMode should the model be solved by preconditioned conjugate
#' gradient where that is expected to be faster, as it is for large training
#' populations. The cost of both ways of fitting the model is estimated from
#' the numbers of records, markers, traits and nSubsample, and fast mode is
#' used only when its cost is less than half that of the default, and never
#' when nSubsample covers every record. Otherwise the model is fitted as it
#' is when fastMode is FALSE, with the same result. Supplying Vu and Ve
#' with fastMode TRUE always uses fast mode. See details.
#' @param Vu the variances of the marker effects, one per random term. For
#' one trait, a vector, and for several traits, a list of covariance
#' matrices among the traits, one per term, or the block diagonal matrix held
#' in the Vu slot of a fit. If value is NULL, this and Ve are estimated by
#' REML. When supplied with Ve, the marker effects are solved for with the
#' variance components held fixed: by fast mode when fastMode is TRUE, and
#' directly otherwise, which is for one trait only. See details.
#' @param Ve the error variance, or for several traits the covariance matrix
#' of the errors among traits, which must be positive definite. If value is
#' NULL, this and Vu are estimated by REML. See Vu and details.
#' @param nSubsample the number of records used to estimate the variance
#' components when fastMode is TRUE. Ignored when Vu and Ve are supplied, and
#' when fastMode is FALSE. A value of zero uses every record, which is only
#' advisable for small training populations.
#' @param varComp an earlier fit of this model, from which Vu and Ve are
#' taken in place of supplying them, so that only the marker effects are
#' solved for. Not used together with Vu and Ve. The fit must be of this
#' model or of its counterpart ending in 2, to the same number of traits, on
#' the same markers. See Vu and the genomic selection vignette.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @details
#' \subsection{The model}{
#' The model is
#' \deqn{y = X\beta + hb + M_au_a + M_du_d + e,}
#' where y holds the records, X is the design matrix of the fixed effects,
#' formed from the levels of the population's fixEff slot, \eqn{M_a} holds
#' the additive marker codings and \eqn{M_d} the dominance codings, with
#' \eqn{u_a} and \eqn{u_d} their effects, and e the errors. h is each
#' individual's sum of dominance codings, a measure of its genome-wide
#' heterozygosity, and its fixed regression b captures directional dominance
#' such as inbreeding depression \insertCite{xiang2016Genomic}{AlphaSimR}.
#' The marker effects are the BLUPs \insertCite{henderson1975Best}{AlphaSimR}
#' given variance components estimated by REML
#' \insertCite{patterson1971Recovery}{AlphaSimR}. The estimated genetic
#' values use the additive and dominance effects. The estimated breeding
#' values use the average effects of allele substitution
#' \insertCite{falconer1996Introduction}{AlphaSimR}, which are worked out
#' from those effects and the genotype frequencies in the training
#' population.
#' }
#'
#' \subsection{One trait}{
#' When traits gives a single trait, or is a function returning a single
#' value, \eqn{u_a \sim N(0, I\sigma^2_a)}, \eqn{u_d \sim N(0, I\sigma^2_d)}
#' and \eqn{e \sim N(0, I\sigma^2_e)}. The three variances are estimated by
#' average information REML
#' \insertCite{gilmour1995Average,johnson1995Restricted}{AlphaSimR}, with
#' EM-REML steps taken where an average information step fails to increase
#' the likelihood, for at most maxIter iterations. The model is fitted in
#' whichever of two equivalent forms is cheaper: on the variance matrix of the
#' records when there are fewer records than marker effects, and on
#' Henderson's mixed model equations otherwise.
#' }
#'
#' \subsection{Several traits}{
#' When traits gives more than one trait, the traits are fitted together as a
#' multiple trait model, which needs every trait recorded on every
#' individual. The additive effects of each marker, its dominance effects and
#' the errors of each record each have an unstructured covariance matrix
#' among the traits. These are estimated by REML, with EM steps followed by
#' Newton steps using the average information matrix
#' \insertCite{gilmour1995Average}{AlphaSimR}, for at most maxIter
#' iterations. The Vu slot of the result holds the two covariance matrices
#' of the marker effects as a block diagonal matrix, additive first, with rows
#' and columns named "a:trait" and "d:trait".
#' }
#'
#' \subsection{Fast mode}{
#' Setting fastMode to TRUE gives a fast and low memory alternative for
#' training populations too large to fit directly. It is used only where it
#' is expected to be faster, see fastMode. The mixed model equations
#' of every random term and trait are then solved together by preconditioned
#' conjugate gradient, iterating over the genotypes rather than over a stored
#' coefficient matrix \insertCite{stranden1999Solving}{AlphaSimR}, so memory
#' use stays close to the size of the genotypes. The variance components are
#' estimated by REML on a random subset of the records, set by nSubsample,
#' with the same method as outside fast mode, unless Vu and Ve are supplied,
#' and the marker effects are always solved for using every record. Given the
#' variance components, the result is the BLUP of the same model that is
#' fitted when fastMode is FALSE.
#' }
#'
#' \subsection{Known variance components}{
#' Retraining a model whose variance components are already known,
#' for example from an earlier fit in the same simulation, needs only
#' the solve for the marker effects. Pass the earlier fit as varComp, or
#' give Vu and Ve. With fastMode TRUE the solve is by fast mode, for any
#' number of traits. With fastMode FALSE the mixed model equations of one
#' trait are solved directly, on whichever of the markers and the records
#' is cheaper, with Vu holding one variance per term. A variance of zero,
#' as REML often estimates for a dominance term, leaves its term out. This
#' takes the place of \code{\link{RRBLUP_D2}} with useEM FALSE, which is
#' deprecated, and Vu is its Va and Vd in that order.
#' }
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
#' #Run GS model in fast mode
#' ans = RRBLUP_D(pop, fastMode=TRUE, simParam=SP)
#'
#' @export
RRBLUP_D = function(pop, traits=1, use="pheno", snpChip=1,
                    useQtl=FALSE, maxIter=NULL,
                    simParam=NULL, nThreads=NULL, ...,
                    fastMode=FALSE, Vu=NULL, Ve=NULL, nSubsample=3000L,
                    varComp=NULL){
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)
  # A function of the traits can return a plain vector, which has no columns
  # to count
  y = as.matrix(y)

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
  nTraits = ncol(y)
  traits = responseNames(traits, nTraits)
  if(nTraits>1){
    checkCompleteRecords(y)
  }
  # Variance components supplied outside fast mode are held fixed by the
  # direct solver below. In fast mode they are passed to the fast route,
  # which supplying them always takes.
  fixedVc = numeric(0)
  if(!fastMode){
    fixedVc = directVarComp(Vu, Ve, varComp, "RRBLUP_D", 2L, nTraits,
                            lociPerChr, lociLoc)
  }
  if(fastMode){
    if(!is.null(varComp)){
      vc = varCompFromFit(varComp, Vu, Ve, "RRBLUP_D", nTraits, lociPerChr,
                          lociLoc)
      Vu = vc$Vu
      Ve = vc$Ve
    }
    fastMode = useFastRoute("RRBLUP_D", Vu, Ve, nrow(y), nLoci, nTraits,
                            length(unique(pop@fixEff)), nSubsample)
  }
  if(is.null(maxIter)){
    maxIter = if(fastMode) 1000L else 100L
  }
  if(fastMode){
    vc = fastMKVarComp(Vu, Ve, 2L, nTraits)
    subset = if(vc$est) fastSubset(nSubsample, nrow(y)) else 1L
    ans = callFastRRBLUP_D(y, fixEff, pop, lociPerChr, lociLoc,
                           vc$Vu, vc$Ve, maxIter, vc$est,
                           as.integer(subset), nThreads)
  }else{
    ans = callRRBLUP_D(y, fixEff, pop, lociPerChr,
                       lociLoc, maxIter, nThreads, fixedVc)
  }

  alpha = as.matrix(ans$alpha)
  a = as.matrix(ans$a)
  d = as.matrix(ans$d)
  bv = gv = vector("list",nTraits)

  for(i in seq_len(nTraits)){
    bv[[i]] = new("TraitA",
                  nLoci=nLoci,
                  lociPerChr=lociPerChr,
                  lociLoc=lociLoc,
                  addEff=alpha[,i],
                  intercept=c(ans$beta)[i],
                  name=paste0("est_BV_",traits[i]))

    gv[[i]] = new("TraitAD",
                  nLoci=nLoci,
                  lociPerChr=lociPerChr,
                  lociLoc=lociLoc,
                  addEff=a[,i],
                  domEff=d[,i],
                  intercept=c(ans$mu)[i],
                  name=paste0("est_GV_",traits[i]))
  }

  if(nTraits>1){
    Vu = blockVu(ans$Vu, c("a","d"), traits)
  }else{
    Vu = as.matrix(ans$Vu)
  }

  output = new("RRsol",
               bv = bv,
               gv = gv,
               female = as.list(NULL),
               male = as.list(NULL),
               Vu = Vu,
               Ve = traitMatrix(ans$Ve, traits),
               model = "RRBLUP_D")

  return(output)
}


#' @title RR-BLUP with Dominance Model 2
#'
#' @description
#' This function is deprecated. Use \code{\link{RRBLUP_D}}, which takes
#' known variance components as Vu and Ve or as an earlier fit in varComp.
#'
#' Fits the single trait model of \code{\link{RRBLUP_D}} with variance
#' components estimated by the EM algorithm, or taken as given. See details.
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
#' @param maxIter maximum number of iterations of the EM algorithm.
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
#' @details
#' \subsection{The model}{
#' The model is the single trait model of \code{\link{RRBLUP_D}},
#' \deqn{y = X\beta + hb + M_au_a + M_du_d + e,}
#' with \eqn{u_a \sim N(0, I\sigma^2_a)}, \eqn{u_d \sim N(0, I\sigma^2_d)}
#' and \eqn{e \sim N(0, I\sigma^2_e)}, where h is each individual's sum of
#' dominance codings, fitted as a covariate for inbreeding depression
#' \insertCite{xiang2016Genomic}{AlphaSimR}. The marker effects are the BLUPs
#' \insertCite{henderson1975Best}{AlphaSimR} given the variance components.
#' When useEM is TRUE, the variance components are estimated by EM-REML
#' \insertCite{dempster1977Maximum,patterson1971Recovery}{AlphaSimR},
#' starting from Va, Vd and Ve and iterating until the change is below tol or
#' maxIter is reached. When useEM is FALSE, Va, Vd and Ve are taken as the
#' true values. Only one trait can be modeled. The estimated breeding values
#' use average effects of allele substitution, as described for
#' \code{\link{RRBLUP_D}}.
#' }
#'
#' \subsection{Replacing it}{
#' This function was for retraining a model cheaply when the variance
#' components are already known. \code{\link{RRBLUP_D}} now does the same:
#' called with an earlier fit as varComp, or with Vu=c(Va, Vd) and Ve,
#' and fastMode FALSE, it solves the mixed model equations directly with
#' the solver of its own REML fit and gives the same result, to rounding,
#' as this function with useEM FALSE. With fastMode TRUE it does so for
#' several traits as well. When the variance components are unknown, use
#' \code{\link{RRBLUP_D}} without them: the EM algorithm here is generally
#' far slower than its average information REML, and the default maxIter
#' is typically too small for convergence. When the algorithm fails to
#' converge a warning is displayed and the last iteration is returned.
#' }
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
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Run GS model, then refit with its variance components held fixed,
#' #which replaces RRBLUP_D2 with useEM=FALSE
#' ans = RRBLUP_D(pop, simParam=SP)
#' ans = RRBLUP_D(pop, varComp=ans, simParam=SP)
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
  .Deprecated(msg=numberedDeprecation("RRBLUP_D2", "RRBLUP_D", "c(Va, Vd)"))
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  # varA and varE are unnamed and indexed by position, so the numeric
  # index is taken before traits is replaced by its name.
  traitIndex = convertTraitsToIndex(traits, simParam)
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
  if(is.null(traitIndex)){
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
    stopifnot("this model fits a single trait, so traits must give one trait"=
                length(traits)==1)
    if(is.null(Va)){
      Va = 2*simParam$varA[traitIndex]/nLoci
      if(is.na(Va)){
        Va = var(y)/nLoci
      }
    }
    if(is.null(Vd)){
      Vd = 2*simParam$varA[traitIndex]/nLoci
      if(is.na(Vd)){
        Vd = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traitIndex]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }

  #Fit model
  stopifnot("this model fits a single trait, so the response must have one column"=
              ncol(y)==1)
  ans = callRRBLUP_D2(y, fixEff, pop, lociPerChr,
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
               Ve = as.matrix(ans$Ve),
               model = "RRBLUP_D")

  return(output)
}

#' @title RR-BLUP GCA Model
#'
#' @description
#' Fits an RR-BLUP model that estimates separate marker effects for
#' females and males, for predicting the GCA of parents of single cross
#' hybrids. One trait is fitted with a single trait model and several traits
#' with a multiple trait model; see details for both.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait or traits to model, a vector
#' of trait names, or a function of the traits returning a single value.
#' Several traits are fitted together as a multiple trait model, which needs
#' every trait recorded on every individual. See details.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. In fast mode, it limits the
#' iterations of the conjugate gradient solver, and NULL gives 1000.
#' Otherwise it limits the iterations of the REML fit, and NULL gives 100.
#' @param fastMode should the model be solved by preconditioned conjugate
#' gradient where that is expected to be faster, as it is for large training
#' populations. The cost of both ways of fitting the model is estimated from
#' the numbers of records, markers, traits and nSubsample, and fast mode is
#' used only when its cost is less than half that of the default, and never
#' when nSubsample covers every record. Otherwise the model is fitted as it
#' is when fastMode is FALSE, with the same result. Supplying Vu and Ve
#' with fastMode TRUE always uses fast mode. See details.
#' @param Vu the variances of the marker effects, one per random term. For
#' one trait, a vector, and for several traits, a list of covariance
#' matrices among the traits, one per term, or the block diagonal matrix held
#' in the Vu slot of a fit. If value is NULL, this and Ve are estimated by
#' REML. When supplied with Ve, the marker effects are solved for with the
#' variance components held fixed: by fast mode when fastMode is TRUE, and
#' directly otherwise, which is for one trait only. See details.
#' @param Ve the error variance, or for several traits the covariance matrix
#' of the errors among traits, which must be positive definite. If value is
#' NULL, this and Vu are estimated by REML. See Vu and details.
#' @param nSubsample the number of records used to estimate the variance
#' components when fastMode is TRUE. Ignored when Vu and Ve are supplied, and
#' when fastMode is FALSE. A value of zero uses every record, which is only
#' advisable for small training populations.
#' @param varComp an earlier fit of this model, from which Vu and Ve are
#' taken in place of supplying them, so that only the marker effects are
#' solved for. Not used together with Vu and Ve. The fit must be of this
#' model or of its counterpart ending in 2, to the same number of traits, on
#' the same markers. See Vu and the genomic selection vignette.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @details
#' \subsection{The model}{
#' The model is
#' \deqn{y = X\beta + M_fu_f + M_mu_m + e,}
#' where y holds the records, X is the design matrix of the fixed effects,
#' formed from the levels of the population's fixEff slot, \eqn{M_f} and
#' \eqn{M_m} hold the marker dosages an individual inherited from its female
#' and male parents, \eqn{u_f} and \eqn{u_m} their effects, and e the errors.
#' Fitting separate effects for each parent lets an allele be worth something
#' different depending on the heterotic pool it came from
#' \insertCite{massman2013Genomewide,dejong2023Comparison}{AlphaSimR}. The
#' marker effects are the BLUPs \insertCite{henderson1975Best}{AlphaSimR}
#' given variance components estimated by REML
#' \insertCite{patterson1971Recovery}{AlphaSimR}.
#' }
#'
#' \subsection{One trait}{
#' When traits gives a single trait, or is a function returning a single
#' value, \eqn{u_f \sim N(0, I\sigma^2_f)}, \eqn{u_m \sim N(0, I\sigma^2_m)}
#' and \eqn{e \sim N(0, I\sigma^2_e)}. The three variances are estimated by
#' average information REML
#' \insertCite{gilmour1995Average,johnson1995Restricted}{AlphaSimR}, with
#' EM-REML steps taken where an average information step fails to increase
#' the likelihood, for at most maxIter iterations. The model is fitted in
#' whichever of two equivalent forms is cheaper: on the variance matrix of the
#' records when there are fewer records than marker effects, and on
#' Henderson's mixed model equations otherwise.
#' }
#'
#' \subsection{Several traits}{
#' When traits gives more than one trait, the traits are fitted together as a
#' multiple trait model, which needs every trait recorded on every
#' individual. The female effects of each marker, its male effects and the
#' errors of each record each have an unstructured covariance matrix among
#' the traits. These are estimated by REML, with EM steps followed by Newton
#' steps using the average information matrix
#' \insertCite{gilmour1995Average}{AlphaSimR}, for at most maxIter
#' iterations. The Vu slot of the result holds the covariance matrices of the
#' marker effects as a block diagonal matrix, with rows and columns named
#' "female:trait" and "male:trait".
#' }
#'
#' \subsection{Using the solution}{
#' The solution can be applied with \code{\link{setEBV}}:
#' \code{value="female"} or \code{value="male"} gives parental
#' predictions, while \code{value="gv"} predicts total hybrid genetic value.
#' These model predictions obviously differ from the observed cross means
#' returned by \code{\link{calcGCA}}.
#' }
#'
#' \subsection{Fast mode}{
#' Setting fastMode to TRUE gives a fast and low memory alternative for
#' training populations too large to fit directly. It is used only where it
#' is expected to be faster, see fastMode. The mixed model equations
#' of every random term and trait are then solved together by preconditioned
#' conjugate gradient, iterating over the genotypes rather than over a stored
#' coefficient matrix \insertCite{stranden1999Solving}{AlphaSimR}, so memory
#' use stays close to the size of the genotypes. The variance components are
#' estimated by REML on a random subset of the records, set by nSubsample,
#' with the same method as outside fast mode, unless Vu and Ve are supplied,
#' and the marker effects are always solved for using every record. Given the
#' variance components, the result is the BLUP of the same model that is
#' fitted when fastMode is FALSE.
#' }
#'
#' \subsection{Known variance components}{
#' Retraining a model whose variance components are already known,
#' for example from an earlier fit in the same simulation, needs only
#' the solve for the marker effects. Pass the earlier fit as varComp, or
#' give Vu and Ve. With fastMode TRUE the solve is by fast mode, for any
#' number of traits. With fastMode FALSE the mixed model equations of one
#' trait are solved directly, on whichever of the markers and the records
#' is cheaper, with Vu holding one variance per term. A variance of zero,
#' as REML often estimates for a dominance term, leaves its term out. This
#' takes the place of \code{\link{RRBLUP_GCA2}} with useEM FALSE, which is
#' deprecated, and Vu is its VuF and VuM in that order.
#' }
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
#' ans = RRBLUP_GCA(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' #Run GS model in fast mode
#' ans = RRBLUP_GCA(pop, fastMode=TRUE, simParam=SP)
#'
#' @export
RRBLUP_GCA = function(pop, traits=1, use="pheno", snpChip=1,
                      useQtl=FALSE, maxIter=NULL,
                      simParam=NULL, nThreads=NULL, ...,
                      fastMode=FALSE, Vu=NULL, Ve=NULL, nSubsample=3000L,
                      varComp=NULL){
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)
  # A function of the traits can return a plain vector, which has no columns
  # to count
  y = as.matrix(y)

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
  nTraits = ncol(y)
  traits = responseNames(traits, nTraits)
  if(nTraits>1){
    checkCompleteRecords(y)
  }
  # Variance components supplied outside fast mode are held fixed by the
  # direct solver below. In fast mode they are passed to the fast route,
  # which supplying them always takes.
  fixedVc = numeric(0)
  if(!fastMode){
    fixedVc = directVarComp(Vu, Ve, varComp, "RRBLUP_GCA", 2L, nTraits,
                            lociPerChr, lociLoc)
  }
  if(fastMode){
    if(!is.null(varComp)){
      vc = varCompFromFit(varComp, Vu, Ve, "RRBLUP_GCA", nTraits, lociPerChr,
                          lociLoc)
      Vu = vc$Vu
      Ve = vc$Ve
    }
    fastMode = useFastRoute("RRBLUP_GCA", Vu, Ve, nrow(y), nLoci, nTraits,
                            length(unique(pop@fixEff)), nSubsample)
  }
  if(is.null(maxIter)){
    maxIter = if(fastMode) 1000L else 100L
  }
  if(fastMode){
    vc = fastMKVarComp(Vu, Ve, 2L, nTraits)
    subset = if(vc$est) fastSubset(nSubsample, nrow(y)) else 1L
    ans = callFastRRBLUP_GCA(y, fixEff, pop, lociPerChr, lociLoc,
                             vc$Vu, vc$Ve, maxIter, vc$est,
                             as.integer(subset), nThreads)
  }else{
    ans = callRRBLUP_GCA(y, fixEff, pop,
                         lociPerChr, lociLoc, maxIter,
                         nThreads, fixedVc)
  }

  alpha1 = as.matrix(ans$alpha1)
  alpha2 = as.matrix(ans$alpha2)
  gv = female = male = vector("list",nTraits)

  for(i in seq_len(nTraits)){
    gv[[i]] = new("TraitA2",
                  nLoci=nLoci,
                  lociPerChr=lociPerChr,
                  lociLoc=lociLoc,
                  addEff=alpha1[,i],
                  addEffMale=alpha2[,i],
                  intercept=c(ans$mu)[i],
                  name=paste0("est_GV_",traits[i]))

    female[[i]] = new("TraitA",
                      nLoci=nLoci,
                      lociPerChr=lociPerChr,
                      lociLoc=lociLoc,
                      addEff=alpha1[,i],
                      intercept=c(ans$beta1)[i],
                      name=paste0("est_female_",traits[i]))

    male[[i]] = new("TraitA",
                    nLoci=nLoci,
                    lociPerChr=lociPerChr,
                    lociLoc=lociLoc,
                    addEff=alpha2[,i],
                    intercept=c(ans$beta2)[i],
                    name=paste0("est_male_",traits[i]))
  }

  if(nTraits>1){
    Vu = blockVu(ans$Vu, c("female","male"), traits)
  }else{
    Vu = as.matrix(ans$Vu)
  }

  output = new("RRsol",
               gv = gv,
               bv = as.list(NULL),
               female = female,
               male = male,
               Vu = Vu,
               Ve = traitMatrix(ans$Ve, traits),
               model = "RRBLUP_GCA")

  return(output)
}

#' @title RR-BLUP GCA Model 2
#'
#' @description
#' This function is deprecated. Use \code{\link{RRBLUP_GCA}}, which takes
#' known variance components as Vu and Ve or as an earlier fit in varComp.
#'
#' Fits the single trait model of \code{\link{RRBLUP_GCA}} with variance
#' components estimated by the EM algorithm, or taken as given. See details.
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
#' @details
#' \subsection{The model}{
#' The model is the single trait model of \code{\link{RRBLUP_GCA}},
#' \deqn{y = X\beta + M_fu_f + M_mu_m + e,}
#' with \eqn{u_f \sim N(0, I\sigma^2_f)}, \eqn{u_m \sim N(0, I\sigma^2_m)}
#' and \eqn{e \sim N(0, I\sigma^2_e)}, where \eqn{M_f} and \eqn{M_m} hold
#' the marker dosages an individual inherited from its female and male
#' parents. The marker effects are the BLUPs
#' \insertCite{henderson1975Best}{AlphaSimR} given the variance components.
#' When useEM is TRUE, the variance components are estimated by EM-REML
#' \insertCite{dempster1977Maximum,patterson1971Recovery}{AlphaSimR},
#' starting from VuF, VuM and Ve and iterating until the change is below tol
#' or maxIter is reached. When useEM is FALSE, VuF, VuM and Ve are taken as
#' the true values. Only one trait can be modeled.
#' }
#'
#' \subsection{Replacing it}{
#' This function was for retraining a model cheaply when the variance
#' components are already known. \code{\link{RRBLUP_GCA}} now does the same:
#' called with an earlier fit as varComp, or with Vu=c(VuF, VuM) and Ve,
#' and fastMode FALSE, it solves the mixed model equations directly with
#' the solver of its own REML fit and gives the same result, to rounding,
#' as this function with useEM FALSE. With fastMode TRUE it does so for
#' several traits as well. When the variance components are unknown,
#' use \code{\link{RRBLUP_GCA}} without them: the EM algorithm here is
#' generally far slower than its average information REML, and the default
#' maxIter is typically too small for convergence. When the algorithm fails
#' to converge a warning is displayed and the last iteration is returned.
#' }
#'
#' \subsection{Using the solution}{
#' The solution can be applied with \code{\link{setEBV}}:
#' \code{value="female"} or \code{value="male"} gives parental
#' predictions, while \code{value="gv"} predicts total hybrid genetic value.
#' These model predictions obviously differ from the observed cross means
#' returned by \code{\link{calcGCA}}.
#' }
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
#' #Run GS model, then refit with its variance components held fixed,
#' #which replaces RRBLUP_GCA2 with useEM=FALSE
#' ans = RRBLUP_GCA(pop, simParam=SP)
#' ans = RRBLUP_GCA(pop, varComp=ans, simParam=SP)
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
  .Deprecated(msg=numberedDeprecation("RRBLUP_GCA2", "RRBLUP_GCA",
                                      "c(VuF, VuM)"))
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  # varA and varE are unnamed and indexed by position, so the numeric
  # index is taken before traits is replaced by its name.
  traitIndex = convertTraitsToIndex(traits, simParam)
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
  if(is.null(traitIndex)){
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
    stopifnot("this model fits a single trait, so traits must give one trait"=
                length(traits)==1)
    if(is.null(VuF)){
      VuF = 2*simParam$varA[traitIndex]/nLoci
      if(is.na(VuF)){
        VuF = var(y)/nLoci
      }
    }
    if(is.null(VuM)){
      VuM = 2*simParam$varA[traitIndex]/nLoci
      if(is.na(VuM)){
        VuM = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traitIndex]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }

  #Fit model
  stopifnot("this model fits a single trait, so the response must have one column"=
              ncol(y)==1)

  ans = callRRBLUP_GCA2(y, fixEff, pop,
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
               Ve = as.matrix(ans$Ve),
               model = "RRBLUP_GCA")

  return(output)
}

#' @title RR-BLUP SCA Model
#'
#' @description
#' An extension of \code{\link{RRBLUP_GCA}} that adds dominance effects.
#' One trait is fitted with a single trait model and several traits with a
#' multiple trait model; see details for both. Only diploids are supported.
#' Note that we have not seen any consistent benefit of this model over
#' \code{\link{RRBLUP_GCA}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait or traits to model, a vector
#' of trait names, or a function of the traits returning a single value.
#' Several traits are fitted together as a multiple trait model, which needs
#' every trait recorded on every individual. See details.
#' @param use train model using phenotypes "pheno", genetic values "gv",
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip.
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. In fast mode, it limits the
#' iterations of the conjugate gradient solver, and NULL gives 1000.
#' Otherwise it limits the iterations of the REML fit, and NULL gives 100.
#' @param fastMode should the model be solved by preconditioned conjugate
#' gradient where that is expected to be faster, as it is for large training
#' populations. The cost of both ways of fitting the model is estimated from
#' the numbers of records, markers, traits and nSubsample, and fast mode is
#' used only when its cost is less than half that of the default, and never
#' when nSubsample covers every record. Otherwise the model is fitted as it
#' is when fastMode is FALSE, with the same result. Supplying Vu and Ve
#' with fastMode TRUE always uses fast mode. See details.
#' @param Vu the variances of the marker effects, one per random term. For
#' one trait, a vector, and for several traits, a list of covariance
#' matrices among the traits, one per term, or the block diagonal matrix held
#' in the Vu slot of a fit. If value is NULL, this and Ve are estimated by
#' REML. When supplied with Ve, the marker effects are solved for with the
#' variance components held fixed: by fast mode when fastMode is TRUE, and
#' directly otherwise, which is for one trait only. See details.
#' @param Ve the error variance, or for several traits the covariance matrix
#' of the errors among traits, which must be positive definite. If value is
#' NULL, this and Vu are estimated by REML. See Vu and details.
#' @param nSubsample the number of records used to estimate the variance
#' components when fastMode is TRUE. Ignored when Vu and Ve are supplied, and
#' when fastMode is FALSE. A value of zero uses every record, which is only
#' advisable for small training populations.
#' @param varComp an earlier fit of this model, from which Vu and Ve are
#' taken in place of supplying them, so that only the marker effects are
#' solved for. Not used together with Vu and Ve. The fit must be of this
#' model or of its counterpart ending in 2, to the same number of traits, on
#' the same markers. See Vu and the genomic selection vignette.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' traits
#'
#' @details
#' \subsection{The model}{
#' The model is
#' \deqn{y = X\beta + hb + M_fu_f + M_mu_m + M_du_d + e,}
#' where y holds the records, X is the design matrix of the fixed effects,
#' formed from the levels of the population's fixEff slot, \eqn{M_f} and
#' \eqn{M_m} hold the marker dosages an individual inherited from its female
#' and male parents, \eqn{M_d} holds the dominance codings, \eqn{u_f},
#' \eqn{u_m} and \eqn{u_d} are their effects, and e the errors. h is each
#' individual's sum of dominance codings, and its fixed regression b
#' captures directional dominance, as in \code{\link{RRBLUP_D}}. The
#' dominance effects allow cross-specific deviations from the parents'
#' GCA \insertCite{massman2013Genomewide}{AlphaSimR}. The marker effects are
#' the BLUPs \insertCite{henderson1975Best}{AlphaSimR} given variance
#' components estimated by REML \insertCite{patterson1971Recovery}{AlphaSimR}.
#' }
#'
#' \subsection{One trait}{
#' When traits gives a single trait, or is a function returning a single
#' value, \eqn{u_f \sim N(0, I\sigma^2_f)}, \eqn{u_m \sim N(0, I\sigma^2_m)},
#' \eqn{u_d \sim N(0, I\sigma^2_d)} and \eqn{e \sim N(0, I\sigma^2_e)}. The
#' four variances are estimated by average information REML
#' \insertCite{gilmour1995Average,johnson1995Restricted}{AlphaSimR}, with
#' EM-REML steps taken where an average information step fails to increase
#' the likelihood, for at most maxIter iterations. The model is fitted in
#' whichever of two equivalent forms is cheaper: on the variance matrix of the
#' records when there are fewer records than marker effects, and on
#' Henderson's mixed model equations otherwise.
#' }
#'
#' \subsection{Several traits}{
#' When traits gives more than one trait, the traits are fitted together as a
#' multiple trait model, which needs every trait recorded on every
#' individual. The female, male and dominance effects of each marker and the
#' errors of each record each have an unstructured covariance matrix among
#' the traits. These are estimated by REML, with EM steps followed by Newton
#' steps using the average information matrix
#' \insertCite{gilmour1995Average}{AlphaSimR}, for at most maxIter
#' iterations. The Vu slot of the result holds the covariance matrices of the
#' marker effects as a block diagonal matrix, with rows and columns named
#' "female:trait", "male:trait" and "d:trait".
#' }
#'
#' \subsection{Using the solution}{
#' The solution can be applied with \code{\link{setEBV}}:
#' \code{value="female"} or \code{value="male"} gives parental
#' predictions, while \code{value="gv"} predicts total hybrid genetic value.
#' These model predictions obviously differ from the observed cross means
#' returned by \code{\link{calcGCA}}. The dominance effects contribute to
#' the predicted genetic value, but SCA's are not returned.
#' }
#'
#' \subsection{Fast mode}{
#' Setting fastMode to TRUE gives a fast and low memory alternative for
#' training populations too large to fit directly. It is used only where it
#' is expected to be faster, see fastMode. The mixed model equations
#' of every random term and trait are then solved together by preconditioned
#' conjugate gradient, iterating over the genotypes rather than over a stored
#' coefficient matrix \insertCite{stranden1999Solving}{AlphaSimR}, so memory
#' use stays close to the size of the genotypes. The variance components are
#' estimated by REML on a random subset of the records, set by nSubsample,
#' with the same method as outside fast mode, unless Vu and Ve are supplied,
#' and the marker effects are always solved for using every record. Given the
#' variance components, the result is the BLUP of the same model that is
#' fitted when fastMode is FALSE.
#' }
#'
#' \subsection{Known variance components}{
#' Retraining a model whose variance components are already known,
#' for example from an earlier fit in the same simulation, needs only
#' the solve for the marker effects. Pass the earlier fit as varComp, or
#' give Vu and Ve. With fastMode TRUE the solve is by fast mode, for any
#' number of traits. With fastMode FALSE the mixed model equations of one
#' trait are solved directly, on whichever of the markers and the records
#' is cheaper, with Vu holding one variance per term. A variance of zero,
#' as REML often estimates for a dominance term, leaves its term out. This
#' takes the place of \code{\link{RRBLUP_SCA2}} with useEM FALSE, which is
#' deprecated, and Vu is its VuF, VuM and VuD in that order.
#' }
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
#' ans = RRBLUP_SCA(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#'
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#'
#' #Run GS model in fast mode
#' ans = RRBLUP_SCA(pop, fastMode=TRUE, simParam=SP)
#'
#' @export
RRBLUP_SCA = function(pop, traits=1, use="pheno", snpChip=1,
                      useQtl=FALSE, maxIter=NULL,
                      simParam=NULL, nThreads=NULL, ...,
                      fastMode=FALSE, Vu=NULL, Ve=NULL, nSubsample=3000L,
                      varComp=NULL){
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  # The model splits each individual's genome into one maternal and one
  # paternal haplotype, which is the design of a diploid single cross
  if(pop@ploidy!=2L){
    stop("RRBLUP_SCA: only diploids are supported")
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)
  # A function of the traits can return a plain vector, which has no columns
  # to count
  y = as.matrix(y)

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
  nTraits = ncol(y)
  traits = responseNames(traits, nTraits)
  if(nTraits>1){
    checkCompleteRecords(y)
  }
  # Variance components supplied outside fast mode are held fixed by the
  # direct solver below. In fast mode they are passed to the fast route,
  # which supplying them always takes.
  fixedVc = numeric(0)
  if(!fastMode){
    fixedVc = directVarComp(Vu, Ve, varComp, "RRBLUP_SCA", 3L, nTraits,
                            lociPerChr, lociLoc)
  }
  if(fastMode){
    if(!is.null(varComp)){
      vc = varCompFromFit(varComp, Vu, Ve, "RRBLUP_SCA", nTraits, lociPerChr,
                          lociLoc)
      Vu = vc$Vu
      Ve = vc$Ve
    }
    fastMode = useFastRoute("RRBLUP_SCA", Vu, Ve, nrow(y), nLoci, nTraits,
                            length(unique(pop@fixEff)), nSubsample)
  }
  if(is.null(maxIter)){
    maxIter = if(fastMode) 1000L else 100L
  }
  if(fastMode){
    vc = fastMKVarComp(Vu, Ve, 3L, nTraits)
    subset = if(vc$est) fastSubset(nSubsample, nrow(y)) else 1L
    ans = callFastRRBLUP_SCA(y, fixEff, pop, lociPerChr, lociLoc,
                             vc$Vu, vc$Ve, maxIter, vc$est,
                             as.integer(subset), nThreads)
  }else{
    ans = callRRBLUP_SCA(y, fixEff, pop,
                         lociPerChr, lociLoc, maxIter,
                         nThreads, fixedVc)
  }

  a1 = as.matrix(ans$a1)
  a2 = as.matrix(ans$a2)
  d = as.matrix(ans$d)
  alpha1 = as.matrix(ans$alpha1)
  alpha2 = as.matrix(ans$alpha2)
  gv = female = male = vector("list",nTraits)

  for(i in seq_len(nTraits)){
    gv[[i]] = new("TraitA2D",
                  nLoci=nLoci,
                  lociPerChr=lociPerChr,
                  lociLoc=lociLoc,
                  addEff=a1[,i],
                  addEffMale=a2[,i],
                  domEff=d[,i],
                  intercept=c(ans$mu)[i],
                  name=paste0("est_GV_",traits[i]))

    female[[i]] = new("TraitA",
                      nLoci=nLoci,
                      lociPerChr=lociPerChr,
                      lociLoc=lociLoc,
                      addEff=alpha1[,i],
                      intercept=c(ans$beta1)[i],
                      name=paste0("est_female_",traits[i]))

    male[[i]] = new("TraitA",
                    nLoci=nLoci,
                    lociPerChr=lociPerChr,
                    lociLoc=lociLoc,
                    addEff=alpha2[,i],
                    intercept=c(ans$beta2)[i],
                    name=paste0("est_male_",traits[i]))
  }

  if(nTraits>1){
    Vu = blockVu(ans$Vu, c("female","male","d"), traits)
  }else{
    Vu = as.matrix(ans$Vu)
  }

  output = new("RRsol",
               gv = gv,
               bv = as.list(NULL),
               female = female,
               male = male,
               Vu = Vu,
               Ve = traitMatrix(ans$Ve, traits),
               model = "RRBLUP_SCA")

  return(output)
}

#' @title RR-BLUP SCA Model 2
#'
#' @description
#' This function is deprecated. Use \code{\link{RRBLUP_SCA}}, which takes
#' known variance components as Vu and Ve or as an earlier fit in varComp.
#'
#' Fits the single trait model of \code{\link{RRBLUP_SCA}} with variance
#' components estimated by the EM algorithm, or taken as given. Only
#' diploids are supported. See details.
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
#' @details
#' \subsection{The model}{
#' The model is the single trait model of \code{\link{RRBLUP_SCA}},
#' \deqn{y = X\beta + hb + M_fu_f + M_mu_m + M_du_d + e,}
#' with \eqn{u_f \sim N(0, I\sigma^2_f)}, \eqn{u_m \sim N(0, I\sigma^2_m)},
#' \eqn{u_d \sim N(0, I\sigma^2_d)} and \eqn{e \sim N(0, I\sigma^2_e)}, where
#' \eqn{M_f} and \eqn{M_m} hold the marker dosages an individual inherited
#' from its female and male parents, \eqn{M_d} the dominance codings, and h
#' is each individual's sum of dominance codings. The marker effects are the
#' BLUPs \insertCite{henderson1975Best}{AlphaSimR} given the variance
#' components. When useEM is TRUE, the variance components are estimated by
#' EM-REML \insertCite{dempster1977Maximum,patterson1971Recovery}{AlphaSimR},
#' starting from VuF, VuM, VuD and Ve and iterating until the change is below
#' tol or maxIter is reached. When useEM is FALSE, VuF, VuM, VuD and Ve are
#' taken as the true values. Only one trait can be modeled.
#' }
#'
#' \subsection{Replacing it}{
#' This function was for retraining a model cheaply when the variance
#' components are already known. \code{\link{RRBLUP_SCA}} now does the same:
#' called with an earlier fit as varComp, or with Vu=c(VuF, VuM, VuD) and Ve,
#' and fastMode FALSE, it solves the mixed model equations directly with
#' the solver of its own REML fit and gives the same result, to rounding,
#' as this function with useEM FALSE. With fastMode TRUE it does so for
#' several traits as well. When the variance components are unknown,
#' use \code{\link{RRBLUP_SCA}} without them: the EM algorithm here is
#' generally far slower than its average information REML, and the default
#' maxIter is typically too small for convergence. When the algorithm fails
#' to converge a warning is displayed and the last iteration is returned.
#' }
#'
#' \subsection{Using the solution}{
#' The solution can be applied with \code{\link{setEBV}}:
#' \code{value="female"} or \code{value="male"} gives parental
#' predictions, while \code{value="gv"} predicts total hybrid genetic value.
#' These model predictions obviously differ from the observed cross means
#' returned by \code{\link{calcGCA}}. The dominance effects contribute to
#' the predicted genetic value, but SCA's are not returned.
#' }
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
#' #Run GS model, then refit with its variance components held fixed,
#' #which replaces RRBLUP_SCA2 with useEM=FALSE
#' ans = RRBLUP_SCA(pop, simParam=SP)
#' ans = RRBLUP_SCA(pop, varComp=ans, simParam=SP)
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
  .Deprecated(msg=numberedDeprecation("RRBLUP_SCA2", "RRBLUP_SCA",
                                      "c(VuF, VuM, VuD)"))
  # The model splits each individual's genome into one maternal and one
  # paternal haplotype, which is the design of a diploid single cross
  if(pop@ploidy!=2L){
    stop("RRBLUP_SCA2: only diploids are supported")
  }
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,nThreads=nThreads,...)

  # varA and varE are unnamed and indexed by position, so the numeric
  # index is taken before traits is replaced by its name.
  traitIndex = convertTraitsToIndex(traits, simParam)
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
  if(is.null(traitIndex)){
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
    stopifnot("this model fits a single trait, so traits must give one trait"=
                length(traits)==1)
    if(is.null(VuF)){
      VuF = 2*simParam$varA[traitIndex]/nLoci
      if(is.na(VuF)){
        VuF = var(y)/nLoci
      }
    }
    if(is.null(VuM)){
      VuM = 2*simParam$varA[traitIndex]/nLoci
      if(is.na(VuM)){
        VuM = var(y)/nLoci
      }
    }
    if(is.null(VuD)){
      VuD = simParam$varA[traitIndex]/nLoci
      if(is.na(VuD)){
        VuD = var(y)/nLoci/2
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traitIndex]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }

  #Fit model
  stopifnot("this model fits a single trait, so the response must have one column"=
              ncol(y)==1)
  ans = callRRBLUP_SCA2(y, fixEff, pop,
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
               Ve = as.matrix(ans$Ve),
               model = "RRBLUP_SCA")

  return(output)
}

#' @title Set estimated breeding values (EBV)
#'
#' @description
#' Adds estimated values to a population's EBV slot using
#' output from a genomic selection model or from
#' \code{\link{pedBLUP}}.
#' The estimated values can be either estimated
#' breeding values, estimated genetic values, or
#' estimated general combining values.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param solution an object of \code{\link{RRsol-class}} or
#' \code{\link{PedSol-class}}
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
#' @param onlyValues should only the estimated values be returned,
#' see return. When TRUE, append is ignored and the population is
#' not changed.
#'
#' @details
#' For output from \code{\link{pedBLUP}}, value "gv" gives the estimated
#' breeding value plus the estimate for the first level of the fixed
#' effects, and value "bv" gives the estimated breeding value centered on
#' the training population, or on targetPop when it is supplied. Values
#' "female" and "male" are not available. Individuals that were not in
#' the fitted pedigree get the average of their parents' estimates, which
#' is their best linear unbiased prediction.
#'
#' @return Returns an object of \code{\link{Pop-class}} if
#' onlyValues=FALSE, if onlyValues=TRUE a matrix is returned
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
#' #Return the estimates without changing the population
#' est = setEBV(pop, ans, onlyValues=TRUE, simParam=SP)
#'
#' @export
setEBV = function(pop, solution, value="gv", targetPop=NULL,
                  append=FALSE, simParam=NULL, nThreads=NULL,
                  onlyValues=FALSE){
  simParam = getSimParam(simParam)
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  if(is(solution, "PedSol")){
    return(.setEBVPed(pop=pop, solution=solution, value=value,
                      targetPop=targetPop, append=append,
                      simParam=simParam, onlyValues=onlyValues))
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

  if(onlyValues){
    return(ebv)
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
#' on what the whole simulation needs. Every model reads its markers
#' straight from the population's genotypes, so fitting one does not copy
#' the genotypes of the loci that are not markers.
#'
#' @param nInd the number of individuals in the training population
#' @param nMarker the number of markers per individual
#' @param model the model being fitted, given as the name of the function
#' that fits it. One of "RRBLUP", "RRBLUP2", "RRBLUP_D", "RRBLUP_D2",
#' "RRBLUP_GCA", "RRBLUP_GCA2", "RRBLUP_SCA" or "RRBLUP_SCA2".
#' The older values "REG", "GCA" and "SCA" are still accepted and are read as
#' "RRBLUP", "RRBLUP_GCA" and "RRBLUP_SCA", and "fastRRBLUP" is read as
#' "RRBLUP" with fastMode TRUE.
#' @param nTraits the number of traits fitted at once.
#' \code{\link{RRBLUP}}, \code{\link{RRBLUP_D}},
#' \code{\link{RRBLUP_GCA}} and \code{\link{RRBLUP_SCA}} fit more than one;
#' the numbered models fit one.
#' @param nFixEff the number of fixed effect levels, which is the number of
#' distinct values in the population's fixEff slot.
#' @param fastMode whether to estimate the memory of the fast mode of
#' \code{\link{RRBLUP}}, \code{\link{RRBLUP_D}}, \code{\link{RRBLUP_GCA}}
#' or \code{\link{RRBLUP_SCA}}. Ignored by every other model. A model given
#' fastMode TRUE uses fast mode only where it is expected to be faster, and
#' otherwise needs the memory estimated with fastMode FALSE.
#' @param nSubsample the number of records used to estimate variance
#' components in fast mode. Ignored when fastMode is FALSE, by the models
#' without a fast mode, and in fast mode when Vu and Ve are supplied.
#' @param useEM whether the variance components will be estimated. FALSE is
#' for a numbered model with useEM FALSE, or for an unnumbered model given Vu
#' and Ve, or varComp, with fastMode FALSE, which solves one trait without
#' estimating them. Ignored in fast mode, where nSubsample and supplying Vu
#' and Ve decide it.
#'
#' @details
#' The models differ in what they have to hold in memory, and the differences
#' are large enough to decide which one is usable.
#'
#' \code{\link{RRBLUP}} with fastMode TRUE reads the genotypes as one byte
#' per locus, then frees them and iterates over a copy packed to one bit per
#' allele, so it holds no square matrix at all. Its estimate is
#' dominated by the genotypes themselves and by the subset of records used
#' for the variance components. Fitted to several traits, it holds a few
#' vectors of markers and of records per trait, and the variance components
#' are estimated as for several traits in \code{\link{RRBLUP}}, on the
#' subset of records.
#'
#' \code{\link{RRBLUP}} decomposes a square matrix whose dimensions are the
#' smaller of the number of records and the number of markers, because a
#' matrix and its transpose share their nonzero eigenvalues. Fitted to several
#' traits, it does the same with a thin singular value decomposition of the
#' markers when they are fewer than the records.
#'
#' The numbered models estimate variance components with Henderson's mixed
#' model equations, whose coefficient matrix is square with dimensions equal
#' to the number of fixed effects plus the number of random effects, so they
#' grow with the number of markers rather than the number of records.
#' Estimating variance components needs the trace of part of the inverse of
#' that matrix, which is reached through the inverse of its Cholesky factor,
#' so it costs one more matrix of the same size.
#'
#' With the variance components fixed, which is useEM false, the numbered
#' models write the model on the records instead when there are fewer
#' records than marker effects, so their estimate depends on useEM, and with
#' useEM false it grows with the number of records once the marker effects
#' outnumber them. The unnumbered models given Vu and Ve with fastMode FALSE
#' use the solver of the REML fit of the GCA, SCA and dominance models,
#' with a single term for \code{\link{RRBLUP}}, and are estimated as that
#' fit is, which is an upper bound.
#'
#' The GCA, SCA and dominance models with fastMode TRUE hold the markers
#' packed to one bit per allele, assumed here to be a diploid's two bits per
#' marker, and no square matrix of the markers or the records. Their
#' iteration holds a few matrices of markers by terms and traits, and a block
#' of the preconditioner for every marker, whose size is the square of the
#' number of terms times the number of traits. Before that, the variance
#' components are estimated by the same model on the subset of records,
#' given each term's markers or, with more markers than records in the
#' subset, a factor of their cross product, which holds a square matrix of
#' the subset for each term.
#'
#' The GCA, SCA and dominance models fit more than one random effect. They
#' work with whichever of two forms is cheaper: a square matrix of the
#' records, holding several of them when there are more marker effects than
#' records, or Henderson's equations, whose coefficient matrix grows with the
#' number of marker effects, when there are fewer. They are the most
#' demanding of the set when the marker effects outnumber the records.
#' Fitted to several traits, both forms grow with the number of traits: the
#' square matrix of the records has one row per record and trait, and
#' Henderson's equations one per marker effect and trait, so the memory
#' needed grows with the square of the number of traits.
#'
#' @return Returns an estimate for the required gigabytes of RAM
#'
#' @examples
#' RRBLUPMemUse(nInd=1000, nMarker=5000)
#'
#' # The same training data, fitted several ways
#' RRBLUPMemUse(nInd=5000, nMarker=2000, model="RRBLUP", fastMode=TRUE)
#' RRBLUPMemUse(nInd=5000, nMarker=2000, model="RRBLUP")
#' RRBLUPMemUse(nInd=5000, nMarker=2000, model="RRBLUP_SCA")
#' RRBLUPMemUse(nInd=5000, nMarker=2000, model="RRBLUP_SCA", fastMode=TRUE)
#'
#' # RRBLUP with more markers than records, estimating the variance
#' # components and given them
#' RRBLUPMemUse(nInd=1000, nMarker=10000, model="RRBLUP")
#' RRBLUPMemUse(nInd=1000, nMarker=10000, model="RRBLUP", useEM=FALSE)
#'
#' @export
RRBLUPMemUse = function(nInd, nMarker, model="RRBLUP", nTraits=1L,
                        nFixEff=1L, fastMode=FALSE, nSubsample=3000L,
                        useEM=TRUE){
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
  # fastRRBLUP was a function of its own before it became RRBLUP's fast mode
  if(model=="fastRRBLUP"){
    model = "RRBLUP"
    fastMode = TRUE
  }
  # Outside fast mode, a model given its variance components solves one trait
  # with the solver of its REML fit, which for RRBLUP is the multikernel
  # solver with one term rather than the spectral one, evaluated once. The
  # estimate for the REML fit holds a little more than a single evaluation
  # does, so it serves as an upper bound.
  fixedVc = !fastMode & !useEM &
    (model %in% c("RRBLUP","RRBLUP_D","RRBLUP_GCA","RRBLUP_SCA"))
  if(fixedVc & (nTraits>1)){
    stop("variance components can be supplied for several traits only ",
         "when fastMode is TRUE")
  }

  # Genotypes are read as one byte per locus before being turned into
  # dosages, so both are held while the conversion runs
  genoBytes = n*m

  if((model=="RRBLUP") & fastMode){
    nSub = as.double(nSubsample)
    if((nSub<=0) | (nSub>n)){
      nSub = n
    }
    if(nT>1){
      # Every trait has its own vectors of markers in the iteration, seven of
      # them, its share of the tables the products are taken with, and a few
      # vectors of records
      doubles = (11*m+4*n)*nT+q*n
      # Variance components by the multiple trait REML of RRBLUP on the
      # subset of records, with its optimizer and fixed effects
      optimDoubles = 8*(nT*(nT+1))^2
      fixDoubles = 2*nT*q^2+4*q^2+5*nSub*q
      if(m<nSub){
        doubles = doubles+4*nSub*m+nT^2*nSub+optimDoubles+fixDoubles
      }else{
        doubles = doubles+4*nSub^2+nT^2*nSub+optimDoubles+fixDoubles
      }
    }else{
      # No dosages are stored and no square matrix is formed. The solver
      # holds six vectors of markers plus the column means, and a few of
      # records.
      doubles = 7*m+3*n+q*n
      # Variance components come from a subset of the records. Whichever
      # cross product is smaller is the one decomposed.
      if(m<nSub){
        doubles = doubles+nSub*m+5*m^2
      }else{
        doubles = doubles+3*nSub^2
      }
    }
  }else if((model=="RRBLUP") & !fixedVc){
    if(nT>1){
      # The Newton steps of the REML fit work with square matrices over the
      # t*(t+1) parameters of Vu and Ve, which grow with the fourth power of
      # the number of traits: the average information, the Jacobian and
      # curvature of the Cholesky parameterization, the Hessian and a
      # product formed on the way to it, the Hessian's eigenvectors and the
      # eigensolver's workspace
      optimDoubles = 8*(nT*(nT+1))^2
      # The fixed effects: each likelihood evaluation holds a q by q inverse
      # per trait, and the current and trial points coexist, with a few
      # q by q temporaries and copies of the fixed design on the records,
      # its orthonormal basis, its rotation and its part in the null space
      fixDoubles = 2*nT*q^2+4*q^2+5*n*q
      if(m<n){
        # Rotation by a thin singular value decomposition of the dosages,
        # which holds the left singular vectors twice while the ones that
        # are kept are copied out, plus one small inverse per record
        doubles = 4*n*m+nT^2*n+optimDoubles+fixDoubles
      }else{
        # Rotation by the eigenvectors of the record cross product, which
        # is held, symmetrized and decomposed
        doubles = n*m+4*n^2+nT^2*n+optimDoubles+fixDoubles
      }
    }else if(m<n){
      # Decomposition of the marker cross product, which is copied into the
      # eigensolver, with its eigenvectors and two matrices derived from them
      doubles = n*m+5*m^2
    }else{
      # Decomposition of the record cross product, which is copied into the
      # eigensolver, with its eigenvectors
      doubles = n*m+3*n^2
    }
  }else if(model=="RRBLUP2"){
    if(useEM | (m<=n)){
      # Henderson's equations, whose coefficient matrix is the fixed effects
      # plus the markers. The EM algorithm also holds the inverse of its
      # Cholesky factor, because that is what supplies the trace it needs
      if(useEM){
        doubles = n*m+3*(q+m)^2
      }else{
        doubles = n*m+2*(q+m)^2
      }
    }else{
      # The model written on the records, which holds a square matrix of the
      # records instead. Reached only when the variance components are fixed
      # and there are fewer records than markers
      doubles = n*m+2*n^2
    }
  }else if(model %in% c("RRBLUP_D","RRBLUP_GCA","RRBLUP_SCA") & fastMode){
    nKernel = if(model=="RRBLUP_SCA") 3 else 2
    nSub = as.double(nSubsample)
    if((nSub<=0) | (nSub>n)){
      nSub = n
    }
    # No dosages are decoded. The markers are packed to one bit per allele,
    # a quarter of a byte per marker for a diploid
    genoBytes = n*m/4
    # The iteration: the stacked effects, residuals, preconditioned
    # residuals, search directions and their products, with a row per term
    # and trait; each term's effects and products on the scale of the traits
    # and its lookup tables; the block of the preconditioner and the cross
    # products of the terms for every marker; a few vectors of records per
    # trait; and the fixed effects with their orthonormal basis
    nR = nKernel*nT
    iterDoubles = 5*nR*m+6*nKernel*nT*m+nR^2*m+nKernel^2*m+3*nT*n+
      2*(q+1)*n
    # The variance components, by the direct form of the model on the subset,
    # given each term's codes or their cross product, which is no larger than
    # the codes of as many markers as records. Building a cross product holds
    # a copy of it while it is made symmetric.
    mSub = min(m, nSub)
    vcGB = RRBLUPMemUse(nInd=nSub, nMarker=mSub, model=model, nTraits=nTraits,
                        nFixEff=nFixEff)
    vcDoubles = (vcGB*10^9-nSub*mSub)/8
    if(m>=nSub){
      vcDoubles = vcDoubles+nSub^2
    }
    doubles = max(iterDoubles, vcDoubles)
  }else if(model %in% c("RRBLUP_D","RRBLUP_GCA","RRBLUP_SCA") & (nT>1)){
    nKernel = if(model=="RRBLUP_SCA") 3 else 2
    # The multiple trait form chooses between the same two sides, on the
    # costs of useCoefficientSide with a record and a coefficient per trait
    nRec = n*nT
    nCoef = nT*(q+nKernel*m)
    nPar = (nKernel+1)*nT*(nT+1)/2
    nEval = 10
    recordCost = nEval*(nRec^3+nRec^2*nCoef)
    coefCost = nRec*nCoef^2+nEval*(nCoef^3+nRec*nCoef*(nPar+2))
    # The optimizer's square matrices over the parameters, and the working
    # variates of the derivatives with their products with P, one vector of
    # records per parameter each
    optimDoubles = 8*nPar^2+2*nRec*nPar
    if(coefCost<recordCost){
      # The dosages and a scaled copy of them, the design of all the effects,
      # its cross product and that of each term, and the coefficient matrix
      # with the Cholesky factors of the current and trial points, its
      # inverse and the inverse of a factor
      c = q+nKernel*m
      doubles = 3*nKernel*n*m+(nKernel+1)*c^2+5*nCoef^2+optimDoubles
    }else{
      # The dosages and a scaled copy of them, the kernel of each term, and
      # the variance matrix of the records with the factors of the current
      # and trial points, the inverse of a factor and P
      doubles = 2*nKernel*n*m+nKernel*n^2+5*nRec^2+optimDoubles
    }
  }else if((model %in% c("RRBLUP_D","RRBLUP_GCA","RRBLUP_SCA")) | fixedVc){
    nKernel = switch(model, RRBLUP=1, RRBLUP_SCA=3, 2)
    # The solver works with whichever of the variance matrix of the records
    # and Henderson's equations is cheaper, and chooses by the same costs as
    # useCoefficientSide in src/mme.h
    nCoef = q+nKernel*m
    nEval = 10
    recordCost = nEval*(n^3+n^2*nCoef)
    coefCost = n*nCoef^2+nEval*(nCoef^3+n*nCoef*(nKernel+2))
    if(coefCost<recordCost){
      # Henderson's equations: a copy of the dosages as the design of all
      # the effects, its cross product, and the coefficient matrix with the
      # Cholesky factors of the current and trial points and the inverse of
      # a factor
      doubles = 2*nKernel*n*m+5*nCoef^2
    }else if(m<n){
      # The average information matrix is reached through the dosages
      doubles = 2*nKernel*n*m+(nKernel+3)*n^2
    }else{
      # and through square matrices of the records
      doubles = nKernel*n*m+(3*nKernel+3)*n^2
    }
  }else if(model %in% c("RRBLUP_D2","RRBLUP_GCA2","RRBLUP_SCA2")){
    nKernel = if(model=="RRBLUP_SCA2") 3 else 2
    if(useEM){
      doubles = nKernel*n*m+3*(q+nKernel*m)^2
    }else if(n<(nKernel*m)){
      # With the variance components fixed and fewer records than marker
      # effects, the model is written on the records, which holds the
      # variance matrix of the records, its inverse and a product of a term
      # with its transpose on the way to it
      doubles = nKernel*n*m+3*n^2
    }else{
      doubles = nKernel*n*m+2*(q+nKernel*m)^2
    }
  }else{
    stop(paste0("model=",model," not recognized"))
  }

  # Using base 10 rather than base 2, which is the more conservative of the
  # two and leaves a little room for what is not counted here
  return((8*doubles+genoBytes)/10^9) #GB
}
