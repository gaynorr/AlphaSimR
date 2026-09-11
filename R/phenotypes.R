#' Add residual error to genetic values
#'
#' @param gv matrix of genetic values
#' @param varE residual variances, vector or matrix
#' @param reps number of reps for phenotype
#'
#' @keywords internal
addError = function(gv, varE, reps){
  nTraits = ncol(gv)
  nInd = nrow(gv)
  if(is.matrix(varE)){
    stopifnot(isSymmetric(varE),
              ncol(varE)==nTraits)
    error = matrix(rnorm(nInd*nTraits),
                   ncol=nTraits)%*%transMat(varE)
  }else{
    stopifnot(length(varE)==nTraits)
    error = lapply(varE,function(x){
      if(is.na(x)){
        return(rep(NA_real_,nInd))
      }else{
        return(rnorm(nInd,sd=sqrt(x)))
      }
    })
    error = do.call("cbind",error)
  }
  error = error/sqrt(reps)
  pheno = gv + error

  return(pheno)
}


#' Calculate phenotypes
#'
#' @param pop an object of class Pop
#' @param varE a vector or matrix of residual variances
#' @param reps number of reps for phenotype
#' @param p p-value for environment
#' @param traits number of traits
#' @param simParam simulation parameters. If \code{NULL}, the function uses
#' the object named \code{SP} from the global environment.
#'
#' @keywords internal
calcPheno = function(pop, varE, reps, p, traits, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  nTraits = length(traits)

  if(nTraits==0L){
    return(pop@pheno)
  }

  gv = pop@gv
  for(i in seq_len(nTraits)){
    if(.hasSlot(simParam$traits[[traits[i]]], "envVar")){
      stdDev = sqrt(simParam$traits[[traits[i]]]@envVar)
      gv[,traits[i]] = gv[,traits[i]] +
        pop@gxe[[traits[i]]]*qnorm(p[i], sd=stdDev)
    }
  }
  gv = gv[,traits,drop=FALSE]

  # Calculate new phenotypes
  newPheno = addError(gv=gv, varE=varE, reps=reps)

  # Add to old phenotype
  pheno = pop@pheno
  pheno[,traits] = newPheno

  return(pheno)
}

#' @title Set phenotypes
#'
#' @description
#' Sets phenotypes for all traits by adding random error
#' from a multivariate normal distribution.
#'
#' @param pop an object of \code{\link{Pop-class}} or
#' \code{\link{HybridPop-class}}
#' @param h2 a vector of desired narrow-sense heritabilities for
#' each trait. See details.
#' @param H2 a vector of desired broad-sense heritabilities for
#' each trait. See details.
#' @param varE error (co)variances for traits. See details.
#' @param corE an optional matrix for correlations between errors.
#' See details.
#' @param reps number of replications for phenotype. See details.
#' @param fixEff fixed effect to assign to the population. Used
#' by genomic selection models only.
#' @param p the p-value for the environmental covariate
#' used by GxE traits. If NULL, a value is
#' sampled at random.
#' @param onlyPheno should only the phenotype be returned, see return
#' @param traits an integer vector indicate which traits to set. If NULL,
#' all traits will be set.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param ... additional arguments passed to the \code{finalizePheno}
#' function in simParam
#'
#' @details
#' There are three arguments for setting the error variance of a
#' phenotype: h2, H2, and varE. The user should only use one of these
#' arguments. If the user supplies values for more than one, only one
#' will be used according to order in which they are listed above.
#'
#' The h2 argument allows the user to specify the error variance
#' according to narrow-sense heritability. This calculation uses the
#' additive genetic variance and total genetic variance in the founder
#' population. Thus, the heritability relates to the founder population
#' and not the current population.
#'
#' The H2 argument allows the user to specify the error variance
#' according to broad-sense heritability. This calculation uses the
#' total genetic variance in the founder population. Thus, the heritability
#' relates to the founder population and not the current population.
#'
#' The varE argument allows the user to specify the error variance
#' directly. The user may supply a vector describing the error variance
#' for each trait or supply a matrix that specify the covariance of
#' the errors.
#'
#' The corE argument allows the user to specify correlations for the
#' error covariance matrix. These correlations are be supplied in addition
#' to the h2, H2, or varE arguments. These correlations will be used to
#' construct a covariance matrix from a vector of variances. If the user
#' supplied a covariance matrix to varE, these correlations will supercede
#' values provided in that matrix.
#'
#' The reps parameter is for convenient representation of replicated data.
#' It is intended to represent replicated yield trials in plant
#' breeding programs. In this case, varE is set to the plot error and
#' reps is set to the number of plots per entry. The resulting phenotype
#' represents the entry-means.
#'
#' @return Returns an object of \code{\link{Pop-class}} or
#' \code{\link{HybridPop-class}} if onlyPheno=FALSE, if
#' onlyPheno=TRUE a matrix is returned
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Add phenotype with error variance of 1
#' pop = setPheno(pop, varE=1)
#'
#' @export
setPheno = function(pop, h2=NULL, H2=NULL, varE=NULL, corE=NULL,
                    reps=1, fixEff=1L, p=NULL, onlyPheno=FALSE,
                    traits=NULL, simParam=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }

  # Determine which traits are selected
  if(is.null(traits)){
    if(simParam$nTraits>0L){
      traits = 1:simParam$nTraits
    }else{
      traits = integer()
    }
  }else{
    traits = as.integer(traits)
    stopifnot(all(traits>0L),
              all(!duplicated(traits)),
              max(traits)<=simParam$nTraits)
  }
  nTraits = length(traits)

  # Check for valid length of reps vector
  if(length(reps)==1){
    reps = rep(reps, nTraits)
  }else{
    stopifnot(length(reps)==nTraits)
  }

  # Set p-value for GxE traits
  if(is.null(p)){
    p = rep(runif(1), nTraits)
  }else if(length(p)==1){
    p = rep(p, nTraits)
  }else{
    stopifnot(length(p)==nTraits)
  }

  # Calculate varE if using h2 or H2
  if(!is.null(h2)){
    if(length(h2)==1){
      h2 = rep(h2, nTraits)
    }
    varA = simParam$varA[traits]
    varG = simParam$varG[traits]

    stopifnot(length(h2)==nTraits,
              all(varA>0),
              all(varG>0))
    varE = numeric(nTraits)
    for(i in seq_len(nTraits)){
      tmp = varA[i]/h2[i]-varG[i]
      if(tmp<0){
        stop(paste0("h2=",h2[i]," is not possible for trait ",traits[i]))
      }
      varE[i] = tmp
    }
  }else if(!is.null(H2)){
    if(length(H2)==1){
      H2 = rep(H2, nTraits)
    }
    varG = simParam$varG[traits]

    stopifnot(length(H2)==nTraits)
    varE = numeric(nTraits)
    for(i in seq_len(nTraits)){
      tmp = varG[i]/H2[i]-varG[i]
      varE[i] = tmp
    }
  }else if(!is.null(varE)){
    if(is.matrix(varE)){
      stopifnot(nTraits==nrow(varE),
                isSymmetric(varE))
    }else{
      stopifnot(length(varE)==nTraits)
    }
  }else{
    if(is.matrix(simParam$varE)){
      varE = simParam$varE[traits, traits]
    }else{
      varE = simParam$varE[traits]
    }
  }

  # Set error correlations
  if(!is.null(corE)){
    if(is.matrix(varE)){
      varE = diag(varE)
    }
    stopifnot(length(varE)==nrow(corE),
              isSymmetric(corE))

    varE = diag(sqrt(varE),
                nrow=nTraits,
                ncol=nTraits)
    varE = varE%*%corE%*%varE
  }


  # Use lapply if object is a MultiPop
  # Only passing varE after previous processing
  if(is(pop,"MultiPop")){
    stopifnot(!onlyPheno)
    pop@pops = lapply(pop@pops, setPheno, h2=NULL, H2=NULL,
                      varE=varE, corE=NULL, reps=reps, fixEff=fixEff,
                      p=p, traits=traits, simParam=simParam)
    return(pop)
  }

  # Create phenotypes
  pheno = calcPheno(pop=pop, varE=varE, reps=reps, p=p,
                    traits=traits, simParam=simParam)

  pheno = simParam$finalizePheno(pheno, pop=pop, simParam=simParam, ...)

  colnames(pheno) = colnames(pop@gv)

  if(onlyPheno){
    return(pheno)
  }

  pop@pheno = pheno

  if(is(pop,"Pop")){
    pop@fixEff = rep(as.integer(fixEff), pop@nInd)
  }

  return(pop)
}

#' @title Convert a normal (Gaussian) trait to a log-normal trait
#' @param x matrix, values for one or more traits (if not a matrix,
#'   we cast to a matrix).
#' @param meanLogShift \code{NULL}, numeric or list, additional additive shift(s)
#'  on the latent log scale; when \code{NULL} a shift of 0 is assumed, when
#'  numeric shifts for all traits in \code{x} must be provided, and when list
#'  shifts for all traits in \code{x} must be provided with possibility to pass
#'  a \code{NULL} list node to skip the conversion for the trait (see examples).
#' @details If input trait is normal (Gaussian) then this function generates
#'   a log-normal trait by applying exponential link function on the input.
#'   No sampling happens in this function, which makes it deterministic.
#'
#'   Note the possible terminological confusion, a log-normal trait is expressed
#'   on exponential scale and its underlying (latent) values are on the log scale
#'   (see examples). If the supplied latent values \code{x} have mean
#'   \code{mu} and variance \code{sigma2}, then the recoded trait has expected
#'   value \code{exp(meanLogShift + mu + sigma2/2)} and variance
#'   \code{E(y)^2 * (exp(sigma2) - 1)}; therefore, to target an observed mean
#'   \code{M}, set \code{meanLogShift = log(M) - mu - sigma2/2}. See
#'   \code{\link{rlnorm}} for the same expressions in the standard log-normal
#'   parameterization. Note that \code{asLogNormal} does not provide the
#'   \code{sdlog} argument because latent trait variation is already controlled
#'   by other parameters.
#'
#'   The name \code{meanLogShift} is used to emphasize that this argument
#'   is an additional mean shift applied during transformation, not the
#'   primary way to set the latent trait mean. In normal AlphaSimR workflow,
#'   the latent mean is usually already set via
#'   \code{SP$addTrait*(..., mean = ...)} in founding population and
#'   \code{meanLogShift} should be left at its default unless an extra
#'   transformation-specific shift on the latent (log) scale is needed.
#'   One example is to control the mean of the observed values as shown below.
#'   This value should be established at the start of simulation and kept
#'   constant for most use cases.
#' @return matrix of log-normal values.
#' @seealso \code{finalizePop} and \code{finalizePheno} functions in
#'   \code{\link{SimParam}} for automatic conversion (also demonstrated below).
#' @examples
#' #Simulate a founder pop, set latent trait parameters, and create a population
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' trtMeanLog = c(0, 0)
#' trtVarGLog = c(1, 2)
#' SP$addTraitA(nQtlPerChr = 10, mean = trtMeanLog, var = trtVarGLog,
#'              corA = matrix(data = c(1.0, 0.6,
#'                                     0.6, 1.0), ncol = 2))
#' trtVarELog = c(1, 1)
#' trtVarPLog = trtVarGLog + trtVarELog
#' SP$setVarE(varE = trtVarELog)
#' pop = newPop(founderPop)
#' popLarge = randCross(pop, nCrosses = 1000)
#' 
#' meanVarFun = function(x) list(mean = mean(x), var = var(x))
#' 
#' #Latent phenotypes and parameters
#' (phenoLog = pheno(pop))
#' phenoLogLarge = pheno(popLarge)
#' apply(X = phenoLog, MARGIN = 2, FUN = meanVarFun)
#' apply(X = phenoLogLarge, MARGIN = 2, FUN = meanVarFun)
#'
#' #Convert a single input trait
#' (phenoExpMeanLog0 = asLogNormal(pheno(pop)[, 1]))
#' meanVarFun(phenoExpMeanLog0)
#'
#' #Demonstrate meanLogShift argument
#' #Here we aim to obtain observed trait values with the mean of 1.
#' #Since E(trtExp)=exp(meanLogShift + mean(trtLog) + var(trtLog)/2),
#' #trait 1 is centered on the latent scale (mean(trtLog)=0), and
#' #trtVarPLog[1] is the latent variance, we set meanLogShift to -trtVarPLog[1]/2.
#' #See also ?rlnorm for the expression of variance.
#' (phenoExpMeanExp1 = asLogNormal(pheno(pop)[, 1], meanLogShift = -trtVarPLog[1]/2))
#' meanVarFun(phenoExpMeanExp1)
#' exp(trtVarPLog[1] - 1)
#' cbind(phenoLog = phenoLog[, 1],
#'   phenoExpMeanLog0 = phenoExpMeanLog0,
#'   phenoExpMeanExp1 = phenoExpMeanExp1)
#' 
#' tmp = cbind(phenoLog = phenoLogLarge[, 1],
#'   phenoExpMeanLog0 = c(asLogNormal(phenoLogLarge[, 1])),
#'   phenoExpMeanExp1 = c(asLogNormal(phenoLogLarge[, 1], meanLogShift = -trtVarPLog[1]/2)))
#' (tmp2 = apply(X = tmp, MARGIN = 2, FUN = meanVarFun))
#' hist(tmp[, "phenoLog"], main = paste0("Mean: ", v=tmp2$phenoLog$mean))
#' abline(v=tmp2$phenoLog$mean, col = "red")
#' hist(tmp[, "phenoExpMeanLog0"], main = paste0("Mean: ", v=tmp2$phenoExpMeanLog0$mean))
#' abline(v = tmp2$phenoExpMeanLog0$mean, col = "red")
#' hist(tmp[, "phenoExpMeanExp1"], main = paste0("Mean: ", v=tmp2$phenoExpMeanExp1$mean))
#' abline(v= tmp2$phenoExpMeanExp1$mean, col = "red")
#'
#' #Convert multiple input traits
#' asLogNormal(pheno(pop))
#' try(asLogNormal(pheno(pop), meanLogShift = 0))
#' asLogNormal(pheno(pop), meanLogShift = c(0, 1))
#' asLogNormal(pheno(pop), meanLogShift = list(0, NULL))
#' 
#' #Store the recoded trait manually
#' pheno(pop)
#' pop@pheno[, 1] = asLogNormal(pheno(pop)[, 1])
#' pheno(pop)
#' 
#' #Apply and store the transformation automatically via SimParam$finalizePop()
#' finalizePopDefault = SP$finalizePop
#' SP$finalizePop = function(pop, simParam = SP, ...) {
#'   pop@pheno[, 1] = asLogNormal(pheno(pop)[, 1])
#'   return(pop)
#' }
#' pop = newPop(founderPop)
#' pheno(pop)
#' 
#' #Apply and store the transformation automatically via SimParam$finalizePheno()
#' SP$finalizePop = finalizePopDefault
#' SP$finalizePheno = function(pheno, pop, simParam = SP, ...) {
#'   pheno[, 1] = asLogNormal(pheno[, 1])
#'   return(pheno)
#' }
#' pop = newPop(founderPop)
#' pheno(pop)
#' @export
asLogNormal <- function(x, meanLogShift = NULL) {
  if (!is.matrix(x)) {
    x = as.matrix(x)
  }
  nTraits = ncol(x)
  if (is.null(meanLogShift)) {
    meanLogShift = rep(x = 0, times = nTraits)
  }
  if (is.numeric(meanLogShift)) {
    if (length(meanLogShift) != nTraits) {
      stop("You must supply meanLogShift for all traits in x!")
    }
    for (trt in 1:nTraits) {
      x[, trt] = exp(meanLogShift[trt] + x[, trt])
    }
  } else if (is.list(meanLogShift)) {
    if (length(meanLogShift) != nTraits) {
      stop("You must supply meanLogShift for all traits in x!")
    }
    for (trt in 1:nTraits) {
      if (!is.null(meanLogShift[[trt]])) {
        x[, trt] = exp(meanLogShift[[trt]] + x[, trt])
      }
    }
  } else {
    stop("meanLogShift must be NULL, numeric, or list!")
  }
  return(x)
}

#' @title Convert a normal (Gaussian) trait to an ordered categorical (threshold)
#'   trait
#' @param x matrix, values for one or more traits (if not a matrix,
#'   we cast to a matrix).
#' @param p \code{NULL}, numeric, or list, when \code{NULL} \code{threshold}
#'   is used; when numeric, provide a vector of category probabilities to
#'   convert continuous values into for a single trait (if probabilities
#'   do not sum to 1, another category is added and a warning is raised);
#'   when list, provide a list of probabilities - list node with \code{NULL}
#'   will skip conversion for a specific trait (see examples).
#'   If \code{p} is provided, it takes precedence over \code{threshold}.
#'   Internally \code{p} is converted to \code{threshold}, and any supplied
#'   \code{threshold} values are ignored.
#' @param mean numeric, assumed latent mean(s) of \code{x}; used only when
#'   \code{p} is given to convert category probabilities to thresholds.
#'   See also details.
#' @param var numeric, assumed latent variance(s) of \code{x}; used only when
#'   \code{p} is given to convert category probabilities to thresholds.
#'   See also details.
#' @param threshold \code{NULL}, numeric or, list, when numeric, provide
#'   a vector of category thresholds to convert continuous values into for
#'   a single trait (the thresholds specify left-closed and right-opened
#'   intervals [t1, t2), which can be changed with \code{include.lowest}
#'   and \code{right}; ensure you add \code{-Inf} and \code{Inf} or min and
#'   max to cover the whole range of values; otherwise you will get
#'   \code{NA} values);
#'   when list, provide a list of numeric thresholds - list node with \code{NULL}
#'   will skip conversion for a specific trait (see examples). The default values
#'   are set somewhat arbitrarily to get a ratio of (0.16, 0.68, 0.16) of records
#'   in each of the categories with most of the outlying individuals scored
#'   differently than individuals close to the "average".
#' @param include.lowest logical, see \code{\link{cut}}.
#' @param right logical, see \code{\link{cut}}.
#' @details If input trait is normal (Gaussian) then this function generates a
#'   categorical trait according to the ordered probit model.
#'   No sampling happens in this function, which makes it deterministic.
#'
#'   When \code{p} is used, \code{mean} and \code{var} describe the latent
#'   distribution of \code{x} and are used only to derive thresholds.
#'   In normal AlphaSimR workflow, this latent mean and variance are usually
#'   set via \code{SP$addTrait*(..., mean = ..., var = ...)} in founding
#'   population and \code{SP$setVarE}. \code{p} or \code{threshold} values
#'   should be established at the start of simulation and kept constant
#'   for most use cases.
#' @return matrix of values with some traits recorded as ordered categories
#'  in the form of \code{1:nC} with \code{nC} being the number of categories.
#' @seealso \code{finalizePop} and \code{finalizePheno} functions in
#'   \code{\link{SimParam}} for automatic conversion (also demonstrated below).
#' @examples
#' #Simulate a founder pop, set latent trait parameters, and create a population
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' trtMean = c(0, 0)
#' trtVarG = c(1, 2)
#' SP$addTraitA(nQtlPerChr = 10, mean = trtMean, var = trtVarG,
#'              corA = matrix(data = c(1.0, 0.6,
#'                                     0.6, 1.0), ncol = 2))
#' trtVarE = c(1, 1)
#' trtVarP = trtVarG + trtVarE
#' SP$setVarE(varE = trtVarE)
#' pop = newPop(founderPop)
#' pheno(pop)
#'
#' #Convert a single input trait
#' asCategorical(pheno(pop)[, 1])
#'
#' #Demonstrate threshold argument (in units of pheno SD)
#' asCategorical(pheno(pop)[, 1], threshold = c(-1, 0, 1) * sqrt(trtVarP[1]))
#' asCategorical(pheno(pop)[, 1], threshold = c(-Inf, -1, 0, 1, Inf) * sqrt(trtVarP[1]))
#' asCategorical(pheno(pop)[, 1], threshold = c(-Inf, 0, Inf))
#'
#' #Demonstrate p argument
#' asCategorical(pheno(pop)[, 1], p = 0.5, var = trtVarP[1])
#' asCategorical(pheno(pop)[, 1], p = c(0.5, 0.5), var = trtVarP[1])
#' asCategorical(pheno(pop)[, 1], p = c(0.25, 0.5, 0.25), var = trtVarP[1])
#'
#' #Convert multiple input traits (via threshold or p argument)
#' try(asCategorical(pheno(pop)))
#' asCategorical(pheno(pop),
#'               threshold = list(c(-Inf, 0, Inf),
#'                                NULL))
#' try(asCategorical(pheno(pop), p = c(0.5, 0.5)))
#' asCategorical(pheno(pop),
#'               p = list(c(0.5, 0.5),
#'                        NULL),
#'               mean = trtMean, var = trtVarP)
#'
#' asCategorical(pheno(pop),
#'               threshold = list(c(-Inf, 0, Inf),
#'                                c(-Inf, -2, -1, 0, 1, 2, Inf) * sqrt(trtVarP[2])))
#' q = c(-2, -1, 0, 1, 2)
#' p = pnorm(q)
#' p = c(p[1], p[2]-p[1], p[3]-p[2], p[4]-p[3], p[5]-p[4], 1-p[5])
#' asCategorical(pheno(pop),
#'               p = list(c(0.5, 0.5),
#'                        p),
#'               mean = trtMean, var = trtVarP)
#' 
#' #Store the recoded trait manually
#' pheno(pop)
#' pop@pheno[, 1] = asCategorical(pheno(pop)[, 1])
#' pheno(pop)
#' 
#' #Apply and store the transformation automatically via SimParam$finalizePop()
#' finalizePopDefault = SP$finalizePop
#' SP$finalizePop = function(pop, simParam = SP, ...) {
#'   pop@pheno[, 1] = asCategorical(pheno(pop)[, 1])
#'   return(pop)
#' }
#' pop = newPop(founderPop)
#' pheno(pop)
#' 
#' #Apply and store the transformation automatically via SimParam$finalizePheno()
#' SP$finalizePop = finalizePopDefault
#' SP$finalizePheno = function(pheno, pop, simParam = SP, ...) {
#'   pheno[, 1] = asCategorical(pheno[, 1])
#'   return(pheno)
#' }
#' pop = newPop(founderPop)
#' pheno(pop)
#' @export
asCategorical = function(x, p = NULL, mean = 0, var = 1,
                         threshold = c(-Inf, -sqrt(var), sqrt(var), Inf),
                         include.lowest = TRUE, right = FALSE) {
  if (!is.matrix(x)) {
    x = as.matrix(x)
  }
  nTraits = ncol(x)
  if (!is.null(p)) {
    if (is.numeric(p)) {
      if (nTraits > 1) {
        stop("When x contains more than one column, you must supply a list of probabilities! See examples.")
      }
      p = list(p)
    }
    if (length(p) != nTraits) {
      stop("You must supply probabilities for all traits in x !")
    }
    if (length(mean) != nTraits) {
      stop("You must supply means for all traits in x !")
    }
    if (length(var) != nTraits) {
      stop("You must supply variances for all traits in x !")
    }
    threshold = p
    for (trt in 1:nTraits) {
      if (!is.null(p[[trt]])) {
        pSum = sum(p[[trt]])
        if (pSum > 1) {
          stop("Probabilities for trait ", trt, " sum to more than 1!")
        } else if (pSum < 1) {
          warning("Probabilities do not sum to 1 for trait ", trt,
            "! Creating one more category!")
          p[[trt]] = c(p[[trt]], 1 - pSum)
        }
        tmp = qnorm(p = cumsum(p[[trt]]), mean = mean[trt], sd = sqrt(var[trt]))
        if (!(-Inf %in% tmp)) {
          tmp = c(-Inf, tmp)
        }
        if (!(Inf %in% tmp)) {
          tmp = c(tmp, Inf)
        }
        threshold[[trt]] = tmp
      }
    }
  }
  if (is.numeric(threshold)) {
    if (nTraits > 1) {
      stop("When x contains more than one column, you must supply a list of thresholds! See examples.")
    }
    threshold = list(threshold)
  }
  if (length(threshold) != nTraits) {
    stop("You must supply thresholds for all traits in x!")
  }
  for (trt in 1:nTraits) {
    if (!is.null(threshold[[trt]])) {
      x[, trt] = as.numeric(cut(x = x[, trt], breaks = threshold[[trt]],
                                include.lowest = include.lowest, right = right))
    }
  }
  return(x)
}

#' @title Convert a normal (Gaussian) trait to a count (Poisson) trait
#' @param x matrix, values for one or more traits (if not a matrix,
#'   we cast to a matrix).
#' @param meanLogShift \code{NULL}, numeric or list, additional additive shift(s)
#'   on the latent log scale; when \code{NULL} a shift of 0 is assumed, when
#'   numeric shifts for all traits in \code{x} must be provided, and when list
#'   shifts for all traits in \code{x} must be provided with the possibility to
#'   pass a \code{NULL} list node to skip the conversion for a trait
#'   (see examples).
#' @return matrix of values with some traits recoded as counts
#' @details If input trait is normal (Gaussian) then this function generates a
#'   count trait by sampling from the Poisson generalized linear model.
#'   As such, this function's output is stochastic.
#' 
#'   Specifically, it generates \code{y | x ~ Poisson(lambda)} with
#'   \code{lambda = exp(meanLogShift + x)}. If the supplied latent values
#'   \code{x} have mean \code{mu} and variance \code{sigma2}, then the
#'   marginal expected value of the counts is
#'   \code{exp(meanLogShift + mu + sigma2/2)}, which is the same mean
#'   formula as in the log-normal case. Therefore, to target an expected
#'   count mean \code{M}, set \code{meanLogShift = log(M) - mu - sigma2/2}.
#'   The marginal variance differs from the log-normal case and is
#'   \code{E(y) + E(y)^2 * (exp(sigma2) - 1)}, because Poisson sampling
#'   adds the extra \code{E(y)} term on top of the latent-scale heterogeneity.
#'   Consequently, latent variance can be used to correct the expected
#'   mean and to induce overdispersion, but it does not fully determine
#'   the observed variance. If \code{x} already contains an added Gaussian
#'   residual term, that latent variance contributes to the overdispersion as well.
#' 
#'   The name \code{meanLogShift} is used to emphasize that this argument
#'   is an additional shift applied during transformation, not the primary
#'   way to set the latent trait mean. In normal AlphaSimR workflow, the
#'   latent mean and variance are usually already set via
#'   \code{SP$addTrait*(..., mean = ..., var = ...)} in founding population
#'   and \code{SP$setVarE}. Hence, \code{meanLogShift} should be left at
#'   its default unless an extra transformation-specific shift is needed.
#'   One example is to control the mean of the observed values as shown below.
#'   This value should be established at the start of simulation and kept
#'   constant for most use cases.
#' @examples
#' #Simulate a founder pop, set latent trait parameters, and create a population
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' trtMeanLog = c(0, 0)
#' trtVarGLog = c(1, 2)
#' SP$addTraitA(nQtlPerChr = 10, mean = trtMeanLog, var = trtVarGLog,
#'              corA = matrix(data = c(1.0, 0.6,
#'                                     0.6, 1.0), ncol = 2))
#' trtVarELog = c(1, 1)
#' trtVarPLog = trtVarGLog + trtVarELog
#' SP$setVarE(varE = trtVarELog)
#' pop = newPop(founderPop)
#' popLarge = randCross(pop, nCrosses = 1000)
#' pheno(pop)
#'
#' meanVarFun = function(x) list(mean = mean(x), var = var(x))
#'
#' #Convert a single input trait
#' (countGv = asPoisson(gv(pop)[, 1]))
#' meanVarFun(countGv)
#'
#' #Demonstrate meanLogShift argument
#' #For Y|x ~ Poisson(exp(meanLogShift + x)), the expected count is
#' #E(Y)=exp(meanLogShift + mean(x) + var(x)/2). Trait 1 is centered on the
#' #latent scale, so we use -variance/2 to target expected mean 1.
#' (countGvMean1 = asPoisson(gv(pop)[, 1], meanLogShift = -trtVarGLog[1]/2))
#' meanVarFun(countGvMean1)
#'
#' #If x already contains Gaussian residual variance, use the corresponding
#' #latent variance in the same correction.
#' (countPhenoMean1 = asPoisson(pheno(pop)[, 1], meanLogShift = -trtVarPLog[1]/2))
#' meanVarFun(countPhenoMean1)
#'
#' #Large population example to inspect means and variances
#' gvLarge = gv(popLarge)[, 1]
#' phenoLarge = pheno(popLarge)[, 1]
#' tmp = cbind(
#'   gv = gvLarge,
#'   countGv = c(asPoisson(gvLarge)),
#'   countGvMean1 = c(asPoisson(gvLarge, meanLogShift = -trtVarGLog[1]/2)),
#'   pheno = phenoLarge,
#'   countPhenoMean1 = c(asPoisson(phenoLarge, meanLogShift = -trtVarPLog[1]/2))
#' )
#' (tmp2 = apply(X = tmp, MARGIN = 2, FUN = meanVarFun))
#' hist(tmp[, "gv"], main = paste0("Mean: ", v = tmp2$gv$mean,
#'   ", Var: ", v = tmp2$gv$var))
#' abline(v = tmp2$gv$mean, col = "red")
#' hist(tmp[, "countGv"], main = paste0("Mean: ", v = tmp2$countGv$mean,
#'   ", Var: ", v = tmp2$countGv$var))
#' abline(v = tmp2$countGv$mean, col = "red")
#' hist(tmp[, "countGvMean1"], main = paste0("Mean: ", v = tmp2$countGvMean1$mean,
#'   ", Var: ", v = tmp2$countGvMean1$var))
#' abline(v = tmp2$countGvMean1$mean, col = "red")
#' hist(tmp[, "pheno"], main = paste0("Mean: ", v = tmp2$pheno$mean,
#' ", Var: ", v = tmp2$pheno$var))
#' abline(v = tmp2$pheno$mean, col = "red")
#' hist(tmp[, "countPhenoMean1"], main = paste0("Mean: ", v = tmp2$countPhenoMean1$mean,
#'   ", Var: ", v = tmp2$countPhenoMean1$var))
#' abline(v = tmp2$countPhenoMean1$mean, col = "red")
#'
#' #Convert multiple input traits
#' try(asPoisson(pheno(pop), meanLogShift = 0))
#' asPoisson(gv(pop), meanLogShift = c(-trtVarGLog[1]/2, -trtVarGLog[2]/2))
#' asPoisson(pheno(pop), meanLogShift = list(-trtVarPLog[1]/2, NULL))
#'
#' #Store the recoded trait manually
#' pheno(pop)
#' pop@pheno[, 1] = asPoisson(pheno(pop)[, 1])
#' pheno(pop)
#' 
#' #Apply and store the transformation automatically via SimParam$finalizePop()
#' finalizePopDefault = SP$finalizePop
#' SP$finalizePop = function(pop, simParam = SP, ...) {
#'   pop@pheno[, 1] = asPoisson(pheno(pop)[, 1])
#'   return(pop)
#' }
#' pop = newPop(founderPop)
#' pheno(pop)
#' 
#' #Apply and store the transformation automatically via SimParam$finalizePheno()
#' SP$finalizePop = finalizePopDefault
#' SP$finalizePheno = function(pheno, pop, simParam = SP, ...) {
#'   pheno[, 1] = asPoisson(pheno[, 1])
#'   return(pheno)
#' }
#' pop = newPop(founderPop)
#' pheno(pop)
#' @export
asPoisson <- function(x, meanLogShift = NULL) {
  if (!is.matrix(x)) {
    x = as.matrix(x)
  }
  nTraits = ncol(x)
  if (is.null(meanLogShift)) {
    meanLogShift = rep(x = 0, times = nTraits)
  }
  if (is.numeric(meanLogShift)) {
    if (length(meanLogShift) != nTraits) {
      stop("You must supply meanLogShift for all traits in x!")
    }
    for (trt in 1:nTraits) {
      x[, trt] = rpois(n = nrow(x), lambda = exp(meanLogShift[trt] + x[, trt]))
    }
  } else if (is.list(meanLogShift)) {
    if (length(meanLogShift) != nTraits) {
      stop("You must supply meanLogShift for all traits in x!")
    }
    for (trt in 1:nTraits) {
      if (!is.null(meanLogShift[[trt]])) {
        x[, trt] = rpois(n = nrow(x), lambda = exp(meanLogShift[[trt]] + x[, trt]))
      }
    }
  } else {
    stop("meanLogShift must be NULL, numeric, or list!")
  }
  return(x)
}
