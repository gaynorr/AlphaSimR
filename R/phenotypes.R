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
#' @param simParam simulation parameters object
#'
#' @keywords internal
calcPheno = function(pop, varE, reps, p, traits, simParam){
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
#' @param simParam an object of \code{\link{SimParam}}
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
                    traits=NULL, simParam=NULL){
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

#' @title Set population-level phenotypes
#'
#' @description
#' Compute and store population-level phenotypic summaries by applying a
#' user-supplied function to the individual-level phenotype matrix.
#' For a \code{MultiPop} object, the function is applied independently to
#' each population.
#'
#' @param x An object of class \code{\link{Pop-class}} or 
#' \code{\link{MultiPop-class}}.
#' @param FUN A function that accepts the individual-phenotype matrix
#'   from a population (e.g. \code{pop@pheno} or 
#'   \code{multiPop@pops[[i]]@miscPop$pheno}) as its first argument and 
#'   returns a numeric vector of length \code{nTraits} (one value per trait).
#'   The function will be called as \code{FUN(pop@pheno, ...)}.
#' @param force Logical. If \code{TRUE}, individual phenotypes will be
#'   (re)generated by calling \code{\link{setPheno}()} before computing
#'   the population summary. When \code{force=TRUE} the other phenotype
#'   generation arguments (see below) are passed to \code{\link{setPheno}()}.
#'   If \code{FALSE}, the function requires an existing \code{pop@pheno} 
#'   (or \code{multiPop@pops[[i]]@miscPop$pheno}) matrix and will stop with an error 
#'   if that matrix is empty.
#' @param h2 Numeric vector of individual-level narrow-sense heritabilities 
#'   (one per trait). Only passed to \code{\link{setPheno}()} when 
#'   \code{force=TRUE}, otherwise not used.
#' @param H2 Numeric vector of broad-sense heritabilities (one per trait).
#'   Only passed to \code{\link{setPheno}()} when \code{force=TRUE}, 
#'   otherwise not used.
#' @param varE Numeric vector or matrix of residual (co)variances. Only passed
#'   to \code{\link{setPheno}()} when \code{force=TRUE}, otherwise not used.
#' @param corE Optional correlation matrix for residual errors. Only passed to
#'   \code{\link{setPheno}()} when \code{force=TRUE}, otherwise not used.
#' @param reps Integer number of phenotype replications. Only passed to
#'   \code{\link{setPheno}()} when \code{force=TRUE}, otherwise not used.
#' @param fixEff Fixed effect value to assign to the population. Used by
#'   some genomic-selection models. Only passed to \code{\link{setPheno}()}
#'   when \code{force=TRUE}, otherwise not used.
#' @param p Numeric or \code{NULL}. The p-value used for environmental
#'   covariates in GxE traits. If \code{NULL}, a value is sampled
#'   randomly. Only passed to \code{\link{setPheno}()} when
#'   \code{force=TRUE}, otherwise not used.
#' @param onlyPhenoPop Logical. If \code{TRUE} the function returns the
#'   population-level phenotype matrix (see Details). If \code{FALSE}
#'   (default) the input \code{x} object is returned with its
#'   \code{pop@miscPop$pheno} or \code{multiPop@pops[[i]]@miscPop$pheno} 
#'   slot(s) updated.
#' @param traits Integer vector indicating which trait columns to set.
#'   If \code{NULL} (the default), all traits are used.
#' @param simParam A \code{\link{SimParam}} object. If \code{NULL}, the
#'   global \code{SP} is used.
#' @param ... Additional arguments passed to \code{FUN()}.
#'
#' @details
#' For a \code{MultiPop} input, \code{setPhenoPop()} is applied to every
#' population \code{pop} inside \code{multiPop@pops} and each population's 
#' result is stored in its \code{pop@miscPop$pheno} slot. 
#' When a population already has \code{pop@miscPop$pheno} 
#' present, only the columns specified by \code{traits} are overwritten; 
#' other columns are preserved.
#'
#' If \code{force=TRUE}, individual-level phenotypes are created (or
#' replaced) by calling \code{\link{setPheno}()} using the supplied
#' \code{h2}, \code{H2}, \code{varE}, \code{corE}, \code{reps}, \code{p}
#' and \code{fixEff} arguments. See \code{\link{setPheno}} documentation
#' for details about those arguments and how to use them.
#'
#' The \code{FUN} argument should return a numeric vector with a single 
#' summary value per trait. The output produced by applying 
#' \code{FUN} to \code{pop@pheno} is transposed and stored as a single-row 
#' matrix with \code{nTraits} columns.
#'
#' @return
#' If \code{onlyPhenoPop=FALSE} (default), the input object
#' (\code{Pop}, or \code{MultiPop}) is returned with
#' \code{pop@miscPop$pheno} (or \code{multiPop@pops[[i]]@miscPop$pheno}) populated 
#' or updated. If \code{onlyPhenoPop=TRUE}, a numeric matrix is returned:
#' for a single \code{Pop}, a 1-by-\code{nTraits} matrix; and
#' for a \code{MultiPop}, a \code{length(multiPop)} (rows) by \code{nTraits}
#' (columns) matrix where each row is the population-level summary for
#' one constituent population (in the same order as in \code{multiPop@pops}).
#'
#' @examples
#' # Create founder haplotypes
#' founderPop = quickHaplo(nInd = 4, nChr = 1, segSites = 10)
#'
#' # Set simulation parameters and a single additive trait
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#'
#' # Create a population and a multi-population object
#' pop = newPop(founderPop, simParam = SP)
#' multiPop = newMultiPop(pop[1:2], pop[3:4])
#'
#' # Ensure individual phenotypes exist (error variance = 1)
#' multiPop = setPheno(multiPop, varE = 1)
#' # TODO: MultiPop: Make accessory functions (pheno, gv, ...) work with MultiPop class #239
#' #       https://github.com/gaynorr/AlphaSimR/issues/239
#' # pheno(multiPop)
#' lapply(multiPop@pops, FUN = function(x) x@pheno)
#'
#' # Compute population means (one mean per trait per population)
#' multiPop = setPhenoPop(multiPop, FUN = colMeans)
#' multiPop@pops[[1]]@miscPop$pheno  # Population 1 means
#' # TODO: MultiPop: Make accessory functions (pheno, gv, ...) work with MultiPop class #239
#' #       https://github.com/gaynorr/AlphaSimR/issues/239
#' # phenoPop(multiPop) # All population means
#' sapply(multiPop@pops, FUN = function(x) x@miscPop$pheno) # All population means
#'
#' # Compute medians and return only the population-level matrix
#' medians = function(x) apply(x, 2, median)
#' setPhenoPop(multiPop, onlyPhenoPop = TRUE, FUN = medians)
#'
#' @seealso \code{\link{setPheno}}
#' @export
setPhenoPop = function(x, FUN=colMeans, force=FALSE, fixEff=1L, 
                       h2=NULL, H2=NULL, varE=NULL, corE=NULL, reps=1, p=NULL, 
                       onlyPhenoPop=FALSE, traits=NULL, simParam=NULL, ...){
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  
  # Determine which traits are selected
  if (is.null(traits)) {
    if (simParam$nTraits > 0L) {
      traits = 1:simParam$nTraits
    } else {
      traits = integer()
    }
  } else {
    traits = as.integer(traits)
    stopifnot(
      all(traits > 0L),
      all(!duplicated(traits)),
      max(traits) <= simParam$nTraits
    )
  }
  
  # Use lapply if object is a MultiPop
  if (isMultiPop(x)) {
    x@pops = lapply(x@pops, function(x) {
      x = setPhenoPop(x = x, h2=h2, H2=H2, varE=varE, corE=corE, reps=reps, 
                      fixEff=fixEff, p=p, traits=traits, onlyPhenoPop=F,
                      force=force, simParam=simParam, FUN=FUN, ...)
      return(x)
    })
  } else if (isHybridPop(x)) {
    stop("HybridPop-class objects are not supported by this function.")
  } else if (isPop(x)) {
    if (force == TRUE) {
      # Create phenotypes
      x = setPheno(pop=x, h2 = h2, H2 = H2,varE = varE, 
                   corE = corE, reps = reps, fixEff = fixEff,
                   p = p, traits = traits, simParam = simParam)
    } else if (sum(is.na(x@pheno)) == prod(dim(x@pheno))) {
      stop("The phenotypic matrix is empty. Use force=TRUE to create it.")
    }
    
    pheno = x@pheno
    phenoPop = t(FUN(pheno, ...))
    stopifnot(
      "FUN returned an object of unexpected dimensions" = dim(phenoPop) ==
        c(1, simParam$nTraits)
    )
    
    # TODO: MultiPop: Store population values into MultiPop@pops[[I]]@miscPop or create new slots in MultiPop #240
    #       https://github.com/gaynorr/AlphaSimR/issues/240
    if (is.null(x@miscPop$pheno)) {
      x@miscPop$pheno = matrix(NA_real_, ncol = x@nTraits)
    }
    x@miscPop$pheno[, traits] = phenoPop[, traits]
    colnames(x@miscPop$pheno) = colnames(x@gv)
  } else {
    stop("x must be an object of Pop or MultiPop class.")
  }
  
  if (onlyPhenoPop) {
    if (isMultiPop(x)) {
      pheno = lapply(x@pops, function(x) x@miscPop$pheno)
      pheno = do.call('rbind', pheno)
    } else {
      pheno = x@miscPop$pheno
    }
    return(pheno)
  }
  return(x)
}

#' @title Convert a normal (Gaussian) trait to an ordered categorical (threshold)
#'   trait
#' @param x matrix, values for one or more traits (if not a matrix,
#'   we cast to a matrix)
#' @param p NULL, numeric or list, when \code{NULL} the \code{threshold} argument
#'   takes precedence; when numeric, provide a vector of probabilities of
#'   categories to convert continuous values into categories for a single trait
#'   (if probabilities do not sum to 1, another category is added and a warning
#'   is raised); when list, provide a list of numeric probabilities - list node
#'   with \code{NULL} will skip conversion for a specific trait (see examples);
#'   internally \code{p} is converted to \code{threshold} hence input
#'   \code{threshold} is overwritten
#' @param mean numeric, assumed mean(s) of the normal (Gaussian) trait(s);
#'   used only when \code{p} is given
#' @param var numeric, assumed variance(s) of the normal (Gaussian) trait(s);
#'   used only when \code{p} is given
#' @param threshold NULL, numeric or list, when numeric, provide a vector of
#'   threshold values to convert continuous values into categories for a single trait
#'   (the thresholds specify left-closed and right-opened intervals [t1, t2),
#'   which can be changed with \code{include.lowest} and \code{right};
#'   ensure you add \code{-Inf} and \code{Inf} or min and max to cover the whole
#'   range of values; otherwise you will get \code{NA} values);
#'   when list, provide a list of numeric thresholds - list node with \code{NULL}
#'   will skip conversion for a specific trait (see examples)
#' @param include.lowest logical, see \code{\link{cut}}
#' @param right logical, see \code{\link{cut}}
#' @details If input trait is normal (Gaussian) then this function generates a
#'   categorical trait according to the ordered probit model.
#' @return matrix of values with some traits recorded as ordered categories
#'  in the form of 1:nC with nC being the number of categories.
#' @examples
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' trtMean = c(0, 0)
#' trtVarG = c(1, 2)
#' SP$addTraitA(nQtlPerChr = 10, mean = trtMean, var = trtVarG,
#'              corA = matrix(data = c(1.0, 0.6,
#'                                     0.6, 1.0), ncol = 2))
#' pop = newPop(founderPop)
#' trtVarE = c(1, 1)
#' trtVarP = trtVarG + trtVarE
#' pop = setPheno(pop, varE = trtVarE)
#' pheno(pop)
#'
#' #Convert a single input trait
#' asCategorical(x = pheno(pop)[, 1])
#'
#' #Demonstrate threshold argument (in units of pheno SD)
#' asCategorical(x = pheno(pop)[, 1], threshold = c(-1, 0, 1) * sqrt(trtVarP[1]))
#' asCategorical(x = pheno(pop)[, 1], threshold = c(-Inf, -1, 0, 1, Inf) * sqrt(trtVarP[1]))
#' asCategorical(x = pheno(pop)[, 1], threshold = c(-Inf, 0, Inf))
#'
#' #Demonstrate p argument
#' asCategorical(x = pheno(pop)[, 1], p = 0.5, var = trtVarP[1])
#' asCategorical(x = pheno(pop)[, 1], p = c(0.5, 0.5), var = trtVarP[1])
#' asCategorical(x = pheno(pop)[, 1], p = c(0.25, 0.5, 0.25), var = trtVarP[1])
#'
#' #Convert multiple input traits (via threshold or p argument)
#' try(asCategorical(x = pheno(pop)))
#' asCategorical(x = pheno(pop),
#'               threshold = list(c(-Inf, 0, Inf),
#'                                NULL))
#' try(asCategorical(x = pheno(pop), p = c(0.5, 0.5)))
#' asCategorical(x = pheno(pop),
#'               p = list(c(0.5, 0.5),
#'                        NULL),
#'               mean = trtMean, var = trtVarP)
#'
#' asCategorical(x = pheno(pop),
#'               threshold = list(c(-Inf, 0, Inf),
#'                                c(-Inf, -2, -1, 0, 1, 2, Inf) * sqrt(trtVarP[2])))
#' q = c(-2, -1, 0, 1, 2)
#' p = pnorm(q)
#' p = c(p[1], p[2]-p[1], p[3]-p[2], p[4]-p[3], p[5]-p[4], 1-p[5])
#' asCategorical(x = pheno(pop),
#'               p = list(c(0.5, 0.5),
#'                        p),
#'               mean = trtMean, var = trtVarP)
#' @export
asCategorical = function(x, p = NULL, mean = 0, var = 1,
                         threshold = c(-Inf, 0, Inf),
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
        if (!(pSum == 1)) { # TODO: how can we do floating point aware comparison?
          warning("Probabilities do not sum to 1 - creating one more category!")
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
    stop("You must supply thresholds for all traits in x !")
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
#'   we cast to a matrix)
#' @param TODO numeric or list, when numeric, provide a vector of TODO
#' @return matrix of values with some traits recoded as counts
#' @details If input trait is normal (Gaussian) then this function generates a
#'   count trait according to the Poisson generalised linear model.
#' @examples
#' founderPop = quickHaplo(nInd=20, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(nQtlPerChr = 10, mean = c(0, 0), var = c(1, 2),
#'              corA = matrix(data = c(1.0, 0.6,
#'                                     0.6, 1.0), ncol = 2))
#' pop = newPop(founderPop)
#' pop = setPheno(pop, varE = c(1, 1))
#' pheno(pop)
#' #Convert a single input trait
#' asCount(x = pheno(pop)[, 2])
#' asCount(x = pheno(pop)[, 2], TODO = c(-1, 0, 1))
#' asCount(x = pheno(pop)[, 2], TODO = c(-Inf, -1, 0, 1, Inf))
#' #Convert multiple input traits
#' try(asCount(x = pheno(pop)))
#' asCount(x = pheno(pop),
#'           TODO = list(NULL,
#'                       ???))
#' TODO export
# asCount = function(x, TODO = 10) {
#   if (!is.matrix(x)) {
#     x = as.matrix(x)
#   }
#   nTraits = ncol(x)
#   if (is.numeric(TODO)) {
#     if (nTraits > 1) {
#       stop("When x contains more than one column, you must supply a list of TODO! See examples.")
#     }
#     TODO = list(TODO)
#   }
#   for (trt in 1:nTraits) {
#     if (!is.null(TODO[[trt]])) {
# TODO: need to think what to do with an intercept and lambda
#       x[, trt] = round(exp(x[, trt]))
#     }
#   }
#   return(x)
# }
