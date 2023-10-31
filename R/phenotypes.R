#Adds random error to a matrix of genetic values
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

#See setPheno documentation
calcPheno = function(pop, varE, reps, p, traits, simParam){
  nTraits = length(traits)

  if(nTraits==0L){
    return(pop@pheno)
  }

  gv = pop@gv
  for(i in 1:nTraits){
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
    for(i in 1:nTraits){
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
    for(i in 1:nTraits){
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


