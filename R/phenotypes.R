#Adds random error to a matrix of genetic values
addError = function(gv,varE,reps=1){
  nTraits = ncol(gv)
  nInd = nrow(gv)
  if(is.matrix(varE)){
    stopifnot(isSymmetric(varE),
              ncol(varE)==nTraits)
    error = matrix(rnorm(nInd*nTraits),
                   ncol=nTraits)%*%rotMat(varE)
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
  error = error/sqrt(rep(reps,nrow(error)))
  pheno = gv + error
  return(pheno)
}

#See setPheno documentation
calcPheno = function(pop,varE,reps,p,simParam){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(simParam$nTraits == 0L){
    return(matrix(NA_real_,
                  nrow=pop@nInd,
                  ncol=0L))
  }
  if(is.null(p)){
    p = rep(runif(1), simParam$nTraits)
  }else if(length(p)==1){
    p = rep(p,simParam$nTraits)
  }
  stopifnot(length(p)==simParam$nTraits)
  gv = pop@gv
  if(is.null(varE)){
    varE = simParam$varE
  }
  for(i in 1:simParam$nTraits){
    if(.hasSlot(simParam$traits[[i]], "envVar")){
      stdDev = sqrt(simParam$traits[[i]]@envVar)
      gv[,i] = gv[,i]+pop@gxe[[i]]*qnorm(p[i],sd=stdDev)
    }
  }
  pheno = addError(gv=gv,varE=varE,reps=reps)
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
#' @param reps number of replications for phenotype. See details.
#' @param fixEff fixed effect to assign to the population. Used 
#' by genomic selection models only.
#' @param p the p-value for the environmental covariate 
#' used by GxE traits. If NULL, a value is
#' sampled at random.
#' @param onlyPheno should only the phenotype be returned, see return
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
#' SP$addTraitA(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Add phenotype with error variance of 1
#' pop = setPheno(pop, varE=1)
#' 
#' @export
setPheno = function(pop,h2=NULL,H2=NULL,varE=NULL,reps=1,
                    fixEff=1L,p=NULL,onlyPheno=FALSE,
                    simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(class(pop)=="MegaPop"){
    stopifnot(!onlyPheno)
    pop@pops = lapply(pop@pops, setPheno, h2=h2, H2=H2,
                      varE=varE, reps=reps, fixEff=fixEff, 
                      p=p, simParam=simParam)
    return(pop)
  }
  
  # Calculate varE if using h2 or H2
  if(!is.null(h2)){
    if(length(h2)==1){
      h2 = rep(h2, simParam$nTraits)
    }
    stopifnot(length(h2)==simParam$nTraits,
              all(simParam$varG>0),
              all(simParam$varA>0))
    varE = numeric(simParam$nTraits)
    for(i in 1:length(h2)){
      tmp = simParam$varA[i]/h2[i]-simParam$varG[i]
      if(tmp<0){
        stop(paste0("h2=",h2[i]," is not possible for trait ",i))
      }
      varE[i] = tmp
    }
  }else if(!is.null(H2)){
    if(length(H2)==1){
      H2 = rep(H2, simParam$nTraits)
    }
    stopifnot(length(H2)==simParam$nTraits)
    varE = numeric(simParam$nTraits)
    for(i in 1:length(H2)){
      tmp = simParam$varG[i]/H2[i]-simParam$varG[i]
      varE[i] = tmp
    }
  }
  
  # Create phenotypes
  pheno = calcPheno(pop=pop,varE=varE,reps=reps,p=p,
                    simParam=simParam)
  if(onlyPheno){
    return(pheno)
  }
  pop@pheno = pheno
  if(is(pop,"Pop")){
    pop@fixEff = rep(as.integer(fixEff),pop@nInd)
    pop@reps = rep(as.numeric(reps),pop@nInd)
  }
  return(pop)
}


