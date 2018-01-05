addError = function(gv,varE,reps=1){
  nTraits = ncol(gv)
  nInd = nrow(gv)
  if(is.matrix(varE)){
    stopifnot(isSymmetric(varE),
              ncol(varE)==nTraits)
    if(any(diag(varE)==0)){
      zeros = which(diag(varE)==0)
      diag(varE)[zeros] = 1
      error = matrix(rnorm(nInd*nTraits),
                     ncol=nTraits)%*%chol(varE)
      error[,zeros] = 0
    }else{
      error = matrix(rnorm(nInd*nTraits),
                     ncol=nTraits)%*%chol(varE)
    }
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

#' @title Calculate phenotypes
#' 
#' @description 
#' Calculates phenotypes for all traits by adding random error 
#' from a multivariate normal distribution. This function is called 
#' by \code{\link{setPheno}}.
#' 
#' @param pop an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param varE error variances for phenotype. A vector of length 
#' nTraits for independent error or a square matrix of dimensions 
#' nTraits for correlated errors. If NULL, value in simParam is used.
#' @param reps number of replications for phenotype. See details.
#' @param p the p-value for the environmental covariate 
#' used by GxE traits.
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @details
#' The reps parameter is for convient representation of replicated data. 
#' It is intended to represent replicated yield trials in plant 
#' breeding programs. In this case, varE is set to the plot error and 
#' reps is set to the number of plots per entry. The resulting phenotype 
#' represents entry means.
#' 
#' @return Returns a matrix of nInd by nTrait phenotypes
#' 
#' @export
calcPheno = function(pop,varE=NULL,reps=1,p=0.5,
                     simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  validObject(pop)
  if(length(p)==1){
    p = rep(p,simParam$nTraits)
  }
  stopifnot(length(p)==simParam$nTraits)
  gv = pop@gv
  if(is.null(varE)){
    varE = simParam$varE
  }
  for(i in 1:simParam$nTraits){
    traitClass = class(simParam$traits[[i]])
    if(traitClass=="TraitAG" | traitClass=="TraitADG"){
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
#' @param varE error variances for phenotype. A vector of length 
#' nTraits for independent error or a square matrix of dimensions 
#' nTraits for correlated errors. If NULL, value in simParam is used.
#' @param reps number of replications for phenotype. See details.
#' @param p the p-value for the environmental covariate 
#' used by GxE traits.
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @details
#' The reps parameter is for convient representation of replicated data. 
#' It is intended to represent replicated yield trials in plant 
#' breeding programs. In this case, varE is set to the plot error and 
#' reps is set to the number of plots per entry. The resulting phenotype 
#' represents entry means.
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' 
#' @export
setPheno = function(pop,varE=NULL,reps=1,p=0.5,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  pop@pheno = calcPheno(pop=pop,varE=varE,reps=reps,p=p,
                        simParam=simParam)
  validObject(pop)
  return(pop)
}
