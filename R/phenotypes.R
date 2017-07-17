addError = function(gv,varE,reps=1){
  nTraits = ncol(gv)
  nInd = nrow(gv)
  if(is.matrix(varE)){
    stopifnot(nrow(varE)==nTraits,
              ncol(varE)==nTraits)
  }else{
    stopifnot(length(varE)==nTraits)
    if(length(varE)==1){
      varE = matrix(varE)
    }else{
      varE = diag(varE)
    }
  }
  error = matrix(mvrnorm(nInd,
                         mu=rep(0,nTraits),
                        Sigma=varE),ncol=nTraits)
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
#' nTraits for correlated errors.
#' @param reps number of replications for phenotype. See details.
#' @param w the environmental covariate used by GxE traits. If pop 
#' is \code{\link{HybridPop-class}} an error will be returned for any 
#' value other than 0.5.
#' @param simParam an object of \code{\link{SimParam-class}}
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
calcPheno = function(pop,varE,reps=1,w=0.5,simParam=SIMPARAM){
  validObject(pop)
  gv = pop@gv
  if(class(pop)=="HybridPop"){
    stopifnot(w==0.5)
  }else{
    stopifnot(class(pop)=="Pop")
    for(i in 1:simParam@nTraits){
      traitClass = class(simParam@traits[[i]])
      if(traitClass=="TraitAG" | traitClass=="TraitADG"){
        gv[,i] = getGv(simParam@traits[[i]],pop=pop,w=w)
      }
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
#' nTraits for correlated errors.
#' @param reps number of replications for phenotype. See details.
#' @param w the environmental covariate used by GxE traits. If pop 
#' is \code{\link{HybridPop-class}} an error will be returned for any 
#' value other than 0.5.
#' @param simParam an object of \code{\link{SimParam-class}}
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
setPheno = function(pop,varE,reps=1,w=0.5,simParam=SIMPARAM){
  pop@pheno = calcPheno(pop=pop,varE=varE,reps=reps,w=w,
                        simParam=simParam)
  validObject(pop)
  return(pop)
}
