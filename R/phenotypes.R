addError = function(gv,varE){
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
  pheno = gv + MASS::mvrnorm(nInd,
                             mu=rep(0,nTraits),
                             Sigma=varE)
  return(pheno)
}

#' @title Calculate phenotypes
#' 
#' @description 
#' Calculates phenotypes for all traits by adding random error 
#' from a multivariate normal distribution.
#' 
#' @param pop an object of \code{\link{Pop-class}} or 
#' \code{\link{HybridPop-class}}
#' @param varE error variances for phenotype. A vector of length 
#' nTraits for independent error or a square matrix of dimensions 
#' nTraits for correlated errors.
#' @param reps number of replications for phenotype. This value 
#' shrinks error variance by setting varE to varE/reps.
#' @param w the environmental covariate used by GxE traits. If pop 
#' is \code{\link{HybridPop-class}} an error will be returned for any 
#' value other than 0.
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns a matrix of phenotypes with dimensions nInd by 
#' nTraits.
#' 
#' @export
calcPheno = function(pop,varE,reps=1,w=0,simParam=SIMPARAM){
  validObject(pop)
  gv = pop@gv
  if(class(pop)=="HybridPop"){
    stopifnot(w==0)
  }else{
    stopifnot(class(pop)=="Pop")
    for(i in 1:simParam@nTraits){
      traitClass = class(simParam@traits[[i]])
      if(traitClass=="TraitAG" | traitClass=="TraitADG"){
        gv[,i] = getGv(simParam@traits[[i]],pop=pop,w=w)
      }
    }
  }
  pheno = addError(gv,varE/reps)
  return(pheno)
}