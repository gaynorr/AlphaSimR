#' @title Mean genetic values
#' 
#' @description Returns the mean genetic values for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#' 
#' @export
meanG = function(pop){
  colMeans(pop@gv)
}

#' @title Mean phenotypic values
#' 
#' @description Returns the mean phenotypic values for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#' 
#' @export
meanP = function(pop){
  colMeans(pop@pheno)
}

#' @title Total genetic variance
#' 
#' @description Returns total genetic variance for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#' 
#' @export
varG = function(pop){
  popVar(pop@gv)
}

#' @title Phenotypic variance
#' 
#' @description Returns phenotypic variance for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#' 
#' @export
varP = function(pop){
  popVar(pop@pheno)
}

#' @title Additive and dominance variances
#' 
#' @description 
#' Calculates genetic and genic additive and dominance variances 
#' for an object of \code{\link{Pop-class}}
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param retGenParam should genetic values for breeding values, 
#' dominance deviations and allele effects be returned
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return
#' \describe{
#' \item{varA}{an nTrait by nTrait matrix of additive genetic variances}
#' \item{varD}{an nTrait by nTrait matrix of dominance genetic variances}
#' \item{varG}{an nTrait by nTrait matrix of total genetic variances}
#' \item{genicVarA}{an nTrait vector of additive genic variances}
#' \item{genicVarD}{an nTrait vector of dominance genic variances}
#' \item{genicVarG}{an nTrait vector of total genic variances}
#' \item{bv}{an nInd by nTrait matrix of breeding values with dimensions nInd by nTraits}
#' \item{dd}{an nInd by nTrait matrix of dominance deviations with dimensions nInd by nTraits}
#' \item{alpha}{an nTrait list of allele subsitution effects}
#' }
#' 
#' @export
varAD = function(pop,retGenParam=FALSE,simParam=SIMPARAM){
  stopifnot(class(pop)=="Pop")
  bv=NULL
  dd=NULL
  genicVarA=NULL
  genicVarD=NULL
  alpha=list()
  #Loop through bv and dd calculations
  for(i in 1:simParam@nTraits){
    trait = simParam@traits[[i]]
    tmp = calcGenParam(trait,pop)
    genicVarA = c(genicVarA,tmp$genicVarA)
    genicVarD = c(genicVarD,tmp$genicVarD)
    bv = cbind(bv,tmp$bv)
    dd = cbind(dd,tmp$dd)
    alpha[[i]] = tmp$alpha
  }
  output = list(varA=popVar(bv),
                varD=popVar(dd),
                varG=varG(pop),
                genicVarA=genicVarA,
                genicVarD=genicVarD,
                genicVarG=genicVarA+genicVarD)
  if(retGenParam){
    output$bv = bv
    output$dd = dd
    output$alpha = alpha
  }
  return(output)
}

