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
#' Calculates additive and dominance variances for an object of 
#' \code{\link{Pop-class}}
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param retGenParam should genetic values for breeding values, 
#' dominance deviations and allele effects be returned
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return
#' \describe{
#' \item{varA}{a vector of additive genetic variances}
#' \item{varD}{a vector of dominance genetic variances}
#' \item{bv}{a matrix of breeding values with dimensions nInd by nTraits}
#' \item{dd}{a matrix of dominance deviations with dimensions nInd by nTraits}
#' \item{alpha}{a list of allele subsitution effects for each trait}
#' }
#' 
#' @export
varAD = function(pop,retGenParam=FALSE,simParam=SIMPARAM){
  stopifnot(class(pop)=="Pop")
  bv=NULL
  dd=NULL
  alpha=list()
  #Loop through bv and dd calculations
  for(i in 1:simParam@nTraits){
    trait = simParam@traits[[i]]
    if(class(trait)=="TraitA"){
      tmp = list()
      tmp$bv = scale(pop@gv[,i],scale=F)
      tmp$dd = rep(0,pop@nInd)
      tmp$alpha = trait@addEff
    }else if(class(trait)=="TraitAD"){
      tmp = calcGenParam(trait,pop)
    }else if(class(trait)=="TraitAG"){
      
    }else if(class(trait)=="TraitADG"){
      
    }else{
      stop("No method for trait class",class(simParam@traits[[i]]))
    }
    bv = cbind(bv,tmp$bv)
    dd = cbind(dd,tmp$dd)
    alpha[[i]] = tmp$alpha
  }
  output = list(varA=popVar(bv),
                varD=popVar(dd))
  if(retGenParam){
    output$bv = bv
    output$dd = dd
    output$alpha = alpha
  }
  return(output)
}
