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

#' @title Sumarize genetic parameters
#' 
#' @description 
#' Calculates genetic and genic additive and dominance variances 
#' for an object of \code{\link{Pop-class}}
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param indValues should breeding values, dominance deviations 
#' and allele subsitution effects be returned
#' @param simParam an object of \code{\link{SimParam}}
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
genParam = function(pop,indValues=FALSE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  stopifnot(class(pop)=="Pop")
  bv=NULL
  dd=NULL
  genicVarA=NULL
  genicVarD=NULL
  alpha=list()
  #Loop through bv and dd calculations
  for(i in 1:simParam$nTraits){
    trait = simParam$traits[[i]]
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
  if(indValues){
    output$bv = bv
    output$dd = dd
    output$alpha = alpha
  }
  return(output)
}

#' @title Additive variance
#' 
#' @description Returns additive variance for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
varA = function(pop,simParam=NULL){
  genParam(pop,FALSE,simParam=simParam)$varA
}

#' @title Dominance variance
#' 
#' @description Returns dominance variance for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
varD = function(pop,simParam=NULL){
  genParam(pop,FALSE,simParam=simParam)$varD
}

#' @title Breeding value
#' 
#' @description Returns breeding values for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
bv = function(pop,simParam=NULL){
  genParam(pop,TRUE,simParam=simParam)$bv
}

#' @title Dominance deviations
#' 
#' @description Returns dominance deviations for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
dd = function(pop,simParam=NULL){
  genParam(pop,TRUE,simParam=simParam)$dd
}

#' @title Additive genic variance
#' 
#' @description Returns additive genic variance for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
genicVarA = function(pop,simParam=NULL){
  genParam(pop,FALSE,simParam=simParam)$genicVarA
}

#' @title Dominance genic variance
#' 
#' @description Returns dominance genic variance for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
genicVarD = function(pop,simParam=NULL){
  genParam(pop,FALSE,simParam=simParam)$genicVarD
}

#' @title Total genic variance
#' 
#' @description Returns total genic variance for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
genicVarG = function(pop,simParam=NULL){
  genParam(pop,FALSE,simParam=simParam)$genicVarG
}

#' @title Additive variance
#' 
#' @description Returns additive variance for all traits
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @export
varA = function(pop,simParam=NULL){
  genParam(pop,FALSE,simParam=simParam)$varA
}

#' @title Genetic value
#' 
#' @description A wrapper for accessing the gv slot
#' 
#' @param pop a \code{\link{Pop-class}} or similar object
#' 
#' @export
gv = function(pop){
  pop@gv
}

#' @title Phenotype
#' 
#' @description A wrapper for accessing the pheno slot
#' 
#' @param pop a \code{\link{Pop-class}} or similar object
#' 
#' @export
pheno = function(pop){
  pop@pheno
}

#' @title Estimated breeding value
#' 
#' @description A wrapper for accessing the ebv slot
#' 
#' @param pop a \code{\link{Pop-class}} or similar object
#' 
#' @export
ebv = function(pop){
  pop@ebv
}

#' @title Number of individuals
#' 
#' @description A wrapper for accessing the nInd slot
#' 
#' @param pop a \code{\link{Pop-class}} or similar object
#' 
#' @export
nInd = function(pop){
  pop@nInd
}
