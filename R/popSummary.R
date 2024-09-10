
#' @title Mean genetic values
#'
#' @description Returns the mean genetic values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' meanG(pop)
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
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' meanP(pop)
#'
#' @export
meanP = function(pop){
  colMeans(pop@pheno)
}

#' @title Mean estimated breeding values
#'
#' @description Returns the mean estimated breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' trtH2 = 0.5
#' SP$setVarE(h2=trtH2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@ebv = trtH2 * (pop@pheno - meanP(pop)) #ind performance based EBV
#' meanEBV(pop)
#'
#' @export
meanEBV = function(pop){
  colMeans(pop@ebv)
}

#' @title Total genetic variance
#'
#' @description Returns total genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop)
#'
#' @export
varG = function(pop){
  G = popVar(pop@gv)
  rownames(G) = colnames(G) = colnames(pop@gv)
  return(G)
}

#' @title Phenotypic variance
#'
#' @description Returns phenotypic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varP(pop)
#'
#' @export
varP = function(pop){
  P = popVar(pop@pheno)
  rownames(P) = colnames(P) = colnames(pop@pheno)
  return(P)
}

#' @title Variance of estimated breeding values
#'
#' @description Returns variance of estimated breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' trtH2 = 0.5
#' SP$setVarE(h2=trtH2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@ebv = trtH2 * (pop@pheno - meanP(pop)) #ind performance based EBV
#' varA(pop)
#' varEBV(pop)
#'
#' @export
varEBV = function(pop){
  ebv = popVar(pop@ebv)
  rownames(ebv) = colnames(ebv) = colnames(pop@ebv)
  return(ebv)
}

#' @title Sumarize genetic parameters
#'
#' @description
#' Calculates genetic and genic additive and dominance variances
#' for an object of \code{\link{Pop-class}}
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return
#' \describe{
#' \item{varA}{an nTrait by nTrait matrix of additive genetic variances}
#' \item{varD}{an nTrait by nTrait matrix of dominance genetic variances}
#' \item{varAA}{an nTrait by nTrait matrix of additive-by-additive genetic variances}
#' \item{varG}{an nTrait by nTrait matrix of total genetic variances}
#' \item{genicVarA}{an nTrait vector of additive genic variances}
#' \item{genicVarD}{an nTrait vector of dominance genic variances}
#' \item{genicVarAA}{an nTrait vector of additive-by-additive genic variances}
#' \item{genicVarG}{an nTrait vector of total genic variances}
#' \item{covA_HW}{an nTrait vector of additive covariances due to non-random mating}
#' \item{covD_HW}{an nTrait vector of dominance covariances due to non-random mating}
#' \item{covAA_HW}{an nTrait vector of additive-by-additive covariances due to non-random mating}
#' \item{covG_HW}{an nTrait vector of total genic covariances due to non-random mating}
#' \item{covA_L}{an nTrait vector of additive covariances due to linkage disequilibrium}
#' \item{covD_L}{an nTrait vector of dominance covariances due to linkage disequilibrium}
#' \item{covAA_L}{an nTrait vector of additive-by-additive covariances due to linkage disequilibrium}
#' \item{covAD_L}{an nTrait vector of additive by dominance covariances due to linkage disequilibrium}
#' \item{covAAA_L}{an nTrait vector of additive by additive-by-additive covariances due to linkage disequilibrium}
#' \item{covDAA_L}{an nTrait vector of dominance by additive-by-additive covariances due to linkage disequilibrium}
#' \item{covG_L}{an nTrait vector of total genic covariances due to linkage disequilibrium}
#' \item{mu}{an nTrait vector of trait means}
#' \item{mu_HW}{an nTrait vector of expected trait means under random mating}
#' \item{gv}{a matrix of genetic values with dimensions nInd by nTraits}
#' \item{bv}{a matrix of breeding values with dimensions nInd by nTraits}
#' \item{dd}{a matrix of dominance deviations with dimensions nInd by nTraits}
#' \item{aa}{a matrix of additive-by-additive epistatic deviations with dimensions nInd by nTraits}
#' \item{gv_mu}{an nTrait vector of intercepts with dimensions nInd by nTraits}
#' \item{gv_a}{a matrix of additive genetic values with dimensions nInd by nTraits}
#' \item{gv_d}{a matrix of dominance genetic values with dimensions nInd by nTraits}
#' \item{gv_aa}{a matrix of additive-by-additive genetic values with dimensions nInd by nTraits}
#' \item{alpha}{a list of average allele subsitution effects with length nTraits}
#' \item{alpha_HW}{a list of average allele subsitution effects at Hardy-Weinberg equilibrium with length nTraits}
#' }
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' ans = genParam(pop, simParam=SP)
#'
#' @export
genParam = function(pop,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  nInd = nInd(pop)
  nTraits = simParam$nTraits
  traitNames = simParam$traitNames

  # Blank nInd x nTrait matrices
  gv = matrix(NA_real_, nrow=nInd, ncol=nTraits)
  colnames(gv) = traitNames
  bv = dd = aa = gv_a = gv_d = gv_aa = gv

  # Blank nTrait vectors
  genicVarA = rep(NA_real_, nTraits)
  names(genicVarA) = traitNames
  genicVarD = genicVarAA = covA_HW = covD_HW = covAA_HW =
    covG_HW = mu = mu_HW = gv_mu = covAAA_L = covDAA_L =
    covAD_L = genicVarA

  # Average effect of an allele substitution
  alpha = vector("list", length=nTraits)
  names(alpha) = traitNames
  alpha_HW = alpha

  #Loop through trait calculations
  for(i in seq_len(nTraits)){
    trait = simParam$traits[[i]]
    tmp = calcGenParam(trait,pop,simParam$nThreads)
    genicVarA[i] = tmp$genicVarA2
    covA_HW[i] = tmp$genicVarA-tmp$genicVarA2
    gv[,i] = tmp$gv
    bv[,i] = tmp$bv
    mu[i] = tmp$mu
    mu_HW[i] = tmp$mu_HWE
    gv_a[,i] = tmp$gv_a
    gv_mu[i] = tmp$gv_mu
    if(.hasSlot(trait,"domEff")){
      genicVarD[i] = tmp$genicVarD2
      covD_HW[i] = tmp$genicVarD-tmp$genicVarD2
      dd[,i] = tmp$dd
      gv_d[,i] = tmp$gv_d
    }else{
      genicVarD[i] = 0
      covD_HW[i] = 0
      dd[,i] = rep(0,pop@nInd)
      gv_d[,i] = rep(0,pop@nInd)
    }
    if(.hasSlot(trait,"epiEff")){
      genicVarAA[i] = tmp$genicVarAA2
      covAA_HW[i] = tmp$genicVarAA-tmp$genicVarAA2
      aa[,i] = tmp$aa
      gv_aa[,i] = tmp$gv_aa
    }else{
      genicVarAA[i] = 0
      covAA_HW[i] = 0
      aa[,i] = rep(0,pop@nInd)
      gv_aa[,i] = rep(0,pop@nInd)
    }
    if(nInd==1){
      covAD_L[i] = 0
      covAAA_L[i] = 0
      covDAA_L[i] = 0
    } else {
      covAD_L[i] = popVar(cbind(bv[,i],dd[,i]))[1,2]
      covAAA_L[i] = popVar(cbind(bv[,i],aa[,i]))[1,2]
      covDAA_L[i] = popVar(cbind(dd[,i],aa[,i]))[1,2]
    }
    alpha[[i]] = tmp$alpha
    alpha_HW[[i]] = tmp$alpha_HW
  }

  varA = popVar(bv)
  rownames(varA) = colnames(varA) = traitNames

  varD = popVar(dd)
  rownames(varD) = colnames(varD) = traitNames

  varAA = popVar(aa)
  rownames(varAA) = colnames(varAA) = traitNames

  varG = popVar(gv)
  rownames(varG) = colnames(varG) = traitNames

  genicVarG = genicVarA + genicVarD + genicVarAA
  covG_HW = covA_HW + covD_HW + covAA_HW

  output = list(varA=varA,
                varD=varD,
                varAA=varAA,
                varG=varG,
                genicVarA=genicVarA,
                genicVarD=genicVarD,
                genicVarAA=genicVarAA,
                genicVarG=genicVarG,
                covA_HW=covA_HW,
                covD_HW=covD_HW,
                covAA_HW=covAA_HW,
                covG_HW=covG_HW,
                covA_L=diag(varA)-genicVarA-covA_HW,
                covD_L=diag(varD)-genicVarD-covD_HW,
                covAA_L=diag(varAA)-genicVarAA-covAA_HW,
                covAD_L=covAD_L,
                covAAA_L=covAAA_L,
                covDAA_L=covDAA_L,
                covG_L=diag(varG)-genicVarG-covG_HW,
                mu=mu,
                mu_HW=mu_HW,
                gv=gv,
                bv=bv,
                dd=dd,
                aa=aa,
                gv_mu=gv_mu,
                gv_a=gv_a,
                gv_d=gv_d,
                gv_aa=gv_aa,
                alpha=alpha,
                alpha_HW=alpha_HW)
  return(output)
}

#' @title Additive variance
#'
#' @description Returns additive variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varA(pop, simParam=SP)
#'
#' @export
varA = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$varA
}

#' @title Dominance variance
#'
#' @description Returns dominance variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varD(pop, simParam=SP)
#'
#' @export
varD = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$varD
}

#' @title Additive-by-additive epistatic variance
#'
#' @description Returns additive-by-additive epistatic
#' variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varAA(pop, simParam=SP)
#'
#' @export
varAA = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$varAA
}

#' @title Breeding value
#'
#' @description Returns breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' bv(pop, simParam=SP)
#'
#' @export
bv = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$bv
}

#' @title Dominance deviations
#'
#' @description Returns dominance deviations for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' dd(pop, simParam=SP)
#'
#' @export
dd = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$dd
}

#' @title Additive-by-additive epistatic deviations
#'
#' @description Returns additive-by-additive epistatic
#' deviations for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' aa(pop, simParam=SP)
#'
#' @export
aa = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$aa
}

#' @title Additive genic variance
#'
#' @description Returns additive genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#'
#' @export
genicVarA = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$genicVarA
}

#' @title Dominance genic variance
#'
#' @description Returns dominance genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarD(pop, simParam=SP)
#'
#' @export
genicVarD = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$genicVarD
}

#' @title Additive-by-additive genic variance
#'
#' @description Returns additive-by-additive epistatic
#' genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarAA(pop, simParam=SP)
#'
#' @export
genicVarAA = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$genicVarAA
}

#' @title Total genic variance
#'
#' @description Returns total genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#'
#' @export
genicVarG = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$genicVarG
}

#' @title Genetic value
#'
#' @description A wrapper for accessing the gv slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' gv(pop)
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
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pheno(pop)
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
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@ebv = matrix(rnorm(pop@nInd), nrow=pop@nInd, ncol=1)
#' ebv(pop)
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
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' nInd(pop)
#'
#' @export
nInd = function(pop){
  pop@nInd
}
