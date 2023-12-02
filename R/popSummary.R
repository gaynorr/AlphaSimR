# Internal function for calculating mean EBV of populations
# Used selectPop MultiPop-class
meanEBV = function(pop){
  colMeans(pop@ebv)
}

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

#' @title Total genetic variance
#'
#' @description Returns total genetic variance (=variance of genetic values)
#' for all traits
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
#' @description Returns phenotypic variance (=variance of phenotypic values)
#' for all traits
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

#' @title Sumarize genetic parameters
#'
#' @description
#' Calculates genetic values and variances for an object of \code{\link{Pop-class}}
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return
#' \describe{
#' \item{varA}{an nTrait by nTrait matrix of additive genetic variances}
#' \item{varD}{an nTrait by nTrait matrix of dominance genetic variances}
#' \item{varI}{an nTrait by nTrait matrix of imprinting genetic variances}
#' \item{varAA}{an nTrait by nTrait matrix of additive-by-additive genetic variances}
#' \item{varG}{an nTrait by nTrait matrix of total genetic variances}
#' \item{genicVarA}{an nTrait vector of additive genic variances}
#' \item{genicVarD}{an nTrait vector of dominance genic variances}
#' \item{genicVarI}{an nTrait vector of imprinting genic variances}
#' \item{genicVarAA}{an nTrait vector of additive-by-additive genic variances}
#' \item{genicVarG}{an nTrait vector of total genic variances}
#' \item{covA_HW}{an nTrait vector of additive covariances due to non-random mating}
#' \item{covD_HW}{an nTrait vector of dominance covariances due to non-random mating}
#' \item{covI_HW}{an nTrait vector of imprinting covariances due to non-random mating}
#' \item{covAA_HW}{an nTrait vector of additive-by-additive covariances due to non-random mating}
#' \item{covG_HW}{an nTrait vector of total genic covariances due to non-random mating}
#' \item{covA_L}{an nTrait vector of additive covariances due to linkage disequilibrium}
#' \item{covD_L}{an nTrait vector of dominance covariances due to linkage disequilibrium}
#' \item{covI_L}{an nTrait vector of imprinting covariances due to linkage disequilibrium}
#' \item{covAA_L}{an nTrait vector of additive-by-additive covariances due to linkage disequilibrium}
#' \item{covAD_L}{an nTrait vector of additive by dominance covariances due to linkage disequilibrium}
#' \item{covAI_L}{an nTrait vector of additive by imprinting covariances due to linkage disequilibrium}
#' \item{covAAA_L}{an nTrait vector of additive by additive-by-additive covariances due to linkage disequilibrium}
#' \item{covDAA_L}{an nTrait vector of dominance by additive-by-additive covariances due to linkage disequilibrium}
#' \item{covIAA_L}{an nTrait vector of imprinting by additive-by-additive covariances due to linkage disequilibrium}
#' \item{covG_L}{an nTrait vector of total genic covariances due to linkage disequilibrium}
#' \item{mu}{an nTrait vector of trait means}
#' \item{mu_HW}{an nTrait vector of expected trait means under random mating}
#' \item{gv}{a matrix of genetic values with dimensions nInd by nTraits}
#' \item{bv}{a matrix of breeding values with dimensions nInd by nTraits}
#' \item{dd}{a matrix of dominance deviations with dimensions nInd by nTraits}
#' \item{idM}{a matrix of maternal imprinting deviations with dimensions nInd by nTraits
#'            (paternal imprinting deviations are \code{-idM})}
#' \item{aa}{a matrix of additive-by-additive epistatic deviations with dimensions nInd by nTraits}
#' \item{gv_mu}{an nTrait vector of intercepts with dimensions nInd by nTraits}
#' \item{gv_a}{a matrix of additive genetic values with dimensions nInd by nTraits}
#' \item{gv_d}{a matrix of dominance genetic values with dimensions nInd by nTraits}
#' \item{gv_i}{a matrix of imprinting genetic values with dimensions nInd by nTraits}
#' \item{gv_aa}{a matrix of additive-by-additive genetic values with dimensions nInd by nTraits}
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
  stopifnot(class(pop)=="Pop")
  nInd = nInd(pop)
  nTraits = simParam$nTraits
  traitNames = simParam$traitNames

  # Blank nInd x nTrait matrices
  gv = matrix(NA_real_, nrow=nInd, ncol=nTraits)
  colnames(gv) = traitNames
  bv = bvM = bvP = dd = idM = aa = gv_a = gv_d = gv_i = gv_aa = gv

  # Blank nTrait vectors
  genicVarA = rep(NA_real_, nTraits)
  names(genicVarA) = traitNames
  genicVarD = genicVarI = genicVarAA = covA_HW = covD_HW = covI_HW = covAA_HW =
    covG_HW = mu = mu_HW = gv_mu = covAAA_L = covDAA_L = covIAA_L =
    covAD_L = covAI_L = covDI_L = genicVarA

  #Loop through trait calculations
  for(i in 1:nTraits){
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
    if(.hasSlot(trait,"impEff")){
      genicVarI[i] = tmp$genicVarI2
      covI_HW[i] = tmp$genicVarI-tmp$genicVarI2
      idM[,i] = tmp$idM # taking just mat dev since pat dev = - mat dev
      gv_i[,i] = tmp$gv_i
      bvM[,i] = tmp$bvM
      bvP[,i] = tmp$bvP
    }else{
      genicVarI[i] = 0
      covI_HW[i] = 0
      idM[,i] = rep(0,pop@nInd)
      gv_i[,i] = rep(0,pop@nInd)
      bvM[,i] = rep(0,pop@nInd)
      bvP[,i] = rep(0,pop@nInd)
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
      covAI_L[i] = 0
      covDI_L[i] = 0
      covAAA_L[i] = 0
      covDAA_L[i] = 0
      covIAA_L[i] = 0
    } else {
      covAD_L[i] = popVar(cbind(bv[,i],dd[,i]))[1,2]
      covAI_L[i] = popVar(cbind(bv[,i],idM[,i]))[1,2]
      covDI_L[i] = popVar(cbind(dd[,i],idM[,i]))[1,2]
      covAAA_L[i] = popVar(cbind(bv[,i],aa[,i]))[1,2]
      covDAA_L[i] = popVar(cbind(dd[,i],aa[,i]))[1,2]
      covIAA_L[i] = popVar(cbind(idM[,i],aa[,i]))[1,2]
      # TODO: do we need covIMA and covIPA etc.?
      # TODO: do we need covIMA and covIPA etc.?
    }
  }

  # TODO: do we need varAM and varAP?
  varA = popVar(bv)
  rownames(varA) = colnames(varA) = traitNames

  varD = popVar(dd)
  rownames(varD) = colnames(varD) = traitNames

  # TODO: do we need varIM and varIP (they are the same in a "random" setting)?
  varI = popVar(idM)
  rownames(varI) = colnames(varI) = traitNames

  varAA = popVar(aa)
  rownames(varAA) = colnames(varAA) = traitNames

  varG = popVar(gv)
  rownames(varG) = colnames(varG) = traitNames

  genicVarG = genicVarA + genicVarD + genicVarI + genicVarAA
  covG_HW = covA_HW + covD_HW + covI_HW + covAA_HW

  output = list(varA=varA,
                varD=varD,
                varI=varI,
                varAA=varAA,
                varG=varG,
                genicVarA=genicVarA,
                genicVarD=genicVarD,
                genicVarI=genicVarI,
                genicVarAA=genicVarAA,
                genicVarG=genicVarG,
                covA_HW=covA_HW,
                covD_HW=covD_HW,
                covI_HW=covI_HW,
                covAA_HW=covAA_HW,
                covG_HW=covG_HW,
                covA_L=diag(varA)-genicVarA-covA_HW,
                covD_L=diag(varD)-genicVarD-covD_HW,
                covI_L=diag(varI)-genicVarI-covI_HW,
                covAA_L=diag(varAA)-genicVarAA-covAA_HW,
                covAD_L=covAD_L,
                covAI_L=covAI_L,
                covAAA_L=covAAA_L,
                covDAA_L=covDAA_L,
                covIAA_L=covIAA_L,
                covG_L=diag(varG)-genicVarG-covG_HW,
                mu=mu,
                mu_HW=mu_HW,
                gv=gv,
                bv=bv,
                bvM=bvM,
                bvP=bvP,
                dd=dd,
                idM=idM,
                aa=aa,
                gv_mu=gv_mu,
                gv_a=gv_a,
                gv_d=gv_d,
                gv_i=gv_i,
                gv_aa=gv_aa)
  return(output)
}

#' @title Additive variance
#'
#' @description Returns additive variance (=variance of breeding values)
#' for all traits
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
#' @description Returns dominance variance (=variance of dominance deviations)
#' for all traits
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

#' @title Imprinting variance
#'
#' @description Returns imprinting variance (=variance of imprinting deviations)
#' for all traits
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
#' SP$addTraitAI(10, meanID=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varI(pop, simParam=SP)
#'
#' @export
varI = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$varI
}

#' @title Additive-by-additive epistatic variance
#'
#' @description Returns additive-by-additive epistatic variance (=variance of
#' additive-by-additive epistatic deviations) for all traits
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

#' @title Breeding values
#'
#' @description Returns breeding values for all traits.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @details
#' With imprinting the output depends on sex of individuals.
#' If \code{sex == "H"} (in many plants) then mean breeding value is returned,
#' which is mean of maternal and paternal breeding value.
#' If \code{sex == "F"} (female) then maternal breeding value is returned.
#' If \code{sex == "M"} (male)   then paternal breeding value is returned.
#' If you must get maternal or paternal breeding value irrespective of the
#' sex of individuals, use \code{bvM()} and \code{bvP()} - see examples.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' # --- Additive & Dominance effects ---
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
#' # --- Additive & Imprinting effects - Bisexual individuals ---
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAI(10, meanID=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@sex
#' # Mean breeding values
#' # (individuals are bisexual)
#' (tmpBv = bv(pop, simParam=SP))
#'
#' # Maternal breeding values
#' # (as if individuals are females or for their female sexual organ)
#' (tmpMBv = bvM(pop, simParam=SP))
#' (tmpMId = tmpMBv - tmpBv)
#' idM(pop, simParam=SP)
#'
#' # Paternal breeding values
#' # (as if individuals are males or for their male sexual organ)
#' (tmpPBv = bvP(pop, simParam=SP))
#' (tmpPId = tmpPBv - tmpBv)
#' idP(pop, simParam=SP)
#'
#' # Compare the values
#' data.frame(sex = pop@sex, bv = tmpBv[, 1], bvM = tmpMBv[, 1], bvP = tmpPBv[, 1],
#'                                            idM = tmpMId[, 1], idP = tmpPId[, 1])
#'
#' # --- Male and female individuals ---
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$setSexes("yes_sys")
#' SP$addTraitAI(10, meanID=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@sex
#' # Individual breeding values
#' # (individuals are either females or males)
#' (tmpBv = bv(pop, simParam=SP))
#'
#' # Maternal breeding values
#' # (as if individuals are females)
#' (tmpMBv = bvM(pop, simParam=SP))
#'
#' # Paternal breeding values
#' # (as if individuals are males)
#' (tmpPBv = bvP(pop, simParam=SP))
#'
#' # Compare the values
#' data.frame(sex = pop@sex, bv = tmpBv[, 1], bvM = tmpMBv[, 1], bvP = tmpPBv[, 1])
#'
#' @export
bv = function(pop,simParam=NULL){
  ret = genParam(pop,simParam=simParam)$bv
  impTest = sapply(SP$traits, FUN = function(x) .hasSlot(x,"impEff"))
  if(any(impTest)){
    retM = bvM(pop,simParam=simParam)
    retP = bvP(pop,simParam=simParam)
    sel = pop@sex %in% "F"
    if(any(sel)){
      ret[sel, ] = retM[sel, ] # for females use mat breed val
    }
    sel = pop@sex %in% "M"
    if(any(sel)){
      ret[sel, ] = retP[sel, ] # for males use pat breed val
    }
  }
  return(ret)
}

#' @describeIn id Maternal breeding values (with imprinting)
#' @export
bvM = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$bvM
}

#' @describeIn id Paternal breeding values (with imprinting)
#' @export
bvP = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$bvP
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

#' @title Imprinting deviations
#'
#' @description Returns imprinting deviations for all traits according to the
#'   sex of individuals.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @details
#' If \code{sex == "H"} (in many plants) then mean imprinting deviation is returned,
#' which is 0 be definition.
#' If \code{sex == "F"} (female) then maternal imprinting deviation is returned.
#' If \code{sex == "M"} (male)   then paternal imprinting deviation is returned,
#' which is -(maternal imprinting deviation).
#'
#' If you must get maternal or paternal imprinting deviation irrespective of the
#' sex of individuals, use \code{idM()} and \code{idP()} - see examples.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' # --- Bisexual individuals ---
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAI(10, meanID=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@sex
#' # Mean imprinting deviations
#' # (individuals are bisexual)
#' (tmp = id(pop, simParam=SP))
#'
#' # Maternal imprinting deviations
#' # (as if individuals are females or for their female sexual organ)
#' (tmpM = idM(pop, simParam=SP))
#'
#' # Paternal imprinting deviations
#' # (as if individuals are males or for their male sexual organ)
#' (tmpP = idP(pop, simParam=SP))
#'
#' # Compare the values
#' data.frame(sex = pop@sex, id = tmp[, 1], idM = tmpM[, 1], idP = tmpP[, 1])
#'
#' # --- Male and female individuals ---
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$setSexes("yes_sys")
#' SP$addTraitAI(10, meanID=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' pop@sex
#' # Individual imprinting deviations
#' # (individuals are either females or males)
#' (tmp = id(pop, simParam=SP))
#'
#' # Maternal imprinting deviations
#' # (as if individuals are females)
#' (tmpM = idM(pop, simParam=SP))
#'
#' # Paternal imprinting deviations
#' # (as if individuals are males)
#' (tmpP = idP(pop, simParam=SP))
#'
#' # Compare the values
#' data.frame(sex = pop@sex, id = tmp[, 1], idM = tmpM[, 1], idP = tmpP[, 1])
#'
#' @export
id = function(pop,simParam=NULL){
  ret = idM(pop,simParam=simParam)
  sel = pop@sex %in% "H"
  if(any(sel)){
    ret[sel, ] = 0 # for bisexual use mean imp dev = 0
  }
  sel = pop@sex %in% "M" # for males use pat imp dev
  if(any(sel)){
    ret[sel, ] = -ret[sel, ] # pat imp dev = - mat imp dev
  }
  return(ret)
}

#' @describeIn id Maternal imprinting deviations
#' @export
idM = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$idM
}

#' @describeIn id Paternal imprinting deviations
#' @export
idP = function(pop,simParam=NULL){
  -genParam(pop,simParam=simParam)$idM
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
#' @description Returns additive genic variance (=sum of variances of breeding
#' values at individual loci) for all traits
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
#' @description Returns dominance genic variance (=sum of variances of dominance
#' deviations at individual loci) for all traits
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

#' @title Imprinting genic variance
#'
#' @description Returns imprinting genic variance (=sum of variances of
#' imprinting deviations at individual loci) for all traits
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
#' SP$addTraitAI(10, meanID=0.5)
#' SP$setVarE(h2=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarI(pop, simParam=SP)
#'
#' @export
genicVarI = function(pop,simParam=NULL){
  genParam(pop,simParam=simParam)$genicVarI
}

#' @title Additive-by-additive epistatic genic variance
#'
#' @description Returns additive-by-additive epistatic genic variance (=sum of
#' variances of additive-by-additive epistatic deviations at individual loci)
#' for all traits
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
#' @description Returns total genic variance (=sum of variances of genetic
#' values at individual loci) for all traits
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
