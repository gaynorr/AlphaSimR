# fmt: skip file

#' @title Mean genetic values
#'
#' @description Returns the mean genetic values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}} or \code{\link{HybridPop-class}}
#'
#' @details
#' The mean is calculated from the genetic values in \code{pop}.
#' It includes the trait intercept and can differ from both the intercept
#' and the requested founder mean.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
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
#' @details
#' This function calculated the mean of values stored in the \code{ebv} slot;
#' supplied by the user or \code{\link{setEBV}}.
#' Depending on how it was populated, it can contain estimated breeding
#' values, estimated genetic values, or other estimated/predicted values.
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
#' @details
#' This is the variance of genetic values in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
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
#' @details
#' This function calculated the variance of values stored in the \code{ebv} slot;
#' supplied by the user or \code{\link{setEBV}}.
#' Depending on how it was populated, it can contain estimated breeding
#' values, estimated genetic values, or other estimated/predicted values.
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

#' @title Calculate quantitative genetic quantities
#'
#' @description
#' Calculates quantitative genetic quantities and their variances
#' for an object of \code{\link{Pop-class}}
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details See \code{vignette("traits", package = "AlphaSimR")} and
#' \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory on trait and quantitative genetics implementation in AlphaSimR.
#' The text below is a very short summary of the theory.
#'
#' In quantitative genetics and AlphaSimR we use
#' different parameterizations of genetic values.
#' While these are equivalent
#' in the sense that they all parameterize the variation of genetic values,
#' the terms within each parameterization can have different meanings.
#'
#' To generate trait genetic values, AlphaSimR uses a trait model with
#' the genotypic parameterization of genetic values into
#' the trait intercept (\code{gv_mu}; fixed for a given trait definition) and
#' contributions from additive (\code{gv_a}), dominance (\code{gv_d}), and
#' additive-by-additive epistatic (\code{gv_aa}) genotypic effects.
#' For trait \code{t},
#' \code{gv[, t] = gv_mu[t] + gv_a[, t] + gv_d[, t] + gv_aa[, t] =
#' gv_mu[t] + gv_a[, t] + gv_n[, t]}.
#' Here, \code{gv_n = gv_d + gv_aa}. These contributions need not have
#' the population mean of zero. The intercept is calibrated to achieve the
#' requested founder mean and need not equal the population mean genetic value.
#' This genotypic parameterization is also called the "functional" or "biological" model,
#' though it is a simplified representation of the underlying trait biology.
#'
#' To generate classic quantitative genetic quantities, AlphaSimR uses
#' the breeding value parameterization of genetic values into
#' the population mean genetic value (\code{mu}; variable),
#' breeding value (\code{bv}), dominance deviation (\code{dd}), and
#' additive-by-additive epistatic deviation (\code{aa}).
#' For trait \code{t},
#' \code{gv[, t] = mu[t] + bv[, t] + dd[, t] + aa[, t] =
#' mu[t] + bv[, t] + nd[, t]}.
#' Here, \code{nd = dd + aa}.
#' The breeding values and all deviations have mean zero in \code{pop}.
#' The reference quantities are calculated from \code{pop}, so the same
#' individual can have different breeding values and deviations when
#' evaluated in different populations.
#' This breeding value parameterization is also called the "statistical" model
#' because it depends on statistical decomposition of genetic values.
#'
#' Note, the output \code{gv_mu} is the trait intercept \eqn{\mu_0},
#' whereas \code{mu} is the population mean of genetic value \eqn{\bar{G}}.
#' Despite its name, \code{gv_mu} is not the mean of \code{gv}.
#' For the returned quantities, the following identity holds
#' \code{mu = gv_mu + colMeans(gv_a + gv_d + gv_aa)}.
#' After scaling the genotypic effects,
#' the trait intercept is chosen as
#' the requested founder mean minus
#' the mean contribution of the genotypic effects in the founder population.
#' With the trait definition held fixed,
#' the intercept remains constant while the population mean can change.
#'
#' Each \code{covX_HW} is the difference between the corresponding genic
#' variance calculated using observed and Hardy-Weinberg genotype frequencies.
#' Each \code{covX_L}, for \code{X} in \code{G}, \code{A}, \code{D},
#' \code{AA}, and \code{N}, is calculated as \code{diag(varX) - genicVarX - covX_HW}.
#' These are within-trait variance contributions, which can be negative,
#' rather than covariances between traits.
#' The cross-component terms \code{covAD_L}, \code{covAAA_L}, \code{covDAA_L}, and
#' \code{covAN_L} are covariances.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than sample divisor \code{nInd - 1}.
#'
#' @return
#' \describe{
#' \item{varG}{an nTraits by nTraits matrix of total genetic variances and covariances between traits}
#' \item{varA}{an nTraits by nTraits matrix of additive genetic variances and covariances between traits}
#' \item{varD}{an nTraits by nTraits matrix of dominance genetic variances and covariances between traits}
#' \item{varAA}{an nTraits by nTraits matrix of additive-by-additive epistatic genetic variances and covariances between traits}
#' \item{varN}{an nTraits by nTraits matrix of non-additive genetic variances and covariances between traits}
#' \item{genicVarG}{an nTraits vector of total genic variances}
#' \item{genicVarA}{an nTraits vector of additive genic variances}
#' \item{genicVarD}{an nTraits vector of dominance genic variances}
#' \item{genicVarAA}{an nTraits vector of additive-by-additive epistatic genic variances}
#' \item{genicVarN}{an nTraits vector of non-additive genic variances}
#' \item{covG_HW}{an nTraits vector of adjustments to total genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covA_HW}{an nTraits vector of adjustments to additive genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covD_HW}{an nTraits vector of adjustments to dominance genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covAA_HW}{an nTraits vector of adjustments to additive-by-additive epistatic genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covN_HW}{an nTraits vector of adjustments to non-additive genic variances
#'   due to departures from Hardy-Weinberg genotype frequencies}
#' \item{covG_L}{an nTraits vector of contributions to total genetic variances
#'   due to linkage disequilibrium}
#' \item{covA_L}{an nTraits vector of contributions to additive genetic variances
#'   due to linkage disequilibrium}
#' \item{covD_L}{an nTraits vector of contributions to dominance genetic variances
#'   due to linkage disequilibrium}
#' \item{covAA_L}{an nTraits vector of contributions to additive-by-additive epistatic genetic variances
#'   due to linkage disequilibrium}
#' \item{covN_L}{an nTraits vector of contributions to non-additive genetic variances
#'   due to linkage disequilibrium}
#' \item{covAD_L}{an nTraits vector of within-trait covariances
#'   between breeding values and dominance deviations}
#' \item{covAAA_L}{an nTraits vector of within-trait covariances
#'   between breeding values and additive-by-additive epistatic deviations}
#' \item{covDAA_L}{an nTraits vector of within-trait covariances
#'   between dominance deviations and additive-by-additive epistatic deviations}
#' \item{covAN_L}{an nTraits vector of within-trait covariances
#'   between breeding values and non-additive deviations}
#' \item{gv}{a matrix of genetic values with dimensions nInd by nTraits}
#' \item{mu}{an nTraits vector of population mean of genetic values \eqn{\bar{G}} (breeding value parameterization)}
#' \item{mu_HW}{an nTraits vector of genetic value means expected under
#'   Hardy-Weinberg and linkage equilibrium at the observed allele frequencies}
#' \item{bv}{an nInd by nTraits matrix of breeding values (breeding value parameterization)}
#' \item{dd}{an nInd by nTraits matrix of dominance deviations (breeding value parameterization)}
#' \item{aa}{an nInd by nTraits matrix of additive-by-additive epistatic deviations (breeding value parameterization)}
#' \item{nd}{an nInd by nTraits matrix of non-additive deviations (breeding value parameterization)}
#' \item{gv_mu}{an nTraits vector of trait intercepts \eqn{\mu_0} (genotypic parameterization)}
#' \item{gv_a}{an nInd by nTraits matrix of additive genotypic effect
#'   contributions to genetic values (genotypic parameterization)}
#' \item{gv_d}{an nInd by nTraits matrix of dominance genotypic effect
#'   contributions to genetic values (genotypic parameterization)}
#' \item{gv_aa}{an nInd by nTraits matrix of additive-by-additive epistatic genotypic effect
#'   contributions to genetic values (genotypic parameterization)}
#' \item{gv_n}{an nInd by nTraits matrix of non-additive genotypic effects
#'   contributions to genetic values (genotypic parameterization)}
#' \item{alpha}{a list of average effects of allele substitution
#'   calculated using observed marginal genotype frequencies, with length nTraits}
#' \item{alpha_HW}{a list of average effects of allele substitution
#'   calculated using Hardy-Weinberg genotype frequencies at the observed
#'   allele frequencies, with length nTraits}
#' }
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genParam(pop, simParam=SP)
#'
#' @export
genParam = function(pop,simParam=NULL,nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }

  nInd = nInd(pop)
  nTraits = simParam$nTraits
  traitNames = simParam$traitNames

  # Blank nInd x nTrait matrices
  gv = matrix(NA_real_, nrow=nInd, ncol=nTraits)
  colnames(gv) = traitNames
  bv = dd = aa = nd = gv_a = gv_d = gv_aa = gv_n = gv

  # Blank nTrait vectors
  genicVarA = rep(NA_real_, nTraits)
  names(genicVarA) = traitNames
  genicVarD = genicVarAA = genicVarN =
    covG_HW = covA_HW = covD_HW = covAA_HW = covN_HW =
    covAD_L = covAAA_L = covDAA_L = covAN_L =
    mu = mu_HW = gv_mu = genicVarA

  # Average effect of an allele substitution
  alpha = vector("list", length=nTraits)
  names(alpha) = traitNames
  alpha_HW = alpha

  #Loop through trait calculations
  for(i in seq_len(nTraits)){
    trait = simParam$traits[[i]]
    tmp = calcGenParam(trait,pop,nThreads)
    genicVarA[i] = tmp$genicVarA2
    genicVarN[i] = 0
    covA_HW[i] = tmp$genicVarA-tmp$genicVarA2
    covN_HW[i] = 0
    gv[,i] = tmp$gv
    mu[i] = tmp$mu # mean gv for the trait in the pop
    mu_HW[i] = tmp$mu_HWE
    bv[,i] = tmp$bv
    nd[,i] = rep(0,pop@nInd)
    gv_mu[i] = tmp$gv_mu # trait gv intercept
    gv_a[,i] = tmp$gv_a
    gv_n[,i] = rep(0,pop@nInd)
    if(.hasSlot(trait,"domEff")){
      genicVarD[i] = tmp$genicVarD2
      genicVarN[i] = genicVarN[i] + genicVarD[i]
      covD_HW[i] = tmp$genicVarD-tmp$genicVarD2
      covN_HW[i] = covN_HW[i] + covD_HW[i]
      dd[,i] = tmp$dd
      nd[,i] = nd[,i] + dd[,i]
      gv_d[,i] = tmp$gv_d
      gv_n[,i] = gv_n[,i] + gv_d[,i]
    }else{
      genicVarD[i] = 0
      covD_HW[i] = 0
      dd[,i] = rep(0,pop@nInd)
      gv_d[,i] = rep(0,pop@nInd)
    }
    if(.hasSlot(trait,"epiEff")){
      genicVarAA[i] = tmp$genicVarAA2
      genicVarN[i] = genicVarN[i] + genicVarAA[i]
      covAA_HW[i] = tmp$genicVarAA-tmp$genicVarAA2
      covN_HW[i] = covN_HW[i] + covAA_HW[i]
      aa[,i] = tmp$aa
      nd[,i] = nd[,i] + aa[,i]
      gv_aa[,i] = tmp$gv_aa
      gv_n[,i] = gv_n[,i] + gv_aa[,i]
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
      covAN_L[i] = 0
    } else {
      covAD_L[i] = popVar(cbind(bv[,i],dd[,i]))[1,2]
      covAAA_L[i] = popVar(cbind(bv[,i],aa[,i]))[1,2]
      covDAA_L[i] = popVar(cbind(dd[,i],aa[,i]))[1,2]
      covAN_L[i] = popVar(cbind(bv[,i],nd[,i]))[1,2]
    }
    alpha[[i]] = tmp$alpha
    alpha_HW[[i]] = tmp$alpha_HW
  }

  varG = popVar(gv)
  rownames(varG) = colnames(varG) = traitNames

  varA = popVar(bv)
  rownames(varA) = colnames(varA) = traitNames

  varD = popVar(dd)
  rownames(varD) = colnames(varD) = traitNames

  varAA = popVar(aa)
  rownames(varAA) = colnames(varAA) = traitNames

  varN = popVar(nd)
  rownames(varN) = colnames(varN) = traitNames

  genicVarG = genicVarA + genicVarD + genicVarAA
  covG_HW = covA_HW + covD_HW + covAA_HW

  output = list(varG=varG,
                varA=varA,
                varD=varD,
                varAA=varAA,
                varN=varN,
                genicVarG=genicVarG,
                genicVarA=genicVarA,
                genicVarD=genicVarD,
                genicVarAA=genicVarAA,
                genicVarN=genicVarN,
                covG_HW=covG_HW,
                covA_HW=covA_HW,
                covD_HW=covD_HW,
                covAA_HW=covAA_HW,
                covN_HW=covN_HW,
                covG_L=diag(varG)-genicVarG-covG_HW,
                covA_L=diag(varA)-genicVarA-covA_HW,
                covD_L=diag(varD)-genicVarD-covD_HW,
                covAA_L=diag(varAA)-genicVarAA-covAA_HW,
                covN_L=diag(varN)-genicVarN-covN_HW,
                covAD_L=covAD_L,
                covAAA_L=covAAA_L,
                covDAA_L=covDAA_L,
                covAN_L=covAN_L,
                gv=gv,
                mu=mu, # mean gv for the trait in the pop
                mu_HW=mu_HW,
                bv=bv,
                dd=dd,
                aa=aa,
                nd=nd,
                gv_mu=gv_mu, # trait gv intercept
                gv_a=gv_a,
                gv_d=gv_d,
                gv_aa=gv_aa,
                gv_n=gv_n,
                alpha=alpha,
                alpha_HW=alpha_HW)
  return(output)
}

#' @title Additive genetic variance
#'
#' @description Returns additive genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' This is the variance of additive genetic (breeding) values in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varA(pop, simParam=SP)
#'
#' @export
varA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varA
}

#' @title Dominance genetic variance
#'
#' @description Returns dominance genetic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' This is the variance of dominance deviations in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varD(pop, simParam=SP)
#'
#' @export
varD = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varD
}

#' @title Additive-by-additive epistatic genetic variance
#'
#' @description Returns additive-by-additive epistatic genetic
#' variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' This is the variance of additive-by-additive epistatic deviations in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varAA(pop, simParam=SP)
#'
#' @export
varAA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varAA
}

#' @title Non-additive genetic variance
#'
#' @description Returns non-additive genetic variance for all traits
#'   (includes dominance and epistatic genetic variance)
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' This is the variance of total non-additive deviations in \code{pop}.
#' Variances and covariances use the population divisor \code{nInd},
#' rather than the sample divisor \code{nInd - 1} (see \code{\link{popVar}}).
#' Since \code{nd = dd + aa}, the variance includes twice the covariance
#' between dominance and additive-by-additive epistatic deviations.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nTraits by nTraits matrix of variances and covariances between traits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop)
#' varA(pop, simParam=SP)
#' varN(pop, simParam=SP)
#'
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' varG(pop)
#' varA(pop, simParam=SP)
#' varD(pop, simParam=SP)
#' varAA(pop, simParam=SP)
#' varN(pop, simParam=SP)
#'
#' @export
varN = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$varN
}

#' @title Breeding value
#'
#' @description Returns breeding values for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Breeding values derive from the breeding value parameterization of genetic values.
#' They have mean zero in \code{pop},
#' which supplies the reference quantities used by \code{\link{genParam}}.
#' Evaluating the same individual in a different population can therefore
#' change its returned values.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nInd by nTraits matrix.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' bv(pop, simParam=SP)
#'
#' @export
bv = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$bv
}

#' @title Dominance deviations
#'
#' @description Returns dominance deviations for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Dominance deviations derive from the breeding value parameterization of genetic values.
#' They have mean zero in \code{pop},
#' which supplies the reference quantities used by \code{\link{genParam}}.
#' Evaluating the same individual in a different population can therefore
#' change returned values.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nInd by nTraits matrix.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' dd(pop, simParam=SP)
#'
#' @export
dd = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$dd
}

#' @title Additive-by-additive epistatic deviations
#'
#' @description Returns additive-by-additive epistatic
#' deviations for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Additive-by-additive epistatic deviations derive from the breeding value parameterization of genetic values.
#' They have mean zero in \code{pop},
#' which supplies the reference quantities used by \code{\link{genParam}}.
#' Evaluating the same individual in a different population can therefore
#' change its returned values.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nInd by nTraits matrix.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' aa(pop, simParam=SP)
#'
#' @export
aa = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$aa
}

#' @title Non-additive deviations
#'
#' @description Returns non-additive deviations for all traits
#'   (includes dominance and epistatic deviations)
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Non-additive deviations derive from the breeding value parameterization of genetic values.
#' They have mean zero in \code{pop},
#' which supplies the reference quantities used by \code{\link{genParam}}.
#' Evaluating the same individual in a different population can therefore
#' change its returned values.
#' The total non-additive deviation is \code{nd = dd + aa}.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return An nInd by nTraits matrix.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' cbind(dd(pop, simParam=SP),
#'       aa(pop, simParam=SP),
#'       nd(pop, simParam=SP))
#'
#' @export
nd = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$nd
}

#' @title Total genic variance
#'
#' @description Returns total genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#'
#' @export
genicVarG = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarG
}

#' @title Additive genic variance
#'
#' @description Returns additive genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#'
#' @export
genicVarA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarA
}

#' @title Dominance genic variance
#'
#' @description Returns dominance genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarD(pop, simParam=SP)
#'
#' @export
genicVarD = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarD
}

#' @title Additive-by-additive epistatic genic variance
#'
#' @description Returns additive-by-additive epistatic
#' genic variance for all traits
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarAA(pop, simParam=SP)
#'
#' @export
genicVarAA = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarAA
}

#' @title Non-additive genic variance
#'
#' @description Returns non-additive genic variance for all traits
#'   (includes dominance and epistatic genic variance)
#' @param pop an object of \code{\link{Pop-class}}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @details
#' Genic variances refer to an idealized population with the same allele
#' frequencies as \code{pop}, but with Hardy-Weinberg genotype frequencies
#' and linkage equilibrium.
#' They are theoretical expectations under this reference distribution,
#' rather than variances of values evaluated in the observed population.
#' See \code{\link{genParam}} for the corresponding genetic variances
#' and contributions from departures from equilibrium.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' @return A vector of length nTraits.
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#' genicVarN(pop, simParam=SP)
#'
#' SP = SimParam$new(founderPop)
#' SP$addTraitADE(10, meanDD=0.5, relAA=0.2)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' genicVarG(pop, simParam=SP)
#' genicVarA(pop, simParam=SP)
#' genicVarD(pop, simParam=SP)
#' genicVarAA(pop, simParam=SP)
#' genicVarN(pop, simParam=SP)
#'
#' @export
genicVarN = function(pop,simParam=NULL,nThreads=NULL){
  genParam(pop,simParam=simParam,nThreads=nThreads)$genicVarN
}

#' @title Genetic value
#'
#' @description A wrapper for accessing the gv slot
#'
#' @param pop a \code{\link{Pop-class}} or similar object
#'
#' @details
#' Genetic values include the trait intercept and the contributions of
#' additive, dominance, and additive-by-additive epistatic genotypic effects
#' present in the trait model.
#' For a fixed trait definition, these values do not depend on the population,
#' unlike the breeding values and non-additive deviations.
#' For traits with genotype-by-environment effects, these genetic values
#' refer to the target environment, where the environmental covariate is zero.
#'
#' See \code{vignette("traits", package = "AlphaSimR")} for the trait model.
#'
#' See \code{vignette("QuanGen", package = "AlphaSimR")} for
#' background theory and its demonstration.
#'
#' See \code{vignette("GxE", package = "AlphaSimR")} for
#' the genotype-by-environment model
#'
#' @return An nInd by nTraits matrix of genetic values.
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
#' @details
#' The \code{ebv} slot stores predictions supplied by the user or
#' \code{\link{setEBV}}.
#' Depending on how it was populated, it can contain estimated breeding
#' values, estimated genetic values, or other estimated/predicted values.
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
#' ebv(pop)
#'
#' @export
ebv = function(pop){
  pop@ebv
}

#' @title Calculate parent average
#'
#' @param pop \code{\link{Pop-class}} with individuals whose parent average
#'   will be calculated
#' @param parents \code{\link{Pop-class}} with mothers and fathers of individuals
#'   in \code{pop}; if \code{NULL} must provide \code{mothers} and \code{fathers}
#' @param mothers \code{\link{Pop-class}} with mothers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param fathers \code{\link{Pop-class}} with fathers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param use character, calculate using \code{"\link{gv}"},
#'   \code{"\link{ebv}"}, or \code{"\link{pheno}"}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return a matrix of parent averages with dimensions nInd by nTraits
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
#' pop2 = randCross(pop, nCrosses=10, nProgeny=2)
#' parentAverage(pop2, parents = pop)
#' parentAverage(pop2, mothers = pop, fathers = pop)
#'
#' @export
parentAverage = function(pop, parents = NULL, mothers = NULL, fathers = NULL,
                         use = "gv", simParam = NULL, nThreads=NULL) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if (!is.null(parents)) {
    matchMothers = match(x = pop@mother, table = parents@id)
    matchFathers = match(x = pop@father, table = parents@id)
  } else {
    if (is.null(mothers) | is.null(fathers)) {
      stop("must provide either 'parents' or both 'mothers' and 'fathers'!")
    }
    matchMothers = match(x = pop@mother, table = mothers@id)
    matchFathers = match(x = pop@father, table = fathers@id)
  }
  if (anyNA(matchMothers)) {
    stop("some parents/mothers not found!")
  }
  if (anyNA(matchFathers)) {
    stop("some parents/fathers not found!")
  }
  if (use %in% c("gv", "ebv", "pheno")) {
    if (!is.null(parents)) {
      ret = 0.5 * (slot(object = parents, name = use)[matchMothers, , drop = FALSE] +
                   slot(object = parents, name = use)[matchFathers, , drop = FALSE])
    } else {
      ret = 0.5 * (slot(object = mothers, name = use)[matchMothers, , drop = FALSE] +
                   slot(object = fathers, name = use)[matchFathers, , drop = FALSE])
    }
  # Commented out because bv() uses different reference populations for parents and offspring.
  # TODO: See https://github.com/gaynorr/AlphaSimR/issues/291
  # } else if (use == "bv") {
  #   if (!is.null(parents)) {
  #     ret = 0.5 * (bv(parents, simParam = simParam,
  #                     nThreads=nThreads)[matchMothers, , drop = FALSE] +
  #                  bv(parents, simParam = simParam,
  #                     nThreads=nThreads)[matchFathers, , drop = FALSE])
  #   } else {
  #     ret = 0.5 * (bv(mothers, simParam = simParam,
  #                     nThreads=nThreads)[matchMothers, , drop = FALSE] +
  #                  bv(fathers, simParam = simParam,
  #                     nThreads=nThreads)[matchFathers, , drop = FALSE])
  #   }
  } else {
    stop("use must be one of 'gv', 'ebv', or 'pheno'!")
  }
  return(ret)
}

#' @title Calculate Mendelian sampling
#'
#' @param pop \code{\link{Pop-class}} with individuals whose Mendelian samplings
#'   will be calculated
#' @param parents \code{\link{Pop-class}} with mothers and fathers of individuals
#'   in \code{pop}; if \code{NULL} must provide \code{mothers} and \code{fathers}
#' @param mothers \code{\link{Pop-class}} with mothers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param fathers \code{\link{Pop-class}} with fathers of individuals in \code{pop};
#'   if \code{NULL} must provide \code{parents}
#' @param use character, calculate using \code{"\link{gv}"},
#'   \code{"\link{ebv}"}, or \code{"\link{pheno}"}
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return a matrix of Mendelian samplings with dimensions nInd by nTraits
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
#' pop2 = randCross(pop, nCrosses=10, nProgeny=2)
#' mendelianSampling(pop2, parents = pop)
#' mendelianSampling(pop2, mothers = pop, fathers = pop)
#'
#' @export
mendelianSampling = function(pop, parents = NULL, mothers = NULL, fathers = NULL,
                             use = "gv", simParam = NULL, nThreads=NULL) {
  if (is.null(simParam)) {
    simParam = get("SP", envir = .GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  pa = parentAverage(pop = pop, parents = parents, mothers = mothers, fathers = fathers,
                     use = use, simParam = simParam, nThreads=nThreads)
  if (use %in% c("gv", "ebv", "pheno")) {
    ret = slot(object = pop, name = use) - pa
  # Commented out because bv() uses different reference populations for parents and offspring.
  # TODO: See https://github.com/gaynorr/AlphaSimR/issues/291
  # } else if (use == "bv") {
  #   ret = bv(pop, simParam = simParam, nThreads=nThreads) - pa
  } else {
    stop("use must be one of 'gv', 'ebv', or 'pheno'!")
  }
  return(ret)
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
