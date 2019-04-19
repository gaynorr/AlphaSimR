#' @title Fast RR-BLUP
#'
#' @description
#' Solves an RR-BLUP model for genomic predictions given known variance 
#' components. This implementation is meant as a fast and low memory 
#' alternative to \code{\link{RRBLUP}} or \code{\link{RRBLUP2}}. Unlike 
#' the those functions, the fastRRBLUP does not fit fixed effects (other 
#' than the intercept) or account for unequal replication. 
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model or a
#' function of the traits returning a single value. Only univariate models 
#' are supported.
#' @param use train model using phenotypes "pheno", genetic values "gv", 
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. 
#' @param Vu marker effect variance. If value is NULL, a 
#' reasonable value is chosen automatically.
#' @param Ve error variance. If value is NULL, a 
#' reasonable value is chosen automatically.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = fastRRBLUP(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
fastRRBLUP = function(pop, traits=1, use="pheno", snpChip=1, 
                      useQtl=FALSE, maxIter=1000, Vu=NULL, Ve=NULL, 
                      simParam=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,...)
  #fixEff = as.integer(factor(pop@fixEff))
  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }
  # Sort out Vu and Ve
  if(is.function(traits)){
    if(is.null(Vu)){
      Vu = var(y)/nLoci
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot(length(traits)==1)
    if(is.null(Vu)){
      Vu = 2*simParam$varA[traits]/nLoci
      if(is.na(Vu)){
        Vu = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traits]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }
  #Fit model
  ans = callFastRRBLUP(y,pop@geno,lociPerChr,
                       lociLoc,Vu,Ve,maxIter,
                       simParam$nThreads)
  output = new("RRsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               markerEff=ans$u,
               fixEff=as.matrix(ans$beta),
               Vu=as.matrix(Vu),
               Ve=as.matrix(Ve),
               iter=0)
  return(output)
}



#' @title RR-BLUP Model
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait or traits to model, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv", 
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. Only used 
#' when number of traits is greater than 1.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = RRBLUP(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
RRBLUP = function(pop, traits=1, use="pheno", snpChip=1, 
                  useQtl=FALSE, maxIter=1000L, simParam=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,...)
  fixEff = as.integer(factor(pop@fixEff))
  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }
  #Fit model
  if(ncol(y)>1){
    ans = callRRBLUP_MV(y,fixEff,pop@reps,pop@geno,lociPerChr,
                        lociLoc, maxIter, simParam$nThreads)
  }else{
    ans = callRRBLUP(y,fixEff,pop@reps,pop@geno,lociPerChr,lociLoc,
                     simParam$nThreads)
  }
  markerEff=ans$u
  if(is.null(ans[["iter"]])){
    iter = 0
  }else{
    iter = ans$iter
  }
  output = new("RRsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               markerEff=markerEff,
               fixEff=ans$beta[1,,drop=FALSE],
               Vu=as.matrix(ans$Vu),
               Ve=as.matrix(ans$Ve),
               iter=iter)
  return(output)
}

#' @title RR-BLUP Model 2
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions. This implementation is 
#' meant for situations where \code{\link{RRBLUP}} is too slow. Note that 
#' RRBLUP2 is only faster in certain situations, see details below. Most 
#' users should use \code{\link{RRBLUP}}.
#' 
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model or a
#' function of the traits returning a single value. Unlike \code{\link{RRBLUP}}, 
#' only univariate models are supported.
#' @param use train model using phenotypes "pheno", genetic values "gv", 
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. 
#' @param Vu marker effect variance. If value is NULL, a 
#' reasonable starting point is chosen automatically.
#' @param Ve error variance. If value is NULL, a 
#' reasonable starting point is chosen automatically.
#' @param useEM use EM to solve variance components. If false, 
#' the initial values are considered true.
#' @param tol tolerance for EM algorithm convergence
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#' 
#' @details 
#' The RRBLUP2 function works best when the number of markers is not 
#' too large. This is because it solves the RR-BLUP problem by setting 
#' up and solving Henderson's mixed model equations. Solving these equations 
#' involves a square matrix with dimensions equal to the number of fixed 
#' effects plus the number of random effects (markers). Whereas the \code{\link{RRBLUP}} 
#' function solves the RR-BLUP problem using the EMMA approach. This approach involves 
#' a square matrix with dimensions equal to the number of phenotypic records. This means 
#' that the RRBLUP2 function uses less memory than RRBLUP when the number of markers 
#' is approximately equal to or smaller than the number of phenotypic records. 
#' 
#' The RRBLUP2 function is not recommend for cases where the variance components are 
#' unknown. This is uses the EM algorithm to solve for unknown variance components, 
#' which is generally considerably slower than the EMMA approach of \code{\link{RRBLUP}}. 
#' The number of iterations for the EM algorith is set by maxIter. The default value 
#' is typically too small for convergence. When the algorithm fails to converage a 
#' warning is displayed, but results are given for the last iteration. These results may 
#' be "good enough". However we make no claim to this effect, because we can not generalize 
#' to all possible use cases.
#' 
#' The RRBLUP2 function can quickly solve the mixed model equations without estimating variance 
#' components. The variance components are set by defining Vu and Ve. Estimation of components 
#' is suppressed by setting useEM to false. This may be useful if the model is being retrained 
#' multiple times during the simulation. You could run \code{\link{RRBLUP}} function the first 
#' time the model is trained, and then use the variance components from this output for all 
#' future runs with the RRBLUP2 functions. Again, we can make no claim to the general robustness 
#' of this approach.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = RRBLUP2(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
RRBLUP2 = function(pop, traits=1, use="pheno", snpChip=1, 
                   useQtl=FALSE, maxIter=10, Vu=NULL, Ve=NULL, 
                   useEM=TRUE, tol=1e-6, simParam=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,...)
  fixEff = as.integer(factor(pop@fixEff))
  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }
  # Sort out Vu and Ve
  if(is.function(traits)){
    if(is.null(Vu)){
      Vu = var(y)/nLoci
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot(length(traits)==1)
    if(is.null(Vu)){
      Vu = 2*simParam$varA[traits]/nLoci
      if(is.na(Vu)){
        Vu = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traits]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }
  #Fit model
  ans = callRRBLUP2(y,fixEff,pop@reps,pop@geno,lociPerChr,
                    lociLoc,Vu,Ve,tol,maxIter,useEM,simParam$nThreads)
  output = new("RRsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               markerEff=ans$u,
               fixEff=ans$beta[1,,drop=FALSE],
               Vu=as.matrix(ans$Vu),
               Ve=as.matrix(ans$Ve),
               iter=ans$iter)
  return(output)
}

#' @title RR-BLUP Model with Dominance
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions that includes 
#' dominance effects.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv", 
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. Only used 
#' when number of traits is greater than 1.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = RRBLUP_D(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
RRBLUP_D = function(pop, traits=1, use="pheno", snpChip=1, 
                    useQtl=FALSE, maxIter=40L, simParam=NULL, 
                    ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,...)
  fixEff = as.integer(factor(pop@fixEff))
  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }
  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_D(y,fixEff,pop@reps,pop@geno,lociPerChr,
                     lociLoc,maxIter,simParam$nThreads)
  output = new("RRDsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               markerEff=ans$alpha,
               addEff=ans$ans$u[[1]],
               domEff=ans$d,
               fixEff=ans$ans$fixEff[1,,drop=FALSE],
               Vu=ans$ans$Vu,
               Ve=ans$ans$Ve,
               iter=ans$ans$iter)
  return(output)
}


#' @title RR-BLUP with Dominance Model 2
#'
#' @description
#' Fits an RR-BLUP model for genomic predictions that includes 
#' dominance effects. This implementation is meant for situations where 
#' \code{\link{RRBLUP_D}} is too slow. Note that RRBLUP_D2 
#' is only faster in certain situations. Most users should use 
#' \code{\link{RRBLUP_D}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv", 
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations. Only used 
#' when number of traits is greater than 1.
#' @param Va marker effect variance for additive effects. If value is NULL, 
#' a reasonable starting point is chosen automatically.
#' @param Vd marker effect variance for dominance effects. If value is NULL, 
#' a reasonable starting point is chosen automatically.
#' @param Ve error variance. If value is NULL, a 
#' reasonable starting point is chosen automatically.
#' @param useEM use EM to solve variance components. If false, 
#' the initial values are considered true.
#' @param tol tolerance for EM algorithm convergence
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitAD(10, meanDD=0.5)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = RRBLUP_D2(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
RRBLUP_D2 = function(pop, traits=1, use="pheno", snpChip=1, 
                     useQtl=FALSE, maxIter=10, Va=NULL, Vd=NULL, 
                     Ve=NULL, useEM=TRUE, tol=1e-6, simParam=NULL, 
                     ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,...)
  fixEff = as.integer(factor(pop@fixEff))
  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }
  # Sort out Va, Vd and Ve
  if(is.function(traits)){
    if(is.null(Va)){
      Va = var(y)/nLoci
    }
    if(is.null(Vd)){
      Vd = var(y)/nLoci
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot(length(traits)==1)
    if(is.null(Va)){
      Va = 2*simParam$varA[traits]/nLoci
      if(is.na(Va)){
        Va = var(y)/nLoci
      }
    }
    if(is.null(Vd)){
      Vd = 2*simParam$varA[traits]/nLoci
      if(is.na(Vd)){
        Vd = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traits]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }
  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_D2(y,fixEff,pop@reps,pop@geno,lociPerChr,
                      lociLoc,maxIter,Va,Vd,Ve,tol,useEM,
                      simParam$nThreads)
  output = new("RRDsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               markerEff=ans$alpha,
               addEff=ans$ans$u[,1,drop=FALSE],
               domEff=ans$d,
               fixEff=ans$ans$fixEff[1,,drop=FALSE],
               Vu=ans$ans$Vu,
               Ve=as.matrix(ans$ans$Ve),
               iter=ans$ans$iter)
  return(output)
}

#' @title RR-BLUP GCA Model
#'
#' @description
#' Fits an RR-BLUP model that estimates seperate marker effects for
#' females and males. Useful for predicting GCA of parents
#' in single cross hybrids. Can also predict performance of specific 
#' single cross hybrids.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv", 
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations for convergence.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = RRBLUP_GCA(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
RRBLUP_GCA = function(pop, traits=1, use="pheno", snpChip=1, 
                      useQtl=FALSE, maxIter=40L, simParam=NULL, 
                      ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,...)
  fixEff = as.integer(factor(pop@fixEff))
  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }
  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_GCA(y,fixEff,pop@reps,pop@geno,
                       lociPerChr,lociLoc,maxIter,
                       simParam$nThreads)
  output = new("GCAsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               femaleEff=ans$u[[1]],
               maleEff=ans$u[[2]],
               fixEff=ans$beta,
               Vu=ans$Vu,
               Ve=ans$Ve,
               iter=ans$iter)
  return(output)
}

#' @title RR-BLUP GCA Model 2
#'
#' @description
#' Fits an RR-BLUP model that estimates seperate marker effects for
#' females and males. This implementation is meant for situations where 
#' \code{\link{RRBLUP_GCA}} is too slow. Note that RRBLUP_GCA2 
#' is only faster in certain situations. Most users should use 
#' \code{\link{RRBLUP_GCA}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv", 
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations for convergence.
#' @param VuF marker effect variance for females. If value is NULL, a 
#' reasonable starting point is chosen automatically.
#' @param VuM marker effect variance for males. If value is NULL, a 
#' reasonable starting point is chosen automatically.
#' @param Ve error variance. If value is NULL, a 
#' reasonable starting point is chosen automatically.
#' @param useEM use EM to solve variance components. If false, 
#' the initial values are considered true.
#' @param tol tolerance for EM algorithm convergence
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = RRBLUP_GCA2(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
RRBLUP_GCA2 = function(pop, traits=1, use="pheno", snpChip=1, 
                       useQtl=FALSE, maxIter=10, VuF=NULL, VuM=NULL, 
                       Ve=NULL, useEM=TRUE, tol=1e-6, simParam=NULL, 
                       ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,...)
  fixEff = as.integer(factor(pop@fixEff))
  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }
  # Sort out VuF, VuM and Ve
  if(is.function(traits)){
    if(is.null(VuF)){
      VuF = var(y)/nLoci
    }
    if(is.null(VuM)){
      VuM = var(y)/nLoci
    }
    if(is.null(Ve)){
      Ve = var(y)/2
    }
  }else{
    stopifnot(length(traits)==1)
    if(is.null(VuF)){
      VuF = 2*simParam$varA[traits]/nLoci
      if(is.na(VuF)){
        VuF = var(y)/nLoci
      }
    }
    if(is.null(VuM)){
      VuM = 2*simParam$varA[traits]/nLoci
      if(is.na(VuM)){
        VuM = var(y)/nLoci
      }
    }
    if(is.null(Ve)){
      Ve = simParam$varE[traits]
      if(is.na(Ve)){
        Ve = var(y)/2
      }
    }
  }
  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_GCA2(y,fixEff,pop@reps,pop@geno,
                        lociPerChr,lociLoc,maxIter,
                        VuF,VuM,Ve,tol,useEM,
                        simParam$nThreads)
  output = new("GCAsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               femaleEff=ans$u[,1,drop=FALSE],
               maleEff=ans$u[,2,drop=FALSE],
               fixEff=ans$beta,
               Vu=ans$Vu,
               Ve=as.matrix(ans$Ve),
               iter=ans$iter)
  return(output)
}

#' @title RR-BLUP SCA Model
#'
#' @description
#' An extention of \code{\link{RRBLUP_GCA}} that adds dominance effects. 
#' Note that we have not seen any consistent benefit of this model over 
#' \code{\link{RRBLUP_GCA}}.
#'
#' @param pop a \code{\link{Pop-class}} to serve as the training population
#' @param traits an integer indicating the trait to model, or a
#' function of the traits returning a single value.
#' @param use train model using phenotypes "pheno", genetic values "gv", 
#' estimated breeding values "ebv", breeding values "bv", or randomly "rand"
#' @param snpChip an integer indicating which SNP chip genotype 
#' to use
#' @param useQtl should QTL genotypes be used instead of a SNP chip. 
#' If TRUE, snpChip specifies which trait's QTL to use, and thus these 
#' QTL may not match the QTL underlying the phenotype supplied in traits.
#' @param maxIter maximum number of iterations for convergence.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = RRBLUP_SCA(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
RRBLUP_SCA = function(pop, traits=1, use="pheno", snpChip=1, 
                      useQtl=FALSE, maxIter=40L, simParam=NULL, ...){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  y = getResponse(pop=pop,trait=traits,use=use,
                  simParam=simParam,...)
  fixEff = as.integer(factor(pop@fixEff))
  if(useQtl){
    nLoci = simParam$traits[[snpChip]]@nLoci
    lociPerChr = simParam$traits[[snpChip]]@lociPerChr
    lociLoc = simParam$traits[[snpChip]]@lociLoc
  }else{
    nLoci = simParam$snpChips[[snpChip]]@nLoci
    lociPerChr = simParam$snpChips[[snpChip]]@lociPerChr
    lociLoc = simParam$snpChips[[snpChip]]@lociLoc
  }
  #Fit model
  stopifnot(ncol(y)==1)
  ans = callRRBLUP_SCA(y,fixEff,pop@reps,pop@geno,
                       lociPerChr,lociLoc,maxIter,
                       simParam$nThreads)
  ans = ans$ans
  output = new("SCAsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               femaleEff=ans$u[[1]],
               maleEff=ans$u[[2]],
               d=ans$u[[3]],
               fixEff=ans$beta,
               Vu=ans$Vu,
               Ve=ans$Ve,
               iter=ans$iter)
  return(output)
}

#' @title Set EBV
#'
#' @description
#' Sets a population's EBV with genomic estimated
#' values from \code{\link{RRBLUP}}, \code{\link{RRBLUP_GCA}},
#' or \code{\link{RRBLUP_SCA}}.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param solution an object of \code{\link{RRsol-class}},
#' \code{\link{SCAsol-class}}, or \code{\link{GCAsol-class}}
#' @param gender either NULL, "male" or "female". If 
#' solution is \code{\link{GCAsol-class}} or 
#' \code{\link{SCAsol-class}} the EBV is the GCA if used in 
#' the corresponding pool
#' @param useGV if model is \code{\link{RRDsol-class}}, 
#' setting this parameter to TRUE will give use estimated 
#' genetic values. Otherwise, you get estimated breeding 
#' values that depend on the population's allele frequency.
#' @param append should EBVs be appended to existing EBVs
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#' SP$addSnpChip(10)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Run GS model and set EBV
#' ans = RRBLUP(pop, simParam=SP)
#' pop = setEBV(pop, ans, simParam=SP)
#' 
#' #Evaluate accuracy
#' cor(gv(pop), ebv(pop))
#' 
#' @export
setEBV = function(pop, solution, gender=NULL, useGV=FALSE, 
                  append=FALSE, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(class(solution)=="RRsol"){
    ebv = gebvRR(solution, pop, simParam$nThreads)
  }else if(class(solution)=="RRDsol"){
    if(useGV){
      ebv = gegvRRD(solution, pop, simParam$nThreads)
    }else{
      ebv = gebvRR(solution, pop, simParam$nThreads)
    }
  }else{
    if(is.null(gender)){
      if(class(solution)=="GCAsol"){
        ebv = gegvGCA(solution, pop, simParam$nThreads)
      }else{
        ebv = gegvSCA(solution, pop, simParam$nThreads)
      }
    }else if(toupper(gender)=="FEMALE"){
      ebv = gebvGCA(solution, pop, TRUE, simParam$nThreads)
    }else if(toupper(gender)=="MALE"){
      ebv = gebvGCA(solution, pop, FALSE, simParam$nThreads)
    }else{
      stop(paste0("gender=",gender," is not a valid option"))
    }
  }
  if(append){
    pop@ebv = cbind(pop@ebv,ebv)
  }else{
    pop@ebv = ebv
  }
  return(pop)
}

#' @title RRBLUP Memory Usage
#'
#' @description
#' Estimates the amount of RAM needed to run the \code{\link{RRBLUP}}
#' and its related functions for a given training population size. 
#' Note that this functions may underestimate total usage.
#'
#' @param nInd the number of individuals in the training population
#' @param nMarker the number of markers per individual
#' @param model either "REG", "GCA", or "SCA" for \code{\link{RRBLUP}} 
#' \code{\link{RRBLUP_GCA}} and \code{\link{RRBLUP_SCA}} respectively.
#'
#' @return Returns an estimate for the required gigabytes of RAM
#'
#' @examples 
#' RRBLUPMemUse(nInd=1000, nMarker=5000)
#' 
#' @export
RRBLUPMemUse = function(nInd,nMarker,model="REG"){
  y = nInd
  X = nInd #times fixed effects, assuming 1 here
  M = nInd*nMarker
  u = nMarker
  if(toupper(model)=="REG"){
    S = nInd*nInd
    eigval = nInd
    eigvec = nInd*nInd
    eta = nInd
    Hinv = nInd*nInd
  }else if(toupper(model)=="GCA"){
    M = M*2
    V = nInd*nInd*3
    W = W0 = WQX= nInd*nInd
    WX = ee = nInd
    u = u*2
  }else if(toupper(model)=="SCA"){
    M = M*3
    V = nInd*nInd*4
    W = W0 = WQX= nInd*nInd
    WX = ee = nInd
    u = u*3
  }else{
    stop(paste0("model=",toupper(model)," not recognized"))
  }
  objects = ls()
  objects = objects[objects!="model"]
  bytes = sapply(objects,function(x) get(x))
  bytes = 8*sum(bytes)
  return(bytes/1073741824) #GB
}
