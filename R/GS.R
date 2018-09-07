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
                        lociLoc, maxIter)
  }else{
    ans = callRRBLUP(y,fixEff,pop@reps,pop@geno,lociPerChr,lociLoc)
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
               fixEff=ans$beta,
               Vu=as.matrix(ans$Vu),
               Ve=as.matrix(ans$Ve),
               LL=ans$LL,
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
#' The initial values are considered true.
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
                    lociLoc,Vu,Ve,tol,maxIter,useEM)
  if(is.null(ans[["iter"]])){
    iter = 0
  }else{
    iter = ans$iter
  }
  output = new("RRsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               markerEff=ans$u,
               fixEff=ans$beta,
               Vu=as.matrix(ans$Vu),
               Ve=as.matrix(ans$Ve),
               LL=numeric(),
               iter=iter)
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
#' @param useHetCov should the model include a covariate 
#' for heterozygosity.
#' @param maxIter maximum number of iterations. Only used 
#' when number of traits is greater than 1.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @export
RRBLUP_D = function(pop, traits=1, use="pheno", snpChip=1, 
                    useQtl=FALSE, useHetCov=TRUE, maxIter=40L, 
                    simParam=NULL, ...){
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
                     lociLoc,maxIter,useHetCov)
  p = t(ans$p)
  q = 1-p
  ans = ans$ans
  fixEff=ans$beta
  a = ans$u[[1]]
  d = ans$u[[2]]
  if(useHetCov){
    stopifnot(length(fixEff)>1)
    hetCov = fixEff[length(fixEff)]
    fixEff = matrix(fixEff[-length(fixEff)])
    alpha = a+(q-p)*(d+hetCov/nLoci)
  }else{
    hetCov = 0
    alpha = a+(q-p)*d
  }
  output = new("RRDsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               markerEff=alpha,
               addEff=a,
               domEff=d,
               hetCov=hetCov,
               fixEff=fixEff,
               Vu=ans$Vu,
               Ve=ans$Ve,
               LL=ans$LL,
               iter=ans$iter)
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
                       lociPerChr,lociLoc,maxIter)
  output = new("GCAsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               femaleEff=ans$u[[1]],
               maleEff=ans$u[[2]],
               fixEff=ans$beta,
               Vu=ans$Vu,
               Ve=ans$Ve,
               LL=ans$LL,
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
#' @param useHetCov should the model include a covariate 
#' for heterozygosity.
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' traits
#'
#' @export
RRBLUP_SCA = function(pop, traits=1, use="pheno", snpChip=1, 
                      useQtl=FALSE, maxIter=40L, useHetCov=FALSE, 
                      simParam=NULL, ...){
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
                       lociPerChr,lociLoc,maxIter,useHetCov)
  p1 = t(ans$p1)
  q1 = 1-p1
  p2 = t(ans$p2)
  q2 = 1-p2
  ans = ans$ans
  fixEff=ans$beta
  a1 = ans$u[[1]]
  a2 = ans$u[[2]]
  d = ans$u[[3]]
  if(useHetCov){
    stopifnot(length(fixEff)>1)
    hetCov = fixEff[length(fixEff)]
    fixEff = matrix(fixEff[-length(fixEff)])
    alpha1 = a1+(q2-p2)*(d+hetCov/nLoci)
    alpha2 = a2+(q1-p1)*(d+hetCov/nLoci)
  }else{
    hetCov = 0
    alpha1 = a1+(q2-p2)*d
    alpha2 = a2+(q1-p1)*d
  }
  output = new("SCAsol",
               nLoci=nLoci,
               lociPerChr=lociPerChr,
               lociLoc=lociLoc,
               femaleEff=alpha1,
               maleEff=alpha2,
               a1=ans$u[[1]],
               a2=ans$u[[2]],
               d=ans$u[[3]],
               hetCov=hetCov,
               fixEff=fixEff,
               Vu=ans$Vu,
               Ve=ans$Ve,
               LL=ans$LL,
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
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @export
setEBV = function(pop, solution, gender=NULL, useGV=FALSE, 
                  append=FALSE){
  if(class(solution)=="RRsol"){
    ebv = gebvRR(solution, pop)
  }else if(class(solution)=="RRDsol"){
    if(useGV){
      ebv = gegvRRD(solution, pop)
    }else{
      ebv = gebvRR(solution, pop)
    }
  }else{
    if(is.null(gender)){
      if(class(solution)=="GCAsol"){
        ebv = gegvGCA(solution, pop)
      }else{
        ebv = gegvSCA(solution, pop)
      }
    }else if(toupper(gender)=="FEMALE"){
      ebv = gebvGCA(solution, pop, TRUE)
    }else if(toupper(gender)=="MALE"){
      ebv = gebvGCA(solution, pop, FALSE)
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
  return(bytes*1e-9) #GB
}
