#' @title Hybrid crossing
#'
#' @description
#' A convience function for hybrid plant breeding simulations. Allows for
#' easy specification of a test cross scheme and/or creation of an object
#' of \code{\link{HybridPop-class}}. Note that the \code{\link{HybridPop-class}}
#' should only be used if the parents were created using the \code{\link{makeDH}}
#' function or \code{\link{newPop}} using inbred founders. The id for
#' new individuals is [mother_id]_[father_id]
#'
#' @param females female population, an object of \code{\link{Pop-class}}
#' @param males male population, an object of \code{\link{Pop-class}}
#' @param crossPlan either "testcross" for all possible combinantions
#' or a matrix with two columns for designed crosses
#' @param returnHybridPop should results be returned as
#' \code{\link{HybridPop-class}}. If false returns results as
#' \code{\link{Pop-class}}. Population must be fully inbred if TRUE.
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Make crosses for full diallele
#' pop2 = hybridCross(pop, pop, simParam=SP)
#'
#' @export
hybridCross = function(females, males,
                       crossPlan="testcross",
                       returnHybridPop=FALSE,
                       simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if((females@ploidy%%2L != 0L) |
     (males@ploidy%%2L != 0L)){
    stop("You can not cross indiviuals with odd ploidy levels")
  }
  #crossPlan for test cross
  if(length(crossPlan)==1){
    if(crossPlan=="testcross"){
      crossPlan = cbind(rep(1:females@nInd,each=males@nInd),
                        rep(1:males@nInd,females@nInd))
    }else{
      stop(paste0("crossPlan=",crossPlan," is not a valid option"))
    }
  }

  #Set id
  femaleParents = females@id[crossPlan[,1]]
  maleParents = males@id[crossPlan[,2]]
  id = paste(femaleParents, maleParents, sep="_")

  #Return Pop-class
  if(!returnHybridPop){
    return(makeCross2(females=females,
                      males=males,
                      crossPlan=crossPlan,
                      simParam=simParam))
  }

  #Return HybridPop-class
  gv = matrix(NA_real_,
              nrow=length(id),
              ncol=simParam$nTraits)
  gxe = vector("list",simParam$nTraits)
  i = 0L
  for(trait in simParam$traits){
    i = i+1L
    tmp = getHybridGv(trait=trait,
                      females=females,
                      femaleParents=crossPlan[,1],
                      males=males,
                      maleParents=crossPlan[,2],
                      nThreads=simParam$nThreads)
    gv[,i] = tmp[[1]]
    if(length(tmp)==2){
      gxe[[i]] = tmp[[2]]
    }
  }
  if(simParam$nTraits>0){
    pheno = addError(gv, simParam$varE,
                     reps=rep(1, simParam$nTraits))
  }else{
    pheno = gv
  }
  output = new("HybridPop",
               nInd=length(id),
               id=id,
               mother=femaleParents,
               father=maleParents,
               nTraits=simParam$nTraits,
               gv=gv,
               pheno=pheno,
               gxe=gxe)
  return(output)
}

#' @title Calculate GCA
#'
#' @description
#' Calculate general combining ability of test crosses. Intended for
#' output from hybridCross using the "testcross" option, but will work
#' for any population.
#'
#' @param pop an object of \code{\link{Pop-class}} or
#' \code{\link{HybridPop-class}}
#' @param use tabulate either genetic values "gv", estimated
#' breeding values "ebv", or phenotypes "pheno"
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10, inbred=TRUE)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Make crosses for full diallele
#' pop2 = hybridCross(pop, pop, simParam=SP)
#' GCA = calcGCA(pop2, use="gv")
#'
#' @export
calcGCA = function(pop,use="pheno"){
  if(use=="pheno"){
    y = pop@pheno
  }else if(use=="gv"){
    y = pop@gv
  }else if(use=="ebv"){
    y = pop@ebv
  }else{
    stop(paste0("use=",use," is not a valid option"))
  }
  if(ncol(y)==0){
    stop(paste("No values for",use))
  }
  female = factor(pop@mother,
                  levels=unique(pop@mother))
  male = factor(pop@father,
                levels=unique(pop@father))
  #Check for balance
  if(nlevels(female)==1 | nlevels(male)==1){
    balanced = TRUE
  }else{
    tmp = table(female,male)
    if(all(tmp==tmp[1])){
      balanced = TRUE
    }else{
      balanced = FALSE
    }
  }
  sca = paste(as.character(female),as.character(male),sep="_")
  sca = factor(sca,levels=unique(sca))
  # Female GCA
  if(nlevels(female)==1){
    GCAf = matrix(colMeans(y),nrow=1)
  }else{
    if(nlevels(male)==1){
      GCAf = y
    }else{
      if(balanced){
        #Calculate simple means
        tmp = aggregate(y~female,FUN=mean)
        GCAf = unname(as.matrix(tmp[,-1,drop=F]))
      }else{
        #Calculate population marginal means
        X = model.matrix(~female+male-1,contrasts=list(male="contr.sum"))
        GCAf = calcCoef(X,y)[1:nlevels(female),,drop=FALSE]
      }
    }
  }
  GCAf = data.frame(levels(female),GCAf,
                    stringsAsFactors=FALSE)
  names(GCAf) = c("id",paste0("Trait",1:pop@nTraits))
  # Male GCA
  if(nlevels(male)==1){
    GCAm = matrix(colMeans(y),nrow=1)
  }else{
    if(nlevels(female)==1){
      GCAm = y
    }else{
      if(balanced){
        #Calculate simple means
        tmp = aggregate(y~male,FUN=mean)
        GCAm = unname(as.matrix(tmp[,-1,drop=F]))
      }else{
        #Calculate population marginal means
        X = model.matrix(~male+female-1,contrasts=list(female="contr.sum"))
        GCAm = calcCoef(X,y)[1:nlevels(male),,drop=FALSE]
      }
    }
  }
  GCAm = data.frame(levels(male),GCAm,
                    stringsAsFactors=FALSE)
  names(GCAm) = c("id",paste0("Trait",1:pop@nTraits))
  # SCA
  if(nlevels(sca)==pop@nInd){
    SCA = y
  }else{
    #Calculate simple means
    tmp = aggregate(y~sca,FUN=mean)
    SCA = unname(as.matrix(tmp[,-1,drop=F]))
  }
  SCA = data.frame(levels(sca),SCA,
                   stringsAsFactors=FALSE)
  names(SCA) = c("id",paste0("Trait",1:pop@nTraits))
  return(list(GCAf=GCAf,
              GCAm=GCAm,
              SCA=SCA))
}

#' @title Set GCA as phenotype
#'
#' @description
#' Calculates general combining ability from a set of testers and
#' returns these values as phenotypes for a population.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param testers an object of \code{\link{Pop-class}}
#' @param use true genetic value (\code{gv}) or phenotypes (\code{pheno}, default)
#' @param h2 a vector of desired narrow-sense heritabilities for
#' each trait. See details in \code{\link{setPheno}}.
#' @param H2 a vector of desired broad-sense heritabilities for
#' each trait. See details in \code{\link{setPheno}}.
#' @param varE error (co)variances for traits.
#' See details in \code{\link{setPheno}}.
#' @param corE an optional matrix for correlations between errors.
#' See details in \code{\link{setPheno}}.
#' @param reps number of replications for phenotype.
#' See details in \code{\link{setPheno}}.
#' @param fixEff fixed effect to assign to the population. Used
#' by genomic selection models only.
#' @param p the p-value for the environmental covariate
#' used by GxE traits. If NULL, a value is
#' sampled at random.
#' @param inbred are both pop and testers fully inbred. They are only
#' fully inbred if created by \code{\link{newPop}} using inbred founders
#' or by the \code{\link{makeDH}} function
#' @param onlyPheno should only the phenotype be returned, see return
#' @param simParam an object of \code{\link{SimParam}}
#'
#'
#' @return Returns an object of \code{\link{Pop-class}} or
#' a matrix if onlyPheno=TRUE
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10, inbred=TRUE)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Set phenotype to average per
#' pop2 = setPhenoGCA(pop, pop, use="gv", inbred=TRUE, simParam=SP)
#'
#' @export
setPhenoGCA = function(pop, testers, use="pheno", h2=NULL, H2=NULL,
                       varE=NULL, corE=NULL, reps=1, fixEff=1L, p=NULL,
                       inbred=FALSE, onlyPheno=FALSE, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is(pop,"MultiPop")){
    stopifnot(class(testers)=="Pop", !onlyPheno)
    pop@pops = lapply(pop@pops, setPhenoGCA, testers=testers,
                      use=use, h2=h2, H2=H2, varE=varE, corE=corE,
                      reps=reps, fixEff=fixEff, p=p, inbred=inbred,
                      onlyPheno=FALSE, simParam=simParam)
    return(pop)
  }
  if(any(duplicated(pop@id))){
    stop("This function does not work with duplicate IDs")
  }
  stopifnot(class(pop)=="Pop",class(testers)=="Pop")
  use = tolower(use)
  #Make hybrids
  tmp = hybridCross(females=pop, males=testers, crossPlan="testcross",
                    returnHybridPop=inbred, simParam=simParam)
  #Get response
  if(use=="pheno"){
    y = setPheno(tmp, h2=h2, H2=H2, varE=varE, corE=corE,
                 p=p, reps=reps, onlyPheno=TRUE, simParam=simParam)
  }else if(use=="gv"){
    y = tmp@gv
  }else{
    stop(paste0("use=",use," is not a valid option"))
  }
  if(ncol(y)==0){
    stop(paste("No values for",use))
  }
  female = factor(tmp@mother,levels=unique(tmp@mother))
  if(nlevels(female)==1){
    GCAf = matrix(colMeans(y),nrow=1)
  }else{
    if(testers@nInd==1){
      GCAf = y
    }else{
        #Calculate simple means
        tmp = aggregate(y~female,FUN=mean)
        GCAf = unname(as.matrix(tmp[,-1,drop=F]))
    }
  }
  if(onlyPheno){
    return(GCAf)
  }
  pop@pheno = GCAf
  pop@fixEff = rep(as.integer(fixEff),pop@nInd)
  return(pop)
}

#' @title Set progeny test as phenotype
#'
#' @description
#' Models a progeny test of individuals in 'pop'. Returns 'pop' with a phenotype
#' representing the average performance of their progeny. The phenotype is generated
#' by mating individuals in 'pop' to randomly chosen individuals in testPop a
#' number of times equal to 'nMatePerInd'.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param testPop an object of \code{\link{Pop-class}}
#' @param nMatePerInd number of times an individual in 'pop' is mated to an
#' individual in testPop
#' @param use true genetic value (\code{gv}) or phenotypes (\code{pheno}, default)
#' @param h2 a vector of desired narrow-sense heritabilities for
#' each trait. See details in \code{\link{setPheno}}.
#' @param H2 a vector of desired broad-sense heritabilities for
#' each trait. See details in \code{\link{setPheno}}.
#' @param varE error (co)variances for traits.
#' See details in \code{\link{setPheno}}.
#' @param corE an optional matrix for correlations between errors.
#' See details in \code{\link{setPheno}}.
#' @param reps number of replications for phenotype.
#' See details in \code{\link{setPheno}}.
#' @param fixEff fixed effect to assign to the population. Used
#' by genomic selection models only.
#' @param p the p-value for the environmental covariate
#' used by GxE traits. If NULL, a value is
#' sampled at random.
#' @param onlyPheno should only the phenotype be returned, see return
#' @param simParam an object of \code{\link{SimParam}}
#'
#' @details
#' The reps parameter is for convenient representation of replicated data.
#' It was intended for representation of replicated yield trials in plant
#' breeding programs. In this case, varE is set to the plot error and
#' reps is set to the number plots per entry. The resulting phenotype
#' would reflect the mean of all replications.
#'
#' @return Returns an object of \code{\link{Pop-class}} or
#' a matrix if onlyPheno=TRUE
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10, inbred=TRUE)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#'
#' #Create two populations of 5 individuals
#' pop1 = newPop(founderPop[1:5], simParam=SP)
#' pop2 = newPop(founderPop[6:10], simParam=SP)
#'
#' #Set phenotype according to a progeny test
#' pop3 = setPhenoProgTest(pop1, pop2, use="gv", simParam=SP)
#'
#' @export
setPhenoProgTest = function(pop, testPop, nMatePerInd=1L, use="pheno",
                            h2=NULL, H2=NULL, varE=NULL, corE=NULL,
                            reps=1, fixEff=1L, p=NULL, onlyPheno=FALSE,
                            simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is(pop,"MultiPop")){
    stopifnot(class(testPop)=="Pop", !onlyPheno)
    pop@pops = lapply(pop@pops, setPhenoProgTest, testPop=testPop,
                      nMatePerInd=nMatePerInd, use=use, h2=h2, H2=H2,
                      varE=varE, corE=corE, reps=reps, fixEff=fixEff,
                      p=p, onlyPheno=FALSE, simParam=simParam)
    return(pop)
  }
  if(any(duplicated(pop@id))){
    stop("This function does not work with duplicate IDs")
  }
  stopifnot(class(pop)=="Pop",class(testPop)=="Pop")
  use = tolower(use)
  #Make hybrids
  tmp = randCross2(females=pop, males=testPop, nCrosses=nInd(pop)*nMatePerInd,
                   balance=TRUE, simParam=simParam)
  #Get response
  if(use=="pheno"){
    y = setPheno(tmp, h2=h2, H2=H2, varE=varE, corE=corE,
                 reps=reps, p=p, onlyPheno=TRUE, simParam=simParam)
  }else if(use=="gv"){
    y = tmp@gv
  }else{
    stop(paste0("use=",use," is not a valid option"))
  }
  if(ncol(y)==0){
    stop(paste("No values for",use))
  }
  female = factor(tmp@mother, levels=pop@id)
  #Calculate simple means
  tmp = aggregate(y~female, FUN=mean)
  GCAf = unname(as.matrix(tmp[,-1,drop=F]))

  if(onlyPheno){
    return(GCAf)
  }
  pop@pheno = GCAf
  pop@fixEff = rep(as.integer(fixEff),pop@nInd)
  return(pop)
}
