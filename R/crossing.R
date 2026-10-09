#' @title Make designed crosses
#'
#' @description
#' Makes crosses within a population using a user supplied
#' crossing plan.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param crossPlan a matrix with two columns representing
#' female and male parents. Either integers for the position in
#' population or character strings for the IDs.
#' @param nProgeny number of progeny per cross. May be a single value for all 
#' crosses or a vector with values for each cross.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @family mating functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Cross individual 1 with individual 10 and 2 with 4
#' crossPlan = matrix(c(1,10,
#'                      2,4),
#'                    nrow=2, ncol=2, byrow=TRUE)
#' pop2 = makeCross(pop, crossPlan, simParam=SP)
#' getPed(pop2)
#'
#' #The same but variable nProgeny
#' pop3 = makeCross(pop, crossPlan, nProgeny=c(1,2),simParam=SP)
#' getPed(pop3)
#' @export
makeCross = function(pop, crossPlan, nProgeny=1,
                     simParam=NULL, nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  if(pop@ploidy%%2L != 0L){
    stop("You can not cross indiviuals with odd ploidy levels")
  }
  
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1], pop@id),
                      match(crossPlan[,2], pop@id))
    if(any(is.na(crossPlan))){
      stop("Failed to match supplied IDs")
    }
  }
  
  if((max(crossPlan)>nInd(pop)) |
     (min(crossPlan)<1L)){
    stop("Invalid crossPlan")
  }
  
  # Handle nProgeny
  if(length(nProgeny)==1){
    # Any value other than 1 has to go through rep, including 0, which
    # must yield an empty cross plan rather than one progeny per cross.
    if(nProgeny!=1){
      crossPlan = cbind(rep(crossPlan[,1], each=nProgeny),
                        rep(crossPlan[,2], each=nProgeny))
    }
  }else{
    if(nrow(crossPlan)!=length(nProgeny)){
      stop("Length of nProgeny must equal 1 or nrow(crossPlan)")
    }
    
    crossPlan = cbind(rep(crossPlan[,1], times=nProgeny),
                      rep(crossPlan[,2], times=nProgeny))
  }
  
  tmp = cross(pop@geno,
              crossPlan[,1],
              pop@geno,
              crossPlan[,2],
              simParam$femaleMap,
              simParam$maleMap,
              simParam$isTrackRec,
              pop@ploidy,
              pop@ploidy,
              simParam$v,
              simParam$p,
              simParam$femaleCentromere,
              simParam$maleCentromere,
              simParam$quadProb,
              nThreads)
  
  dim(tmp$geno) = NULL # Account for matrix bug in RcppArmadillo
  
  rPop = new("RawPop",
             nInd=nrow(crossPlan),
             nChr=pop@nChr,
             ploidy=pop@ploidy,
             nLoci=pop@nLoci,
             geno=tmp$geno)
  
  if(simParam$isTrackRec){
    hist = tmp$recHist
  }else{
    hist = NULL
  }
  
  return(.newPop(rawPop=rPop,
                 mother=pop@id[crossPlan[,1]],
                 father=pop@id[crossPlan[,2]],
                 iMother=pop@iid[crossPlan[,1]],
                 iFather=pop@iid[crossPlan[,2]],
                 femaleParentPop=pop,
                 maleParentPop=pop,
                 hist=hist,
                 simParam=simParam,
                 nThreads=nThreads))
}

#' @title Make random crosses
#'
#' @description
#' A wrapper for \code{\link{makeCross}} that randomly
#' selects parental combinations for all possible combinations.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross. May be a single value for all 
#' crosses or a vector with values equal to the number of crosses. If providing 
#' a vector, the values are randomly assigned to each cross.
#' @param balance if using sexes, this option will balance the number
#' of progeny per parent
#' @param parents an optional vector of indices for allowable parents
#' @param ignoreSexes should sexes be ignored
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @family mating functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Make 10 crosses
#' pop2 = randCross(pop, 10, simParam=SP)
#'
#' @export
randCross = function(pop, nCrosses, nProgeny=1,
                     balance=TRUE, parents=NULL,
                     ignoreSexes=FALSE,
                     simParam=NULL, nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  if(is.null(parents)){
    parents = 1:pop@nInd
  }else{
    parents = as.integer(parents)
  }
  
  n = length(parents)
  if(n<=1){
    stop("The population must contain more than 1 individual")
  }
  
  # Handle nProgeny
  if(length(nProgeny)>1){
    if(nCrosses!=length(nProgeny)){
      stop("Length of nProgeny must equal 1 or nCrosses")
    }
    nProgeny = nProgeny[sample(nCrosses, nCrosses)]
  }
  
  if(simParam$sexes=="no" | ignoreSexes){
    crossPlan = sampHalfDialComb(n, nCrosses)
    crossPlan[,1] = parents[crossPlan[,1]]
    crossPlan[,2] = parents[crossPlan[,2]]
  }else{
    female = which(pop@sex=="F" & (1:pop@nInd)%in%parents)
    nFemale = length(female)
    if(nFemale==0){
      stop("population doesn't contain any females")
    }
    male = which(pop@sex=="M" & (1:pop@nInd)%in%parents)
    nMale = length(male)
    if(nMale==0){
      stop("population doesn't contain any males")
    }
    if(balance){
      female = female[sample.int(nFemale, nFemale)]
      female = rep(female, length.out=nCrosses)
      tmp = male[sample.int(nMale, nMale)]
      n = nCrosses%/%nMale + 1
      male = NULL
      for(i in seq_len(n)){
        take = nMale - (i:(nMale+i-1))%%nMale
        male = c(male, tmp[take])
      }
      male = male[1:nCrosses]
      crossPlan = cbind(female,male)
    }else{
      crossPlan = sampAllComb(nFemale,
                              nMale,
                              nCrosses)
      crossPlan[,1] = female[crossPlan[,1]]
      crossPlan[,2] = male[crossPlan[,2]]
    }
  }
  
  return(makeCross(pop=pop, crossPlan=crossPlan, nProgeny=nProgeny,
                   simParam=simParam, nThreads=nThreads))
}

#' @title Select and randomly cross
#'
#' @description
#' This is a wrapper that combines the functionalities of
#' \code{\link{randCross}} and \code{\link{selectInd}}. The
#' purpose of this wrapper is to combine both selection and
#' crossing in one function call that minimizes the amount
#' of intermediate populations created. This reduces RAM usage
#' and simplifies code writing. Note that this wrapper does not
#' provide the full functionality of either function.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param nInd the number of individuals to select. These individuals
#' are selected without regard to sex and it supersedes values
#' for nFemale and nMale. Thus if the simulation uses sexes, it is
#' likely better to leave this value as NULL and use nFemale and nMale
#' instead.
#' @param nFemale the number of females to select. This value is ignored
#' if nInd is set.
#' @param nMale the number of males to select. This value is ignored
#' if nInd is set.
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross
#' @param trait the trait for selection. Either a number indicating
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values "gv", estimated
#' breeding values "ebv", breeding values "bv", phenotypes "pheno",
#' or randomly "rand"
#' @param selectTop selects highest values if true.
#' Selects lowest values if false.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param ... additional arguments if using a function for
#' trait
#' @param balance if using sexes, this option will balance the number
#' of progeny per parent. This argument occurs after ..., so the argument
#' name must be matched exactly.
#' @param restrInbr should selection restrict the expected increase
#' in fixation, using \code{\link{selectOCS}}. If FALSE, parents are
#' selected by truncation on \code{trait} and \code{use}. This argument
#' occurs after ..., so the argument name must be matched exactly.
#' @param inbrTarget the target for expected fixation when
#' \code{restrInbr=TRUE}. Its meaning depends on \code{inbrType}. The
#' default of 0.01 with \code{inbrType="relative"} allows a loss of 1\%
#' of the current heterozygosity per generation, which corresponds to an
#' effective population size of 50 for a diploid. See
#' \code{\link{selectOCS}}. This argument occurs after ..., so the
#' argument name must be matched exactly.
#' @param inbrType either "relative", where \code{inbrTarget} is the
#' allowed increase in expected fixation as a proportion of the
#' remaining heterozygosity in \code{pop}, or "absolute", where
#' \code{inbrTarget} is the maximum allowed expected fixation. This
#' argument occurs after ..., so the argument name must be matched
#' exactly.
#' @param snpChip an integer indicating which SNP chip genotypes are
#' used to measure expected fixation when \code{restrInbr=TRUE}. This
#' argument occurs after ..., so the argument name must be matched
#' exactly.
#' @param useQtl should QTL genotypes be used instead of a SNP chip
#' to measure expected fixation. If TRUE, snpChip specifies which
#' trait's QTL to use. This argument occurs after ..., so the argument
#' name must be matched exactly.
#'
#' @details
#' When \code{restrInbr=TRUE}, parents are selected with
#' \code{\link{selectOCS}}, an approximation to optimal contribution
#' selection that maximizes merit while restricting the expected
#' increase in fixation. They are then crossed at random, exactly as
#' when \code{restrInbr=FALSE}. See \code{\link{selectOCS}} for the
#' method and for how \code{inbrTarget} and \code{inbrType} set the
#' target.
#'
#' The method assumes every selected parent contributes equally to the
#' next generation. Setting \code{nInd} when the simulation uses sexes
#' only approximates this, because crosses must be between sexes, so a
#' warning is given.
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @family mating functions
#' @family selection functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' SP$addTraitA(10)
#' SP$setVarE(h2=0.5)
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Select 4 individuals and make 8 crosses
#' pop2 = selectCross(pop, nInd=4, nCrosses=8, simParam=SP)
#'
#' #Select 4 individuals while restricting the increase in expected
#' #fixation measured at the QTL, and make 8 crosses
#' pop3 = selectCross(pop, nInd=4, nCrosses=8, restrInbr=TRUE,
#'                    inbrTarget=0.5, useQtl=TRUE, simParam=SP)
#'
#' @export
selectCross = function(pop, nInd=NULL, nFemale=NULL, nMale=NULL, nCrosses,
                       nProgeny=1, trait=1, use="pheno", selectTop=TRUE,
                       simParam=NULL, nThreads=NULL, ..., balance=TRUE,
                       restrInbr=FALSE, inbrTarget=0.01,
                       inbrType="relative", snpChip=1, useQtl=FALSE){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  if(restrInbr){
    # randCross crosses only between sexes, so the parents selected
    # with nInd cannot all contribute equally as selectOCS assumes
    if(!is.null(nInd) && simParam$sexes!="no"){
      warning("restrInbr with nInd assumes each parent contributes equally, ",
              "which is only approximate when crosses must be between sexes; ",
              "consider nFemale and nMale instead")
    }
    parents = selectOCS(pop=pop, nInd=nInd, nFemale=nFemale,
                        nMale=nMale, trait=trait, use=use,
                        selectTop=selectTop, inbrTarget=inbrTarget,
                        inbrType=inbrType, snpChip=snpChip,
                        useQtl=useQtl, returnPop=FALSE,
                        simParam=simParam, nThreads=nThreads, ...)
  }else if(!is.null(nInd)){
    parents = selectInd(pop=pop, nInd=nInd, trait=trait, use=use,
                        sex="B", selectTop=selectTop,
                        returnPop=FALSE, simParam=simParam,
                        nThreads=nThreads, ...)
  }else{
    if(simParam$sexes=="no")
      stop("You must specify nInd when simParam$sexes is `no`")
    if(is.null(nFemale))
      stop("You must specify nFemale if nInd is NULL")
    if(is.null(nMale))
      stop("You must specify nMale if nInd is NULL")
    females = selectInd(pop=pop, nInd=nFemale, trait=trait, use=use,
                        sex="F", selectTop=selectTop,
                        returnPop=FALSE, simParam=simParam,
                        nThreads=nThreads, ...)
    males = selectInd(pop=pop, nInd=nMale, trait=trait, use=use,
                      sex="M", selectTop=selectTop,
                      returnPop=FALSE, simParam=simParam,
                      nThreads=nThreads, ...)
    parents = c(females,males)
  }
  
  return(randCross(pop=pop, nCrosses=nCrosses, nProgeny=nProgeny,
                   balance=balance, parents=parents,
                   ignoreSexes=FALSE, simParam=simParam,
                   nThreads=nThreads))
}

#' @title Make designed crosses
#'
#' @description
#' Makes crosses between two populations using a user supplied
#' crossing plan.
#'
#' @param females an object of \code{\link{Pop-class}} for female parents.
#' @param males an object of \code{\link{Pop-class}} for male parents.
#' @param crossPlan a matrix with two columns representing
#' female and male parents. Either integers for the position in
#' population or character strings for the IDs.
#' @param nProgeny number of progeny per cross. May be a single value for all 
#' crosses or a vector with values for each cross.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @family mating functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Cross individual 1 with individual 10 and 2 with 4
#' crossPlan = matrix(c(1,10,
#'                      2,4),
#'                    nrow=2, ncol=2, byrow=TRUE)
#' pop2 = makeCross2(pop, pop, crossPlan, simParam=SP)
#' getPed(pop2)
#'
#' #The same but variable nProgeny
#' pop3 = makeCross2(pop, pop, crossPlan, nProgeny=c(1,2),simParam=SP)
#' getPed(pop3)
#' @export
makeCross2 = function(females, males, crossPlan, nProgeny=1, simParam=NULL,
                      nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  if((females@ploidy%%2L != 0L) |
     (males@ploidy%%2L != 0L)){
    stop("You can not cross indiviuals with odd ploidy levels")
  }
  
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],females@id),
                      match(crossPlan[,2],males@id))
    if(any(is.na(crossPlan))){
      stop("Failed to match supplied IDs")
    }
  }
  
  if((max(crossPlan[,1])>nInd(females)) |
     (max(crossPlan[,2])>nInd(males)) |
     (min(crossPlan)<1L)){
    stop("Invalid crossPlan")
  }
  
  # Handle nProgeny
  if(length(nProgeny)==1){
    # Any value other than 1 has to go through rep, including 0, which
    # must yield an empty cross plan rather than one progeny per cross.
    if(nProgeny!=1){
      crossPlan = cbind(rep(crossPlan[,1], each=nProgeny),
                        rep(crossPlan[,2], each=nProgeny))
    }
  }else{
    if(nrow(crossPlan)!=length(nProgeny)){
      stop("Length of nProgeny must equal 1 or nrow(crossPlan)")
    }
    
    crossPlan = cbind(rep(crossPlan[,1], times=nProgeny),
                      rep(crossPlan[,2], times=nProgeny))
  }
  
  tmp=cross(females@geno,
            crossPlan[,1],
            males@geno,
            crossPlan[,2],
            simParam$femaleMap,
            simParam$maleMap,
            simParam$isTrackRec,
            females@ploidy,
            males@ploidy,
            simParam$v,
            simParam$p,
            simParam$femaleCentromere,
            simParam$maleCentromere,
            simParam$quadProb,
            nThreads)
  
  dim(tmp$geno) = NULL # Account for matrix bug in RcppArmadillo
  
  rPop = new("RawPop",
             nInd=nrow(crossPlan),
             nChr=females@nChr,
             ploidy=as.integer((females@ploidy+males@ploidy)/2),
             nLoci=females@nLoci,
             geno=tmp$geno)
  
  if(simParam$isTrackRec){
    hist = tmp$recHist
  }else{
    hist = NULL
  }
  
  return(.newPop(rawPop=rPop,
                 mother=females@id[crossPlan[,1]],
                 father=males@id[crossPlan[,2]],
                 iMother=females@iid[crossPlan[,1]],
                 iFather=males@iid[crossPlan[,2]],
                 femaleParentPop=females,
                 maleParentPop=males,
                 hist=hist,
                 simParam=simParam,
                 nThreads=nThreads))
}

#' @title Make random crosses
#'
#' @description
#' A wrapper for \code{\link{makeCross2}} that randomly
#' selects parental combinations for all possible combinations between
#' two populations.
#'
#' @param females an object of \code{\link{Pop-class}} for female parents.
#' @param males an object of \code{\link{Pop-class}} for male parents.
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross. May be a single value for all 
#' crosses or a vector with values equal to the number of crosses. If providing 
#' a vector, the values are randomly assigned to each cross.
#' @param balance this option will balance the number
#' of progeny per parent
#' @param femaleParents an optional vector of indices for allowable
#' female parents
#' @param maleParents an optional vector of indices for allowable
#' male parents
#' @param ignoreSexes should sex be ignored
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @family mating functions
#'
#' @examples
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
#'
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#'
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#'
#' #Make 10 crosses
#' pop2 = randCross2(pop, pop, 10, simParam=SP)
#'
#' @export
randCross2 = function(females, males, nCrosses, nProgeny=1,
                      balance=TRUE, femaleParents=NULL,
                      maleParents=NULL, ignoreSexes=FALSE,
                      simParam=NULL, nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  #Set allowable parents
  if(is.null(femaleParents)){
    femaleParents = 1:females@nInd
  }else{
    femaleParents = as.integer(femaleParents)
  }
  
  if(is.null(maleParents)){
    maleParents = 1:males@nInd
  }else{
    maleParents = as.integer(maleParents)
  }
  
  if(simParam$sexes=="no" | ignoreSexes){
    female = femaleParents
    male = maleParents
  }else{
    female = which(females@sex=="F" &
                     (1:females@nInd)%in%femaleParents)
    if(length(female)==0){
      stop("population doesn't contain any females")
    }
    male = which(males@sex=="M" &
                   (1:males@nInd)%in%maleParents)
    if(length(male)==0){
      stop("population doesn't contain any males")
    }
  }
  
  # Handle nProgeny
  if(length(nProgeny)>1){
    if(nCrosses!=length(nProgeny)){
      stop("Length of nProgeny must equal 1 or nCrosses")
    }
    nProgeny = nProgeny[sample(nCrosses, nCrosses)]
  }
  
  nMale = length(male)
  nFemale = length(female)
  
  if(balance){
    female = female[sample.int(nFemale, nFemale)]
    female = rep(female, length.out=nCrosses)
    tmp = male[sample.int(nMale, nMale)]
    n = nCrosses%/%nMale + 1
    male = NULL
    for(i in seq_len(n)){
      take = nMale - (i:(nMale+i-1))%%nMale
      male = c(male, tmp[take])
    }
    male = male[1:nCrosses]
    crossPlan = cbind(female,male)
  }else{
    crossPlan = sampAllComb(nFemale,
                            nMale,
                            nCrosses)
    crossPlan[,1] = female[crossPlan[,1]]
    crossPlan[,2] = male[crossPlan[,2]]
  }
  
  return(makeCross2(females=females, males=males,
                    crossPlan=crossPlan, nProgeny=nProgeny,
                    simParam=simParam, nThreads=nThreads))
}

#' @title Self individuals
#'
#' @description
#' Creates selfed progeny from each individual in a
#' population. Only works when sexes is "no".
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param nProgeny number of selfed progeny per individual. May be a single value 
#' for all or a vector providing values for each individual.
#' @param parents an optional vector of indices for allowable parents
#' @param keepParents should previous parents be used for mother and
#' father.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @family mating functions
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
#' #Self pollinate each individual
#' pop2 = self(pop, simParam=SP)
#'
#' @export
self = function(pop, nProgeny=1, parents=NULL, keepParents=TRUE,
                simParam=NULL, nThreads=NULL){
  
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  if(is(pop,"MultiPop")){
    if(!is.null(parents)) stop("parents must be NULL for a MultiPop")
    pop@pops = lapply(pop@pops, self, nProgeny=nProgeny,
                      parents=NULL, keepParents=keepParents,
                      simParam=simParam, nThreads=nThreads)
    return(pop)
  }
  
  if(is.null(parents)){
    parents = 1:pop@nInd
  }else{
    parents = as.integer(parents)
    # The C++ code indexes the genotypes with these directly, so an index
    # outside the population would read memory beyond them
    if(any(is.na(parents)) ||
       any(parents<1L) ||
       any(parents>pop@nInd)){
      stop("Invalid parents")
    }
  }

  if(pop@ploidy%%2L != 0L){
    stop("You can not self aneuploids")
  }
  
  # Handle nProgeny
  if(length(nProgeny)==1){
    crossPlan = rep(parents, each=nProgeny)
    
  }else{
    # crossPlan repeats parents, which may be a subset of the population,
    # so nProgeny is checked against parents and not the population size.
    if(length(parents)!=length(nProgeny)){
      stop("Length of nProgeny must equal 1 or length(parents)")
    }
    
    crossPlan = rep(parents, times=nProgeny)
  }
  
  crossPlan = cbind(crossPlan,crossPlan)
  
  tmp = cross(pop@geno,
              crossPlan[,1],
              pop@geno,
              crossPlan[,2],
              simParam$femaleMap,
              simParam$maleMap,
              simParam$isTrackRec,
              pop@ploidy,
              pop@ploidy,
              simParam$v,
              simParam$p,
              simParam$femaleCentromere,
              simParam$maleCentromere,
              simParam$quadProb,
              nThreads)
  
  dim(tmp$geno) = NULL # Account for matrix bug in RcppArmadillo
  
  rPop = new("RawPop",
             nInd=nrow(crossPlan),
             nChr=pop@nChr,
             ploidy=pop@ploidy,
             nLoci=pop@nLoci,
             geno=tmp$geno)
  
  if(simParam$isTrackRec){
    hist = tmp$recHist
  }else{
    hist = NULL
  }
  
  if(keepParents){
    return(.newPop(rawPop=rPop,
                   mother=pop@mother[crossPlan[,1]],
                   father=pop@father[crossPlan[,1]],
                   iMother=pop@iid[crossPlan[,1]],
                   iFather=pop@iid[crossPlan[,1]],
                   femaleParentPop=pop,
                   maleParentPop=pop,
                   hist=hist,
                   simParam=simParam,
                   nThreads=nThreads))
  }else{
    return(.newPop(rawPop=rPop,
                   mother=pop@id[crossPlan[,1]],
                   father=pop@id[crossPlan[,1]],
                   iMother=pop@iid[crossPlan[,1]],
                   iFather=pop@iid[crossPlan[,1]],
                   femaleParentPop=pop,
                   maleParentPop=pop,
                   hist=hist,
                   simParam=simParam,
                   nThreads=nThreads))
  }
}

#' @title Generates DH lines
#'
#' @description Creates DH lines from each individual in a population.
#' Only works with diploid individuals. For polyploids, use
#' \code{\link{reduceGenome}} and \code{\link{doubleGenome}}.
#'
#' @param pop an object of 'Pop' superclass
#' @param nDH total number of DH lines per individual. May be a single
#' value for all individuals or a vector with values for each individual.
#' A value of zero produces no DH lines for that individual, and if no
#' individual produces any the function returns an empty population.
#' @param useFemale should female recombination rates be used.
#' @param keepParents should previous parents be used for mother and
#' father.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#'
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @family mating functions
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
#' #Create 1 DH for each individual
#' pop2 = makeDH(pop, simParam=SP)
#'
#' @export
makeDH = function(pop, nDH=1, useFemale=TRUE, keepParents=TRUE,
                  simParam=NULL, nThreads=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(is.null(nThreads)){
    nThreads = simParam$nThreads
  }else{
    nThreads = as.integer(nThreads)
  }
  
  if(is(pop,"MultiPop")){
    # Each population in a MultiPop has its own number of individuals, so
    # a vector of values could only be right for one of them
    if(length(nDH)>1){
      stop("nDH must be a single value for a MultiPop")
    }
    pop@pops = lapply(pop@pops, makeDH, nDH=nDH, useFemale=useFemale,
                      keepParents=keepParents, simParam=simParam,
                      nThreads=nThreads)
    return(pop)
  }
  
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  
  # Handle nDH. It is expanded to one value per individual here so that
  # everything below can treat the single value and the vector alike.
  if(length(nDH)==1){
    nDH = rep(nDH, pop@nInd)
  }else if(length(nDH)!=pop@nInd){
    stop("Length of nDH must equal 1 or nInd(pop)")
  }
  nDH = as.integer(nDH)
  if(anyNA(nDH) || any(nDH<0L)){
    stop("nDH must be a non-negative integer")
  }
  
  if(sum(nDH)==0L){
    return(newEmptyPop(ploidy=pop@ploidy, simParam=simParam))
  }
  
  if(useFemale){
    tmp = createDH2(pop@geno, nDH,
                    simParam$femaleMap,
                    simParam$v,
                    simParam$p,
                    simParam$isTrackRec,
                    nThreads)
  }else{
    tmp = createDH2(pop@geno, nDH,
                    simParam$maleMap,
                    simParam$v,
                    simParam$p,
                    simParam$isTrackRec,
                    nThreads)
  }
  
  dim(tmp$geno) = NULL # Account for matrix bug in RcppArmadillo
  
  rPop = new("RawPop",
             nInd=as.integer(sum(nDH)),
             nChr=pop@nChr,
             ploidy=pop@ploidy,
             nLoci=pop@nLoci,
             geno=tmp$geno)
  
  if(simParam$isTrackRec){
    hist = tmp$recHist
  }else{
    hist = NULL
  }
  
  if(keepParents){
    return(.newPop(rawPop=rPop,
                   mother=rep(pop@mother, times=nDH),
                   father=rep(pop@father, times=nDH),
                   isDH=TRUE,
                   iMother=rep(pop@iid, times=nDH),
                   iFather=rep(pop@iid, times=nDH),
                   femaleParentPop=pop,
                   maleParentPop=pop,
                   hist=hist,
                   simParam=simParam,
                   nThreads=nThreads))
  }else{
    return(.newPop(rawPop=rPop,
                   mother=rep(pop@id, times=nDH),
                   father=rep(pop@id, times=nDH),
                   isDH=TRUE,
                   iMother=rep(pop@iid, times=nDH),
                   iFather=rep(pop@iid, times=nDH),
                   femaleParentPop=pop,
                   maleParentPop=pop,
                   hist=hist,
                   simParam=simParam,
                   nThreads=nThreads))
  }
}


# Assign a generation number to every individual of a pedigree
#
# id, id of individual
# mother, name of individual's mother, NA if unknown
# father, name of individual's father, NA if unknown
# maxCycle, the most passes to make over the pedigree, or NULL to let the
#   pedigree decide
#
# Unknown parents are expected to arrive as NA. pedigreeCross normalizes the
# codes for an unknown parent before calling this.
#
# A generation number is one more than the larger of its parents', counting
# an unknown parent as zero, so a founder is generation one. That makes it
# the longest path back to a founder, and it depends only on the parents and
# not on the order the pedigree came in.
#
# Each pass assigns a number to everyone whose parents already have one, so
# the number of passes needed is the depth of the pedigree rather than its
# size. The work of a pass is vectorised over the whole pedigree, which
# keeps the interpreted part of the cost proportional to the depth alone.
#
# A pass that assigns nothing leaves the state exactly as it was, so no
# later pass could assign anything either. That is a cycle, and it is worth
# telling apart from simply not having been given enough passes.
#
# A pedigree of n individuals is at most n generations deep, which is the
# chain in which everyone is the child of the one before. That is what
# maxCycle of NULL uses, so the bound is never the thing that stops a
# pedigree being sorted. It is kept as a bound on the work rather than as a
# test of the pedigree, because a pass that assigns nobody already detects a
# cycle exactly, and does so on the pass it first stalls on.
#
# Returns an integer vector of generation numbers, one per individual.
sortPed = function(id, mother, father, maxCycle=NULL){
  nInd = length(id)
  if(is.null(maxCycle)){
    maxCycle = nInd
  }else{
    maxCycle = as.integer(maxCycle)
  }
  motherRow = match(mother, id)
  fatherRow = match(father, id)

  # Zero marks an individual that has not been given a generation yet, so a
  # parent still sitting at zero is a parent that is not ready
  gen = integer(nInd)
  gen[is.na(motherRow) & is.na(fatherRow)] = 1L

  hasMother = !is.na(motherRow)
  hasFather = !is.na(fatherRow)
  motherRow[!hasMother] = 1L
  fatherRow[!hasFather] = 1L

  for(pass in seq_len(maxCycle)){
    nLeft = sum(gen==0L)
    if(nLeft==0L){
      break
    }

    # An unknown parent contributes zero, which is what makes a half
    # founder one generation on from its one known parent
    motherGen = fatherGen = integer(nInd)
    motherGen[hasMother] = gen[motherRow[hasMother]]
    fatherGen[hasFather] = gen[fatherRow[hasFather]]

    ready = gen==0L &
      (!hasMother | motherGen>0L) &
      (!hasFather | fatherGen>0L)
    gen[ready] = pmax(motherGen[ready], fatherGen[ready]) + 1L

    if(sum(gen==0L)==nLeft){
      stop("Pedigree contains a cycle involving: ",
           paste(id[gen==0L], collapse=", "))
    }
  }

  if(any(gen==0L)){
    stop("Failed to sort pedigree within maxCycle=", maxCycle,
         " passes. Unsorted individuals: ",
         paste(id[gen==0L], collapse=", "),
         ". Try increasing maxCycle.")
  }

  return(gen)
}

# Check the recombination settings given through ... by a crossing function
# that accepts a map population, such as pedigreeCross or hybridCross
#
# dots, the list of ... arguments
#
# Returns dots unchanged. Only v, p and quadProb are allowed: ... would
# otherwise swallow a misspelled argument name without a word.
checkRecombArgs = function(dots){
  recombArgs = c("v", "p", "quadProb")
  if(length(dots)==0L){
    return(dots)
  }
  nm = names(dots)
  if(is.null(nm)){
    nm = rep("", length(dots))
  }
  bad = c(rep("<unnamed>", sum(nm=="")), setdiff(nm[nm!=""], recombArgs))
  if(length(bad)>0L){
    stop(paste0("Unused arguments: ", paste(bad, collapse=", "),
                ". Only ", paste(recombArgs, collapse=", "),
                " may be passed through ..."))
  }
  for(arg in nm){
    value = dots[[arg]]
    if(!is.numeric(value) | length(value)!=1L){
      stop(paste(arg, "must be a single number"))
    }
    if(!is.finite(value)){
      stop(paste(arg, "must be a single number"))
    }
  }
  if(!is.null(dots[["v"]])){
    if(dots[["v"]]<=0){
      stop("v must be greater than zero")
    }
  }
  if(!is.null(dots[["p"]])){
    if(dots[["p"]]<0 | dots[["p"]]>1){
      stop("p must be between zero and one")
    }
  }
  if(!is.null(dots[["quadProb"]])){
    if(dots[["quadProb"]]<0 | dots[["quadProb"]]>1){
      stop("quadProb must be between zero and one")
    }
  }
  return(dots)
}

# Build the temporary SimParam used when a crossing function is given a map
# population, which is taken to mean that no simulation has been set up yet
#
# mapPop, the MapPop or NamedMapPop defining the genetic map
# dots, the list of ... arguments, already passed through checkRecombArgs
# nThreads, the nThreads argument as the user gave it, possibly NULL
#
# Returns the SimParam. It is private to the calling function: it carries no
# traits, is never written to the global environment, and goes out of scope
# when the caller returns.
mapSimParam = function(mapPop, dots, nThreads){
  simParam = SimParam$new(mapPop)
  # Anything not supplied keeps SimParam's own default, so the defaults
  # are never written down twice
  if(!is.null(dots[["v"]])){
    simParam$v = dots[["v"]]
  }
  if(!is.null(dots[["p"]])){
    simParam$p = dots[["p"]]
  }
  if(!is.null(dots[["quadProb"]])){
    simParam$quadProb = dots[["quadProb"]]
  }
  if(!is.null(nThreads)){
    simParam$nThreads = as.integer(nThreads)
  }
  return(simParam)
}

# Coerce, check and extend a pedigree given to pedigreeCross
#
# id, mother, father, the pedigree as the user supplied it
# DH, nSelf, optional per individual vectors, either may be NULL
# unknownParent, the values in mother and father that mean the parent is
#   unknown, as well as NA
# matchID, whether the pedigree's names will be matched to founderPop
#
# All three pedigree vectors are coerced to character, and so are the
# unknownParent codes, so that a code of 0 matches a parent of "0". A parent
# that is NA or one of those codes comes back as NA_character_, and any other
# value names an individual. An id that is NA or one of those codes is not
# allowed, because it could never be told apart from an unknown parent.
#
# The pedigree is then extended backwards. A name used as a parent but
# without a row of its own is given one, with both of its parents unknown,
# and those rows are placed ahead of the supplied pedigree. With matchID
# every such name is added, because a name needs a row before it can be
# matched. Without it, only a name used twice or more, counting the mother
# and father vectors together, is added: a name used once carries no
# relationship to anything, so it is recoded as unknown instead.
#
# Returns the extended vectors, along with nAdded, the number of rows the
# extension put in front.
checkPedigreeInput = function(id, mother, father, DH, nSelf,
                              unknownParent=NA_character_, matchID=FALSE){
  id = as.character(id)
  mother = as.character(mother)
  father = as.character(father)
  unknownParent = as.character(unknownParent)
  unknownParent = unknownParent[!is.na(unknownParent)]
  if(is.null(DH)){
    DH = logical(length(id))
  }else{
    DH = as.logical(DH)
  }
  if(is.null(nSelf)){
    nSelf = rep(0, length(id))
  }

  # Check input data
  if(anyNA(id)){
    stop("id can not contain NA, because every individual needs a name")
  }
  if(any(id%in%unknownParent)){
    stop("id can not contain a value given in unknownParent: ",
         paste(unique(id[id%in%unknownParent]), collapse=", "))
  }
  if(any(duplicated(id))){
    stop("id contains duplicates")
  }
  if(length(id)!=length(mother)){
    stop("length(id) does not match length(mother)")
  }
  if(length(id)!=length(father)){
    stop("length(id) does not match length(father)")
  }
  if(length(id)!=length(DH)){
    stop("length(id) does not match length(DH)")
  }
  if(length(id)!=length(nSelf)){
    stop("length(id) does not match length(nSelf)")
  }
  if(length(id)==0L){
    stop("The pedigree is empty")
  }
  nSelf = suppressWarnings(as.integer(nSelf))
  if(anyNA(nSelf) | any(nSelf<0L)){
    stop("nSelf must be a non-negative integer for every individual")
  }
  if(anyNA(DH)){
    stop("DH must be TRUE or FALSE for every individual")
  }

  # Every way of writing an unknown parent becomes NA from here on
  mother[mother%in%unknownParent] = NA_character_
  father[father%in%unknownParent] = NA_character_

  # A parent with no row of its own. Uses are counted across both vectors,
  # so a name that is both parents of one individual counts twice and a
  # self is not turned into an outcross.
  parents = c(mother, father)
  missing = parents[!is.na(parents) & !(parents%in%id)]
  if(matchID){
    toAdd = unique(missing)
  }else{
    nUse = table(missing)
    toAdd = names(nUse)[nUse>=2L]
    dropped = names(nUse)[nUse<2L]
    mother[mother%in%dropped] = NA_character_
    father[father%in%dropped] = NA_character_
  }

  # Sorting the added names keeps the result the same however the supplied
  # pedigree happened to be ordered
  toAdd = sort(toAdd)
  nAdded = length(toAdd)
  if(nAdded>0L){
    id = c(toAdd, id)
    mother = c(rep(NA_character_, nAdded), mother)
    father = c(rep(NA_character_, nAdded), father)
    # An added individual is a founder, so it is neither selfed nor doubled
    DH = c(rep(FALSE, nAdded), DH)
    nSelf = c(rep(0L, nAdded), nSelf)
  }

  return(list(id=id, mother=mother, father=father, DH=DH, nSelf=nSelf,
              nAdded=nAdded))
}

# Work out how every individual of an extended pedigree is to be made
#
# id, mother, father, the extended pedigree, unknown parents being NA
# founderIds, the ids of the founder population, used only if matchID
# nFounderInd, the size of the founder population, used only if not matchID
# matchID, should the pedigree's names be matched to founderIds
#
# The pedigree has already been extended, so a parent is either another row
# of it or unknown.
#
# With matchID, an individual whose name is in founderIds takes its genotype
# from there and its own ancestry is not simulated. Everything descended from
# such an individual is simulated, and nothing else can be made at all, so an
# individual that is neither matched nor descended from a match is an error.
#
# Without matchID, the founders are the individuals whose parents are both
# unknown, and every individual with exactly one unknown parent needs a
# founder genome for that parent as well. All of them are drawn at random
# from the founder population.
#
# Returns build, saying which rows the result contains; founderRowFP, the
# position in the founder population to copy an individual from, or NA if it
# is to be crossed; motherPed and fatherPed, the rows its parents are; and
# motherFP and fatherFP, the positions to take a parent of unknown identity
# from. This takes no population and no SimParam so that it can be checked on
# its own.
resolveFounders = function(id, mother, father, founderIds, nFounderInd,
                           matchID){
  n = length(id)
  motherPed = match(mother, id)
  fatherPed = match(father, id)
  
  founderRowFP = rep(NA_integer_, n)
  motherFP = rep(NA_integer_, n)
  fatherFP = rep(NA_integer_, n)
  
  if(matchID){
    matched = id%in%founderIds
    if(!any(matched)){
      stop("matchID=TRUE, but no individual in the pedigree matches an ID in founderPop")
    }
    founderRowFP[matched] = match(id[matched], founderIds)
    
    # An individual can be made if it was matched, or if both of its parents
    # can be made. Ancestors of a matched individual are skipped, so they
    # need no genotype of their own.
    build = matched
    repeat{
      canCross = !build &
                 !is.na(motherPed) & !is.na(fatherPed) &
                 build[ifelse(is.na(motherPed), 1L, motherPed)] &
                 build[ifelse(is.na(fatherPed), 1L, fatherPed)]
      if(!any(canCross)){
        break
      }
      build = build | canCross
    }
    
    # Skipping an ancestor is fine, but an individual that is neither made
    # nor an ancestor of one cannot be accounted for at all
    isAncestor = rep(FALSE, n)
    repeat{
      wanted = build | isAncestor
      parents = unique(c(motherPed[wanted], fatherPed[wanted]))
      parents = parents[!is.na(parents)]
      newAncestor = isAncestor
      newAncestor[parents] = TRUE
      if(all(newAncestor==isAncestor)){
        break
      }
      isAncestor = newAncestor
    }
    orphan = !build & !isAncestor
    if(any(orphan)){
      stop(paste("Not enough individuals in founderPop match the pedigree.",
                 "These individuals are neither matched nor descended from a",
                 "match:", paste(id[orphan], collapse=", ")))
    }
  }else{
    # The founders are the individuals with no known parents, and a single
    # unknown parent needs a founder genome of its own
    build = rep(TRUE, n)
    isFounderRow = is.na(motherPed) & is.na(fatherPed)
    needMother = !isFounderRow & is.na(motherPed)
    needFather = !isFounderRow & is.na(fatherPed)
    
    nAnon = sum(isFounderRow) + sum(needMother) + sum(needFather)
    if(nAnon>nFounderInd){
      stop(paste("Pedigree requires",nAnon,"founders, but only",nFounderInd,"were supplied"))
    }
    
    # Randomly assign individuals as founders. The shuffle is deliberate: it
    # makes repeated gene drop replicates from an imported pedigree easy.
    pool = sample.int(nFounderInd, nAnon)
    
    # Hand out the founder genomes in a fixed order
    taken = 0L
    for(i in which(isFounderRow)){
      taken = taken + 1L
      founderRowFP[i] = pool[taken]
    }
    for(i in which(needMother)){
      taken = taken + 1L
      motherFP[i] = pool[taken]
    }
    for(i in which(needFather)){
      taken = taken + 1L
      fatherFP[i] = pool[taken]
    }
  }
  
  return(list(build=build,
              founderRowFP=founderRowFP,
              motherPed=motherPed,
              fatherPed=fatherPed,
              motherFP=motherFP,
              fatherFP=fatherFP))
}


#' @title Pedigree cross
#'
#' @description
#' Creates a \code{\link{Pop-class}} from a generic
#' pedigree and a set of founder individuals.
#'
#' @param founderPop a \code{\link{Pop-class}}, \code{\link{MapPop-class}}
#' or \code{\link{NamedMapPop-class}}. A map population is taken to mean
#' that no simulation has been set up yet: see details. Matching on ID needs
#' a population that has IDs, so a \code{\link{MapPop-class}} can only be
#' used with matchID=FALSE.
#' @param id a vector of unique identifiers for individuals
#' in the pedigree. The values of these IDs are separate from
#' the IDs in the founderPop if matchID=FALSE.
#' @param mother a vector of identifiers for the mothers
#' of individuals in the pedigree. See details for the treatment of
#' unknown parents.
#' @param father a vector of identifiers for the fathers
#' of individuals in the pedigree. See details for the treatment of
#' unknown parents.
#' @param matchID indicates if the IDs in founderPop should be
#' matched to the id argument. See details.
#' @param maxCycle the maximum number of passes to make over the pedigree
#' while working out which generation each individual belongs to. If
#' \code{NULL}, the number of individuals is used, which is the deepest a
#' pedigree of that size can be and so is always enough. A pedigree that
#' cannot be sorted is reported as such whatever this is set to, because a
#' pass that places nobody is a cycle and is detected as one.
#' @param DH an optional vector indicating if an individual
#' should be made a doubled haploid.
#' @param nSelf an optional vector indicating how many generations an
#' individual should be selfed.
#' @param useFemale If creating DH lines, should female recombination
#' rates be used. This parameter has no effect if recombRatio=1.
#' @param simParam an object of class \code{\link{SimParam}}. If
#' \code{NULL}, the function uses the object named \code{SP} from the
#' global environment. Ignored, with a warning, when founderPop is a map
#' population.
#' @param nThreads number of threads to use if OpenMP is available.
#' If \code{NULL}, the number is obtained from \code{simParam$nThreads}.
#' @param unknownParent the values in mother and father that mean the
#' parent is unknown, such as \code{"0"} or \code{""}. More than one may
#' be given. \code{NA} is always treated as unknown, whatever this is set
#' to, because it cannot name an individual.
#' @param ... the recombination settings \code{v}, \code{p} and
#' \code{quadProb}, used only when founderPop is a map population. Any that
#' are not given keep the \code{\link{SimParam}} defaults. Passing them
#' with a \code{\link{Pop-class}} is an error, because that simulation
#' already defines them.
#'
#' @return a \code{\link{Pop-class}}, or a
#' \code{\link{NamedMapPop-class}} when founderPop is a map population.
#' With matchID=TRUE this holds the matched individuals and their
#' descendants, which may be fewer than the pedigree supplied.
#'
#' @details
#' A map population and a population are handled differently, because they
#' mean different things about where the user is.
#'
#' Passing a \code{\link{Pop-class}} means a simulation is already running,
#' so the pedigree is built with the user's \code{\link{SimParam}} and the
#' result is a \code{\link{Pop-class}} belonging to that simulation.
#'
#' Passing a \code{\link{MapPop-class}} or a
#' \code{\link{NamedMapPop-class}} means no simulation has been set up yet.
#' No \code{\link{SimParam}} is looked for and none is created for the user:
#' the function makes a temporary one of its own, uses it to run the crossing,
#' and discards it. The result is a \code{\link{NamedMapPop-class}} carrying
#' the supplied genetic map and the pedigree's IDs, which the user can then
#' pass to \code{SimParam$new()} to start a simulation from. Because the
#' temporary \code{\link{SimParam}} defines no traits, nothing about traits
#' or SNP chips survives the call. The recombination settings it uses can be
#' given as \code{v}, \code{p} and \code{quadProb} through \code{...}.
#' The id, mother and father vectors are always taken as character vectors,
#' and anything else is coerced to one. A parent is unknown if it is
#' \code{NA} or one of the values given as unknownParent, and any other
#' value names an individual. By default \code{NA} is the only unknown
#' value, so a parent of \code{0} names an individual called "0" unless
#' \code{unknownParent="0"} is given. An id of \code{NA}, or an id that is
#' one of the unknownParent values, is not allowed.
#'
#' The pedigree is extended backwards before it is used, so that every
#' parent is either another row of the pedigree or unknown. An individual
#' named as a parent but without a row of its own is given one, with both of
#' its parents unknown, and those rows are placed ahead of the supplied
#' pedigree. How much is added depends on matchID.
#'
#' What happens next also depends on matchID.
#'
#' If matchID is FALSE, only the missing names that are needed are added. A
#' missing name is added when it is used as a parent twice or more, counting
#' the mother and father vectors together, because leaving it out would
#' break a relationship: two half sibs would come back unrelated, and an
#' individual whose mother and father are the same missing name would come
#' back as an outcross rather than a self. A missing name used once carries
#' no relationship, so it is treated as unknown, and the individual it
#' belonged to becomes a founder. A pedigree of id "3" with mother "2" and
#' father "1" is therefore the single founder "3", because neither "1" nor
#' "2" is used anywhere else. The returned mother and father show what was
#' simulated, so a dropped name comes back as \code{NA}.
#'
#' The founders are then the individuals whose parents are both unknown, and
#' they are randomly sampled from founderPop. An individual with exactly one
#' unknown parent is a half founder, and that parent is given a founder
#' genome of its own, so a pedigree with several half founders needs one
#' extra founder for each of them. Two half founders never share a genome.
#' Every individual in the extended pedigree is returned.
#'
#' If matchID is TRUE, every missing name is added, because a name has to
#' have a row of its own before it can be matched. An individual whose id is
#' found in founderPop takes its genotype from there, and its own ancestry
#' is not simulated. Only the
#' matched individuals and what descends from them are simulated, so the
#' returned population is smaller than the extended pedigree whenever the
#' pedigree reaches back beyond a match. An individual that is neither
#' matched nor descended from a match cannot be made at all, and the
#' function stops rather than guessing.
#'
#' @family mating functions
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
#' #Pedigree for a biparental cross with 7 generations of selfing.
#' #An unknown parent is NA.
#' id = as.character(1:10)
#' mother = c(NA, NA, "1", as.character(3:9))
#' father = c(NA, NA, "2", as.character(3:9))
#' pop2 = pedigreeCross(pop, id, mother, father, simParam=SP)
#'
#' #An incomplete pedigree, with half founders in rows 3 and 4
#' founderPop = quickHaplo(nInd=8, nChr=1, segSites=10)
#' SP = SimParam$new(founderPop)
#' \dontshow{SP$nThreads = 1L}
#' pop = newPop(founderPop, simParam=SP)
#' id = as.character(1:5)
#' mother = c(NA, NA, NA, "1", "1")
#' father = c(NA, NA, "2", NA, "2")
#' pop3 = pedigreeCross(pop, id, mother, father, simParam=SP)
#'
#' #A pedigree naming parents that have no row of their own. "A" is used
#' #twice, so it is added and "11" and "12" come back as half sibs. "B" and
#' #"C" are used once, so they are treated as unknown and "13" comes back as
#' #a founder.
#' id = c("11", "12", "13")
#' mother = c("A", "A", "B")
#' father = c(NA, NA, "C")
#' pop4 = pedigreeCross(pop, id, mother, father, simParam=SP)
#' getPed(pop4)
#'
#' #A pedigree that writes an unknown parent as 0
#' id = as.character(1:3)
#' mother = c("0", "0", "1")
#' father = c("0", "0", "2")
#' pop5 = pedigreeCross(pop, id, mother, father, unknownParent="0",
#'                      simParam=SP)
#'
#' @export
pedigreeCross = function(founderPop, id, mother, father, matchID=FALSE,
                         maxCycle=NULL, DH=NULL, nSelf=NULL, useFemale=TRUE,
                         simParam=NULL, nThreads=NULL,
                         unknownParent=NA_character_, ...){
  dots = checkRecombArgs(list(...))

  # A map population is what the user has before a simulation exists, so it
  # is taken to mean that there is no SimParam to find. The global
  # environment is not consulted and nothing is left behind for the user.
  mapInput = is(founderPop,"MapPop")
  
  if(mapInput){
    if(!is.null(simParam)){
      warning("simParam is ignored when founderPop is a map population, because a temporary SimParam is used")
    }
    if(matchID & !is(founderPop,"NamedMapPop")){
      stop("matchID=TRUE needs a population with IDs. Supply a NamedMapPop or a Pop, or use matchID=FALSE")
    }
    mapPop = founderPop

    # Private to this call; see mapSimParam
    simParam = mapSimParam(mapPop, dots=dots, nThreads=nThreads)
    nThreads = simParam$nThreads
    
    founderPop = newPop(mapPop, simParam=simParam, nThreads=nThreads)
  }else{
    # A SimParam is an R6 object, so setting these would change the user's
    # own simulation for good rather than just for this call
    if(length(dots)>0L){
      stop(paste("v, p and quadProb can only be set when founderPop is a map",
                 "population. Set them on the SimParam instead"))
    }
    
    if(is.null(simParam)){
      simParam = get("SP",envir=.GlobalEnv)
    }
    
    if(is.null(nThreads)){
      nThreads = simParam$nThreads
    }else{
      nThreads = as.integer(nThreads)
    }
  }
  
  if(simParam$sexes!="no"){
    stop("pedigreeCross currently only works with sex='no'")
  }
  
  if(!is(founderPop,"Pop")){
    stop("founderPop must be a Pop, a MapPop or a NamedMapPop")
  }
  
  # Coerce and check the pedigree, then extend it back to unknown parents
  input = checkPedigreeInput(id=id, mother=mother, father=father,
                             DH=DH, nSelf=nSelf,
                             unknownParent=unknownParent,
                             matchID=matchID)
  id = input$id
  mother = input$mother
  father = input$father
  DH = input$DH
  nSelf = input$nSelf
  
  # Give every individual a generation number, which also finds a cycle
  gen = sortPed(id=id, mother=mother, father=father,
                maxCycle=maxCycle)

  # How every individual is to be made
  plan = resolveFounders(id=id, mother=mother, father=father,
                         founderIds=founderPop@id,
                         nFounderInd=founderPop@nInd,
                         matchID=matchID)
  build = plan$build
  founderRowFP = plan$founderRowFP
  motherPed = plan$motherPed
  fatherPed = plan$fatherPed
  motherFP = plan$motherFP
  fatherFP = plan$fatherFP

  # Build the pedigree one generation at a time, making the whole of a
  # generation with a single call. A call into the crossing code costs
  # about the same whether it makes one individual or a thousand, so the
  # cost of a pedigree comes from its depth and not from its size.
  nGen = max(gen)
  genPop = vector("list", length=nGen)  # the individuals made in a generation
  genRow = vector("list", length=nGen)  # the pedigree rows they are
  posOf = integer(length(id))           # where a row sits in its generation

  for(g in seq_len(nGen)){
    idx = which(gen==g & build)
    if(length(idx)==0L){
      # With matchID a whole generation can sit above every match
      next
    }

    # Copies out of the founder population, and crosses, are made
    # separately and put back in the order the pedigree has them
    isCopy = !is.na(founderRowFP[idx])
    copyAt = which(isCopy)
    crossAt = which(!isCopy)
    parts = list()
    partRow = integer(0)

    if(length(copyAt)>0L){
      parts = c(parts, list(founderPop[founderRowFP[idx[copyAt]]]))
      partRow = c(partRow, copyAt)
    }

    if(length(crossAt)>0L){
      cr = idx[crossAt]
      mPed = motherPed[cr]
      fPed = fatherPed[cr]

      # Every parent is gathered once, however many of this generation's
      # individuals use it. Parents that are rows of the pedigree come from
      # the generations already made, one subset per generation, and the
      # rest come from the founder population.
      pedNeed = unique(c(mPed, fPed))
      pedNeed = pedNeed[!is.na(pedNeed)]
      poolList = list()
      poolRow = integer(0)
      for(sourceGen in sort(unique(gen[pedNeed]))){
        take = pedNeed[gen[pedNeed]==sourceGen]
        poolList = c(poolList, list(genPop[[sourceGen]][posOf[take]]))
        poolRow = c(poolRow, take)
      }
      nPed = length(poolRow)

      fpNeed = unique(c(motherFP[cr], fatherFP[cr]))
      fpNeed = fpNeed[!is.na(fpNeed)]
      if(length(fpNeed)>0L){
        poolList = c(poolList, list(founderPop[fpNeed]))
      }

      if(length(poolList)==1L){
        pool = poolList[[1]]
      }else{
        pool = mergePops(poolList)
      }

      # Both parents come out of the same pool, so the cross plan is two
      # columns of positions in it
      femaleRow = ifelse(is.na(mPed),
                         nPed + match(motherFP[cr], fpNeed),
                         match(mPed, poolRow))
      maleRow = ifelse(is.na(fPed),
                       nPed + match(fatherFP[cr], fpNeed),
                       match(fPed, poolRow))

      parts = c(parts, list(makeCross2(pool, pool,
                                       crossPlan=cbind(femaleRow, maleRow),
                                       simParam=simParam,
                                       nThreads=nThreads)))
      partRow = c(partRow, crossAt)
    }

    if(length(parts)==1L){
      genInd = parts[[1]]
    }else{
      genInd = mergePops(parts)
    }
    if(is.unsorted(partRow)){
      genInd = genInd[order(partRow)]
    }

    # Selfing, a round at a time across everyone that still has one to do.
    # Individuals needing fewer rounds are held back rather than left out,
    # so that a generation stays in the order its pedigree rows have.
    nSelfGen = nSelf[idx]
    for(selfRound in seq_len(max(nSelfGen))){
      take = which(nSelfGen>=selfRound)
      if(length(take)==length(idx)){
        genInd = self(genInd, simParam=simParam, nThreads=nThreads)
      }else{
        rest = which(nSelfGen<selfRound)
        genInd = mergePops(list(self(genInd, parents=take,
                                     simParam=simParam, nThreads=nThreads),
                                genInd[rest]))
        genInd = genInd[order(c(take, rest))]
      }
    }

    # Doubled haploids, all of the generation's in one call. A zero in nDH
    # is how an individual that is not a doubled haploid is passed over.
    dhGen = DH[idx]
    if(any(dhGen)){
      if(all(dhGen)){
        genInd = makeDH(genInd, useFemale=useFemale, simParam=simParam,
                        nThreads=nThreads)
      }else{
        made = which(dhGen)
        kept = which(!dhGen)
        genInd = mergePops(list(makeDH(genInd, nDH=as.integer(dhGen),
                                       useFemale=useFemale,
                                       simParam=simParam,
                                       nThreads=nThreads),
                                genInd[kept]))
        genInd = genInd[order(c(made, kept))]
      }
    }

    genPop[[g]] = genInd
    genRow[[g]] = idx
    posOf[idx] = seq_along(idx)
  }

  # Collapse to a population in pedigree order. With matchID the ancestors
  # of a matched individual are skipped, so the result holds the matched
  # individuals and what descends from them rather than every row of the
  # pedigree.
  builtRow = unlist(genRow)
  genPop = genPop[!vapply(genPop, is.null, logical(1))]
  if(length(genPop)==1L){
    output = genPop[[1]]
  }else{
    output = mergePops(genPop)
  }
  if(is.unsorted(builtRow)){
    output = output[order(builtRow)]
  }

  # Copy over names
  output@id = id[build]
  output@mother = mother[build]
  output@father = father[build]
  
  if(mapInput){
    # Give back what the user can start a simulation from: their genetic map,
    # the genotypes the pedigree produced, and the pedigree's own names. The
    # temporary SimParam does not escape with it.
    return(new("NamedMapPop",
               id=output@id,
               mother=output@mother,
               father=output@father,
               nInd=output@nInd,
               nChr=output@nChr,
               ploidy=output@ploidy,
               nLoci=output@nLoci,
               geno=output@geno,
               genMap=mapPop@genMap,
               centromere=mapPop@centromere,
               # Crossing breaks inbreeding, and the flag describes how the
               # haplotypes were built rather than what they now are
               inbred=FALSE))
  }
  
  return(output)
}
