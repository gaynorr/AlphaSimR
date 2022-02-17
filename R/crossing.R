#' @title Make designed crosses
#' 
#' @description
#' Makes crosses within a population using a user supplied 
#' crossing plan.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param crossPlan a matrix with two column representing 
#' female and male parents. Either integers for the position in 
#' population or character strings for the IDs.
#' @param nProgeny number of progeny per cross
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
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Cross individual 1 with individual 10
#' crossPlan = matrix(c(1,10), nrow=1, ncol=2)
#' pop2 = makeCross(pop, crossPlan, simParam=SP)
#'
#' @export
makeCross = function(pop,crossPlan,nProgeny=1,
                     simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(pop@ploidy%%2L != 0L){
    stop("You can not cross indiviuals with odd ploidy levels")
  }
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],pop@id),
                      match(crossPlan[,2],pop@id))
    if(any(is.na(crossPlan))){
      stop("Failed to match supplied IDs")
    }
  }
  if((max(crossPlan)>nInd(pop)) |
     (min(crossPlan)<1L)){
    stop("Invalid crossPlan")
  }
  if(nProgeny>1){
    crossPlan = cbind(rep(crossPlan[,1],each=nProgeny),
                      rep(crossPlan[,2],each=nProgeny))
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
              simParam$nThreads)
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
  return(newPop(rawPop=rPop,
                mother=pop@id[crossPlan[,1]],
                father=pop@id[crossPlan[,2]],
                simParam=simParam,
                iMother=pop@iid[crossPlan[,1]],
                iFather=pop@iid[crossPlan[,2]],
                femaleParentPop=pop,
                maleParentPop=pop,
                hist=hist))
}

#' @title Make random crosses
#' 
#' @description 
#' A wrapper for \code{\link{makeCross}} that randomly 
#' selects parental combinations for all possible combinantions. 
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross
#' @param balance if using sexes, this option will balance the number 
#' of progeny per parent
#' @param parents an optional vector of indices for allowable parents
#' @param ignoreSexes should sexes be ignored
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
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Make 10 crosses
#' pop2 = randCross(pop, 10, simParam=SP)
#' 
#' @export
randCross = function(pop,nCrosses,nProgeny=1,
                     balance=TRUE,parents=NULL,
                     ignoreSexes=FALSE,
                     simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
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
      for(i in 1:n){
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
  return(makeCross(pop=pop,crossPlan=crossPlan,nProgeny=nProgeny,simParam=simParam))
}

#' @title Select and randomly cross
#' 
#' @description 
#' This is a wrapper that combines the functionalities of 
#' \code{\link{randCross}} and \code{\link{selectInd}}. The 
#' purpose of this wrapper is to combine both selection and 
#' crossing in one function call that minimized the amount 
#' of intermediate populations created. This reduces RAM usage 
#' and simplifies code writing. Note that this wrapper does not 
#' provide the full functionality of either function. 
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param nInd the number of individuals to select. These individuals 
#' are selected without regards to sex and it supercedes values 
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
#' @param simParam an object of \code{\link{SimParam}}
#' @param ... additional arguments if using a function for 
#' trait
#' @param balance if using sexes, this option will balance the number 
#' of progeny per parent. This argument occurs after ..., so the argument 
#' name must be matched exactly.
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
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Select 4 individuals and make 8 crosses
#' pop2 = selectCross(pop, nInd=4, nCrosses=8, simParam=SP)
#' 
#' @export
selectCross = function(pop,nInd=NULL,nFemale=NULL,nMale=NULL,nCrosses,
                       nProgeny=1,trait=1,use="pheno",selectTop=TRUE,
                       simParam=NULL,...,balance=TRUE){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(!is.null(nInd)){
    parents = selectInd(pop=pop,nInd=nInd,trait=trait,use=use,
                        sex="B",selectTop=selectTop,
                        returnPop=FALSE,simParam=simParam,...)
  }else{
    if(simParam$sexes=="no")
      stop("You must specify nInd when simParam$sexes is `no`")
    if(is.null(nFemale))
      stop("You must specify nFemale if nInd is NULL")
    if(is.null(nMale))
      stop("You must specify nMale if nInd is NULL")
    females = selectInd(pop=pop,nInd=nFemale,trait=trait,use=use,
                        sex="F",selectTop=selectTop,
                        returnPop=FALSE,simParam=simParam,...)
    males = selectInd(pop=pop,nInd=nMale,trait=trait,use=use,
                      sex="M",selectTop=selectTop,
                      returnPop=FALSE,simParam=simParam,...)
    parents = c(females,males)
  }
  
  return(randCross(pop=pop,nCrosses=nCrosses,nProgeny=nProgeny,
                   balance=balance,parents=parents,
                   ignoreSexes=FALSE,simParam=simParam))
}

#' @title Make designed crosses
#' 
#' @description
#' Makes crosses between two populations using a user supplied 
#' crossing plan.
#'
#' @param females an object of \code{\link{Pop-class}} for female parents.
#' @param males an object of \code{\link{Pop-class}} for male parents.
#' @param crossPlan a matrix with two column representing 
#' female and male parents. Either integers for the position in 
#' population or character strings for the IDs.
#' @param nProgeny number of progeny per cross
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
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Cross individual 1 with individual 10
#' crossPlan = matrix(c(1,10), nrow=1, ncol=2)
#' pop2 = makeCross2(pop, pop, crossPlan, simParam=SP)
#' 
#' @export
makeCross2 = function(females,males,crossPlan,nProgeny=1,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
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
  if(nProgeny>1){
    crossPlan = cbind(rep(crossPlan[,1],each=nProgeny),
                      rep(crossPlan[,2],each=nProgeny))
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
            simParam$nThreads)
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
  return(newPop(rawPop=rPop,
                mother=females@id[crossPlan[,1]],
                father=males@id[crossPlan[,2]],
                simParam=simParam,
                iMother=females@iid[crossPlan[,1]],
                iFather=males@iid[crossPlan[,2]],
                femaleParentPop=females,
                maleParentPop=males,
                hist=hist))
}

#' @title Make random crosses
#' 
#' @description 
#' A wrapper for \code{\link{makeCross2}} that randomly 
#' selects parental combinations for all possible combinantions between 
#' two populations.
#' 
#' @param females an object of \code{\link{Pop-class}} for female parents.
#' @param males an object of \code{\link{Pop-class}} for male parents.
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross
#' @param balance this option will balance the number 
#' of progeny per parent
#' @param femaleParents an optional vector of indices for allowable 
#' female parents
#' @param maleParents an optional vector of indices for allowable 
#' male parents
#' @param ignoreSexes should sex be ignored
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
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Make 10 crosses
#' pop2 = randCross2(pop, pop, 10, simParam=SP)
#' 
#' @export
randCross2 = function(females,males,nCrosses,nProgeny=1,
                      balance=TRUE,femaleParents=NULL,
                      maleParents=NULL,ignoreSexes=FALSE,
                      simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
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
  nMale = length(male)
  nFemale = length(female)
  if(balance){
    female = female[sample.int(nFemale, nFemale)]
    female = rep(female, length.out=nCrosses)
    tmp = male[sample.int(nMale, nMale)]
    n = nCrosses%/%nMale + 1
    male = NULL
    for(i in 1:n){
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
  return(makeCross2(females=females,males=males,
                    crossPlan=crossPlan,nProgeny=nProgeny,
                    simParam=simParam))
}

#' @title Self individuals
#' 
#' @description 
#' Creates selfed progeny from each individual in a 
#' population. Only works when sexes is "no".
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param nProgeny total number of selfed progeny per individual
#' @param parents an optional vector of indices for allowable parents
#' @param keepParents should previous parents be used for mother and 
#' father. 
#' @param simParam an object of \code{\link{SimParam}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Self pollinate each individual
#' pop2 = self(pop, simParam=SP)
#' 
#' @export
self = function(pop,nProgeny=1,parents=NULL,keepParents=TRUE,
                simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(class(pop)=="MegaPop"){
    stopifnot(is.null(parents))
    pop@pops = lapply(pop@pops, self, nProgeny=nProgeny, 
                      parents=NULL, keepParents=keepParents, 
                      simParam=simParam)
    return(pop)
  }
  if(is.null(parents)){
    parents = 1:pop@nInd
  }else{
    parents = as.integer(parents)
  }
  if(pop@ploidy%%2L != 0L){
    stop("You can not self aneuploids")
  }
  crossPlan = rep(parents,each=nProgeny)
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
              simParam$nThreads)
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
    return(newPop(rawPop=rPop,
                  mother=rep(pop@mother,each=nProgeny),
                  father=rep(pop@father,each=nProgeny),
                  simParam=simParam,
                  iMother=rep(pop@iid,each=nProgeny),
                  iFather=rep(pop@iid,each=nProgeny),
                  femaleParentPop=pop,
                  maleParentPop=pop,
                  hist=hist))
  }else{
    return(newPop(rawPop=rPop,
                  mother=rep(pop@id,each=nProgeny),
                  father=rep(pop@id,each=nProgeny),
                  simParam=simParam,
                  iMother=rep(pop@iid,each=nProgeny),
                  iFather=rep(pop@iid,each=nProgeny),
                  femaleParentPop=pop,
                  maleParentPop=pop,
                  hist=hist))
  }
}

#' @title Generates DH lines
#' 
#' @description Creates DH lines from each individual in a population. 
#' Only works with diploid individuals. For polyploids, use 
#' \code{\link{reduceGenome}} and \code{\link{doubleGenome}}.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nDH total number of DH lines per individual
#' @param useFemale should female recombination rates be used. 
#' @param keepParents should previous parents be used for mother and 
#' father. 
#' @param simParam an object of 'SimParam' class
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Create 1 DH for each individual
#' pop2 = makeDH(pop, simParam=SP)
#' 
#' @export
makeDH = function(pop,nDH=1,useFemale=TRUE,keepParents=TRUE,
                  simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  if(class(pop)=="MegaPop"){
    pop@pops = lapply(pop@pops, makeDH, nDH=nDH, useFemale=useFemale,
                      keepParents=keepParents, simParam=simParam)
    return(pop)
  }
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  if(useFemale){
    tmp = createDH2(pop@geno,nDH,
                    simParam$femaleMap,
                    simParam$v,
                    simParam$p,
                    simParam$isTrackRec,
                    simParam$nThreads)
  }else{
    tmp = createDH2(pop@geno,nDH,
                    simParam$maleMap,
                    simParam$v,
                    simParam$p,
                    simParam$isTrackRec,
                    simParam$nThreads)
  }
  dim(tmp$geno) = NULL # Account for matrix bug in RcppArmadillo
  rPop = new("RawPop",
             nInd=as.integer(pop@nInd*nDH),
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
    return(newPop(rawPop=rPop,
                  mother=rep(pop@mother, each=nDH),
                  father=rep(pop@father, each=nDH),
                  simParam=simParam,
                  isDH=TRUE,
                  iMother=rep(pop@iid, each=nDH),
                  iFather=rep(pop@iid, each=nDH),
                  femaleParentPop=pop,
                  maleParentPop=pop,
                  hist=hist))
  }else{
    return(newPop(rawPop=rPop,
                  mother=rep(pop@id, each=nDH),
                  father=rep(pop@id, each=nDH),
                  simParam=simParam,
                  isDH=TRUE,
                  iMother=rep(pop@iid, each=nDH),
                  iFather=rep(pop@iid, each=nDH),
                  femaleParentPop=pop,
                  maleParentPop=pop,
                  hist=hist))
  }
}


# Sort Pedigree
#
# id, id of individual
# mother, name of individual's mother
# father, name of individual's father
# maxCycle, number of loops for attempting to sort the pedigree
sortPed = function(id, mother, father, maxCycle=100){
  nInd = length(id)
  output = data.frame(gen=integer(nInd),
                      id=as.character(id),
                      mother=match(mother, id),
                      father=match(father, id),
                      motherID=as.character(mother),
                      fatherID=as.character(father))
  unsorted = rep(TRUE, nInd)
  for(gen in 1:maxCycle){
    for(i in which(unsorted)){
      if(is.na(output$mother[i])&is.na(output$father[i])){
        # Is a founder
        output$gen[i] = gen
        unsorted[i] = FALSE
      }else if(is.na(output$mother[i])){
        # Mother is a founder
        if(!unsorted[output$father[i]]){
          output$gen[i] = gen
          unsorted[i] = FALSE
        }
      }else if(is.na(output$father[i])){
        # Father is a founder
        if(!unsorted[output$mother[i]]){
          output$gen[i] = gen
          unsorted[i] = FALSE
        }
      }else{
        # Both parents are in the pedigree
        if(!unsorted[output$mother[i]] & !unsorted[output$father[i]]){
          output$gen[i] = gen
          unsorted[i] = FALSE
        }
      }
    }
  }
  if(any(unsorted)){
    stop("Failed to sort pedigree, may contain loops or require a higher maxGen")
  }
  return(output)
}

#' @title Pedigree cross
#' 
#' @description
#' Creates a \code{\link{Pop-class}} from a generic 
#' pedigree and a set of founder individuals. 
#'
#' @param founderPop a \code{\link{Pop-class}}
#' @param id a vector of unique identifiers for individuals 
#' in the pedigree. The values of these IDs are seperate from   
#' the IDs in the founderPop if matchID=FALSE.
#' @param mother a vector of identifiers for the mothers 
#' of individuals in the pedigree. Must match one of the 
#' elements in the id vector or they will be treated as unknown.
#' @param father a vector of identifiers for the fathers 
#' of individuals in the pedigree. Must match one of the 
#' elements in the id vector or they will be treated as unknown.
#' @param matchID indicates if the IDs in founderPop should be 
#' matched to the id argument. See details.
#' @param maxCycle the maximum number of loops to make over the pedigree 
#' to sort it.
#' @param DH an optional vector indicating if an individual 
#' should be made a doubled haploid.
#' @param useFemale If creating DH lines, should female recombination 
#' rates be used. This parameter has no effect if, recombRatio=1.
#' @param simParam an object of 'SimParam' class
#' 
#' @description 
#' The way in which the user supplied pedigree is used depends on 
#' the value of matchID. If matchID is TRUE, the IDs in the user 
#' supplied pedigree are matched against founderNames. If matchID 
#' is FALSE, founder individuals in the user supplied pedigree are 
#' randomly sampled from founderPop.
#' 
#' @examples 
#' #Create founder haplotypes
#' founderPop = quickHaplo(nInd=2, nChr=1, segSites=10)
#' 
#' #Set simulation parameters
#' SP = SimParam$new(founderPop)
#' 
#' #Create population
#' pop = newPop(founderPop, simParam=SP)
#' 
#' #Pedigree for a biparental cross with 7 generations of selfing
#' id = 1:10
#' mother = c(0,0,1,3:9)
#' father = c(0,0,2,3:9)
#' pop2 = pedigreeCross(pop, id, mother, father, simParam=SP)
#' 
#' 
#' @export
pedigreeCross = function(founderPop, id, mother, father, matchID=FALSE, 
                         maxCycle=100, DH=NULL, useFemale=TRUE, simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SP",envir=.GlobalEnv)
  }
  
  if(simParam$sexes!="no"){
    stop("pedigreeCross currently only works with sex='no'")
  }
  
  # Coerce input data
  id = as.character(id)
  mother = as.character(mother)
  father = as.character(father)
  if(is.null(DH)){
    DH = logical(length(id))
  }else{
    DH = as.logical(DH)
  }
  
  # Check input data
  stopifnot(!any(duplicated(id)),
            length(id)==length(mother),
            length(id)==length(father),
            length(id)==length(DH))
  
  # Sort pedigree (identifies potential problems)
  ped = sortPed(id=id, mother=mother, father=father, 
                maxCycle=maxCycle)
  
  # Create list for new population
  output = vector("list", length=length(id))
  
  # Order and assign founders 
  isFounder = is.na(ped$father) & is.na(ped$mother)
  motherIsFounder = is.na(ped$mother) & !is.na(ped$father)
  fatherIsFounder = is.na(ped$father) & !is.na(ped$mother)
  founderNames = c(unique(id[isFounder]),
                   unique(mother[motherIsFounder]),
                   unique(father[fatherIsFounder]))
  nFounder = length(founderNames)
  if(matchID){
    # Check that all founders are present
    founderPresent = founderNames%in%founderPop@id
    if(!all(founderPresent)){
      stop(paste("The following founders are missing:", founderNames[!founderPresent]))
    }
  }else{
    # Check that there are enough founders
    if(nFounder>founderPop@nInd){
      stop(paste("Pedigree requires",nFounder,"founders, but only",founderPop@nInd,"were supplied"))
    }
    
    # Randomly assign individuals as founders
    founderPop = founderPop[sample.int(founderPop@nInd,nFounder)]
    
    # isFounder
    n1 = 1
    n2 = sum(isFounder)
    founderPop@id[n1:n2] = id[isFounder]
    founderPop@mother[n1:n2] = mother[isFounder]
    founderPop@father[n1:n2] = father[isFounder]
    
    # motherIsFounder
    n = sum(motherIsFounder)
    if(n>=1){
      n1 = n2 + 1
      n2 = n2 + n
      founderPop@id[n1:n2] = mother[motherIsFounder]
      founderPop@mother[n1:n2] = rep("0", n2-n1+1)
      founderPop@father[n1:n2] = rep("0", n2-n1+1)
    }
    
    # fatherIsFounder
    n = sum(fatherIsFounder)
    if(n>=1){
      n1 = n2 + 1
      n2 = n2 + n
      founderPop@id[n1:n2] = father[fatherIsFounder]
      founderPop@mother[n1:n2] = rep("0", n2-n1+1)
      founderPop@father[n1:n2] = rep("0", n2-n1+1)
    }
  }
  
  # Create individuals
  crossPlan = matrix(c(1,1),ncol=2)
  for(gen in 1:max(ped$gen)){
    for(i in which(ped$gen==gen)){
      if(isFounder[i]){
        # Copy over founder individual
        output[[i]] = founderPop[id[i]]
      }else{
        if(motherIsFounder[i]){
          # Cross founder to newly created individual
          output[[i]] = makeCross2(founderPop[id[i]],
                                   output[[ped$father[i]]],
                                   crossPlan=crossPlan,
                                   simParam=simParam)
        }else if(fatherIsFounder[i]){
          # Cross newly created individual to founder
          output[[i]] = makeCross2(output[[ped$mother[i]]],
                                   founderPop[id[i]],
                                   crossPlan=crossPlan,
                                   simParam=simParam)
        }else{
          # Cross two newly created individuals
          output[[i]] = makeCross2(output[[ped$mother[i]]],
                                   output[[ped$father[i]]],
                                   crossPlan=crossPlan,
                                   simParam=simParam)
        }
        
        # Make the individual a DH?
        if(DH[i]){
          output[[i]] = makeDH(output[[i]],
                               useFemale=useFemale,
                               simParam=simParam)
        }
      } 
    }
  }
  
  # Collapse list to a population
  output = mergePops(output)
  
  # Copy over names
  output@id = id
  output@mother = mother
  output@father = father
  
  return(output)
}

