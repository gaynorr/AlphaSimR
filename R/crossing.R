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
#' @param rawPop should a \code{\link{RawPop-class}} be returned
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{RawPop-class}}
#'
#' @export
makeCross = function(pop,crossPlan,rawPop=FALSE,
                     simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],pop$id),
                      match(crossPlan[,2],pop$id))
    if(any(is.na(crossPlan))){
      stop("Failed to match supplied IDs")
    }
  }
  rPop = new("RawPop",
             nInd=nrow(crossPlan),
             nChr=pop@nChr,
             ploidy=pop@ploidy,
             nLoci=pop@nLoci,
             geno=cross2(pop@geno,crossPlan[,1],
                         pop@geno,crossPlan[,2],
                         simParam@genMaps,
                         simParam@recombRatio))
  if(rawPop) return(rPop)
  return(newPop(rawPop=rPop,
                mother=pop@id[crossPlan[,1]],
                father=pop@id[crossPlan[,2]],
                simParam=simParam))
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
#' @param balance if using gender, this option will balance the number 
#' of progeny per parent
#' @param parents an optional vector of indices for allowable parents
#' @param rawPop should a \code{\link{RawPop-class}} be returned
#' @param ignoreGender should gender be ignored
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{RawPop-class}}
#' 
#' @export
randCross = function(pop,nCrosses,nProgeny=1,
                     balance=TRUE,parents=NULL,
                     rawPop=FALSE,ignoreGender=FALSE,
                     simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  if(is.null(parents)){
    parents = 1:pop@nInd
  }else{
    parents = as.integer(parents)
  }
  n = length(parents)
  if(simParam@gender=="no" | ignoreGender){
    crossPlan = sampHalfDialComb(n, nCrosses)
    crossPlan[,1] = parents[crossPlan[,1]]
    crossPlan[,2] = parents[crossPlan[,2]]
  }else{
    female = which(pop@gender=="F" & (1:pop@nInd)%in%parents)
    if(length(female)==0){
      stop("population doesn't contain any females")
    }
    male = which(pop@gender=="M" & (1:pop@nInd)%in%parents)
    if(length(male)==0){
      stop("population doesn't contain any males")
    }
    if(balance){
      female = female[sample.int(length(female),length(female))]
      female = rep(female,length.out=nCrosses)
      male = male[sample.int(length(male),length(male))]
      male = rep(male,length.out=nCrosses)
      male = male[sample.int(nCrosses,nCrosses)]
      crossPlan = cbind(female,male)
    }else{
      crossPlan = sampAllComb(length(female),
                              length(male),
                              nCrosses)
      crossPlan[,1] = female[crossPlan[,1]]
      crossPlan[,2] = male[crossPlan[,2]]
    }
  }
  if(nProgeny>1){
    crossPlan = cbind(rep(crossPlan[,1],each=nProgeny),
                      rep(crossPlan[,2],each=nProgeny))
  }
  return(makeCross(pop=pop,crossPlan=crossPlan,
                   rawPop=rawPop,simParam=simParam))
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
#' are selected without regards to gender and it supercedes values 
#' for nFemale and nMale. Thus if the simulation uses gender, it is 
#' likely better to leave this value as NULL and use nFemale and nMale 
#' instead.
#' @param nFemale the number of females to select. This value is ignored 
#' if nInd is set.
#' @param nMale the number of males to select. This value is ignored 
#' if nInd is set.
#' @param trait the trait for selection. Either a number indicating 
#' a single trait or a function returning a vector of length nInd.
#' @param use select on genetic values "gv", estimated
#' breeding values "ebv", breeding values "bv", phenotypes "pheno", 
#' or randomly "rand"
#' @param selectTop selects highest values if true. 
#' Selects lowest values if false.
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross
#' @param balance if using gender, this option will balance the number 
#' of progeny per parent
#' @param simParam an object of \code{\link{SimParam-class}}
#' @param ... additional arguments if using a function for 
#' trait
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
selectCross = function(pop,nInd=NULL,nFemale=NULL,nMale=NULL,trait=1,
                       selectTop=TRUE,use="pheno",nCrosses,nProgeny=1,
                       balance=TRUE,simParam=NULL,...){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  if(!is.null(nInd)){
    parents = selectInd(pop=pop,nInd=nInd,trait=trait,use=use,
                        gender="B",selectTop=selectTop,
                        returnPop=FALSE,simParam=simParam,...)
  }else{
    if(simParam@gender=="no")
      stop("You must specify nInd when simParam@gender is `no`")
    if(is.null(nFemale))
      stop("You must specify nFemale if nInd is NULL")
    if(is.null(nMale))
      stop("You must specify nMale if nInd is NULL")
    females = selectInd(pop=pop,nInd=nFemale,trait=trait,use=use,
                        gender="F",selectTop=selectTop,
                        returnPop=FALSE,simParam=simParam,...)
    males = selectInd(pop=pop,nInd=nMale,trait=trait,use=use,
                      gender="M",selectTop=selectTop,
                      returnPop=FALSE,simParam=simParam,...)
    parents = c(females,males)
  }
  
  return(randCross(pop=pop,nCrosses=nCrosses,nProgeny=nProgeny,
                   balance=balance,parents=parents,
                   rawPop=FALSE,ignoreGender=FALSE,
                   simParam=simParam))
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
#' @param rawPop should a \code{\link{RawPop-class}} be returned
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{RawPop-class}}
#'
#' @export
makeCross2 = function(females,males,crossPlan,
                      rawPop=FALSE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  if(females@ploidy!=2){
    stop("Only works with diploids")
  }
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],females$id),
                      match(crossPlan[,2],males$id))
    if(any(is.na(crossPlan))){
      stop("Failed to match supplied IDs")
    }
  }
  rPop = new("RawPop",
             nInd=nrow(crossPlan),
             nChr=females@nChr,
             ploidy=females@ploidy,
             nLoci=females@nLoci,
             geno=cross2(females@geno,crossPlan[,1],
                         males@geno,crossPlan[,2],
                         simParam@genMaps,
                         simParam@recombRatio))
  if(rawPop) return(rPop)
  return(newPop(rawPop=rPop,
                mother=females@id[crossPlan[,1]],
                father=males@id[crossPlan[,2]],
                simParam=simParam))
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
#' @param balance if using gender, this option will balance the number 
#' of progeny per parent
#' @param femaleParents an optional vector of indices for allowable 
#' female parents
#' @param maleParents an optional vector of indices for allowable 
#' male parents
#' @param rawPop should a \code{\link{RawPop-class}} be returned
#' @param ignoreGender should gender be ignored
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{RawPop-class}}
#' 
#' @export
randCross2 = function(females,males,nCrosses,nProgeny=1,
                      balance=TRUE,femaleParents=NULL,
                      maleParents=NULL,rawPop=FALSE,
                      ignoreGender=FALSE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
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
  nF = length(femaleParents)
  nM = length(maleParents)
  if(simParam@gender=="no" | ignoreGender){
      crossPlan = sampAllComb(nF,nM,nCrosses)
      crossPlan[,1] = femaleParents[crossPlan[,1]]
      crossPlan[,2] = maleParents[crossPlan[,2]]
  }else{
    female = which(females@gender=="F" & 
                     (1:females@nInd)%in%femaleParents)
    if(length(female)==0){
      stop("population doesn't contain any females")
    }
    male = which(males@gender=="M" & 
                   (1:males@nInd)%in%maleParents)
    if(length(male)==0){
      stop("population doesn't contain any males")
    }
    if(balance){
      female = female[sample.int(length(female),length(female))]
      female = rep(female,length.out=nCrosses)
      male = male[sample.int(length(male),length(male))]
      male = rep(male,length.out=nCrosses)
      male = male[sample.int(nCrosses,nCrosses)]
      crossPlan = cbind(female,male)
    }else{
      crossPlan = sampAllComb(length(female),
                              length(male),
                              nCrosses)
      crossPlan[,1] = female[crossPlan[,1]]
      crossPlan[,2] = male[crossPlan[,2]]
    }
  }
  if(nProgeny>1){
    crossPlan = cbind(rep(crossPlan[,1],each=nProgeny),
                      rep(crossPlan[,2],each=nProgeny))
  }
  return(makeCross2(females=females,males=males,
                    crossPlan=crossPlan,
                    rawPop=rawPop,simParam=simParam))
}

#' @title Self individuals
#' 
#' @description 
#' Creates selfed progeny from each individual in a 
#' population. Only works when gender is "no".
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param nProgeny total number of selfed progeny per individual
#' @param parents an optional vector of indices for allowable parents
#' @param rawPop should a \code{\link{RawPop-class}} be returned
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{RawPop-class}}
#' 
#' @export
self = function(pop,nProgeny=1,parents=NULL,
                rawPop=FALSE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  if(is.null(parents)){
    parents = 1:pop@nInd
  }else{
    parents = as.integer(parents)
  }
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  crossPlan = rep(parents,each=nProgeny)
  crossPlan = cbind(crossPlan,crossPlan)
  rPop = new("RawPop",
             nInd=nrow(crossPlan),
             nChr=pop@nChr,
             ploidy=pop@ploidy,
             nLoci=pop@nLoci,
             geno=cross2(pop@geno,crossPlan[,1],
                         pop@geno,crossPlan[,2],
                         simParam@genMaps,
                         simParam@recombRatio))
  if(rawPop) return(rPop)
  return(newPop(rawPop=rPop,
                mother=rep(pop@mother,each=nProgeny),
                father=rep(pop@father,each=nProgeny),
                simParam=simParam))
}

#' @title Generates DH lines
#' 
#' @description Creates DH lines from each individual in a population. 
#' Only works when gender is "no".
#' 
#' @param pop an object of 'Pop' superclass
#' @param nDH total number of DH lines per individual
#' @param useFemale should female recombination rates be used. 
#' This parameter has no effect if, recombRatio=1.
#' @param rawPop should a \code{\link{RawPop-class}} be returned
#' @param simParam an object of 'SimParam' class
#' 
#' @return Returns an object of \code{\link{Pop-class}} or 
#' \code{\link{RawPop-class}}
#' 
#' @export
makeDH = function(pop,nDH=1,useFemale=TRUE,
                  rawPop=FALSE,simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  rPop = new("RawPop",
             nInd=as.integer(pop@nInd*nDH),
             nChr=pop@nChr,
             ploidy=pop@ploidy,
             nLoci=pop@nLoci,
             geno=createDH2(pop@geno,nDH,
                            simParam@genMaps,
                            simParam@recombRatio,
                            useFemale))
  if(rawPop) return(rPop)
  return(newPop(rawPop=rPop,
                mother=rep(pop@mother,each=nDH),
                father=rep(pop@father,each=nDH),
                simParam=simParam))
}

#' @title Make crosses based on a pedigree
#' 
#' @description 
#' To Do
#' 
#' @param pedigree an object of \code{\link{Pedigree-class}}
#' @param founders an object of \code{\link{Pop-class}}
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
pedigreeCross = function(pedigree,founders,
                         simParam=NULL){
  if(is.null(simParam)){
    simParam = get("SIMPARAM",envir=.GlobalEnv)
  }
  stopifnot(class(pedigree)=="Pedigree")
  stopifnot(class(founders)=="MapPop")
  
  sortedped = sortPed(pedigree)
  rawPop = new("RawPop",
               nInd=sortedped@nInd,
               nChr=founders@nChr,
               ploidy=founders@ploidy,
               nLoci=founders@nLoci,
               geno=crossPedigree(founders@geno,
                                  sortedped@mother,
                                  sortedped@father,
                                  simParam@genMaps,
                                  simParam@recombRatio))
  return(newPop(rawPop=rawPop,
                mother=as.character(sortedped@mother),
                father=as.character(sortedped@father),
                simParam=simParam))
}
