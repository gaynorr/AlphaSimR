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
#' @param id optional ids to give to progeny
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @export
makeCross = function(pop,crossPlan,id=NULL,simParam){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  stopifnot(class(pop)=="Pop")
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],pop$id),
                      match(crossPlan[,2],pop$id))
    if(any(is.na(crossPlan))){
      stop("Failed to match supplied IDs")
    }
  }
  rawPop = new("RawPop",
               nInd=nrow(crossPlan),
               nChr=pop@nChr,
               ploidy=pop@ploidy,
               nLoci=pop@nLoci,
               geno=cross2(pop@geno,crossPlan[,1],
                           pop@geno,crossPlan[,2],
                           simParam@genMaps,
                           simParam@recombRatio))
  return(newPop(rawPop=rawPop,id=id,
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
#' @param id optional id to assign to progeny
#' @param balance if using gender, this option will balance the number 
#' of progeny per parent
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
randCross = function(pop,nCrosses,nProgeny=1,
                     id=NULL,balance=TRUE,simParam){
  if(simParam@gender=="no"){
    crossPlan = sampHalfDialComb(pop@nInd, nCrosses)
  }else{
    female = which(pop@gender=="F")
    if(length(female)==0){
      stop("population doesn't contain any females")
    }
    male = which(pop@gender=="M")
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
  return(makeCross(pop=pop,crossPlan=crossPlan,id=id,
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
#' @param id optional ids to give to progeny
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @export
makeCross2 = function(females,males,crossPlan,id=NULL,simParam){
  if(females@ploidy!=2){
    stop("Only works with diploids")
  }
  stopifnot(class(males)=="Pop",class(females)=="Pop")
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],females$id),
                      match(crossPlan[,2],males$id))
    if(any(is.na(crossPlan))){
      stop("Failed to match supplied IDs")
    }
  }
  rawPop = new("RawPop",
               nInd=nrow(crossPlan),
               nChr=females@nChr,
               ploidy=females@ploidy,
               nLoci=females@nLoci,
               geno=cross2(females@geno,crossPlan[,1],
                           males@geno,crossPlan[,2],
                           simParam@genMaps,
                           simParam@recombRatio))
  return(newPop(rawPop=rawPop,id=id,
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
#' @param id optional id to assign to progeny
#' @param balance if using gender, this option will balance the number 
#' of progeny per parent
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
randCross2 = function(females,males,nCrosses,nProgeny=1,
                     id=NULL,balance=TRUE,simParam){
  if(simParam@gender=="no"){
      crossPlan = sampAllComb(females@nInd,males@nInd,nCrosses)
  }else{
    female = which(females@gender=="F")
    if(length(female)==0){
      stop("population doesn't contain any females")
    }
    male = which(males@gender=="M")
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
                    id=id,simParam=simParam))
}

#' @title Self individuals
#' 
#' @description 
#' Creates selfed progeny from each individual in a 
#' population. Only works when gender is "no".
#' 
#' @param pop an object of \code{\link{Pop-class}}
#' @param nProgeny total number of selfed progeny per individual
#' @param id optional id to give to progeny
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
self = function(pop,nProgeny,id=NULL,simParam){
  stopifnot(simParam@gender=="no")
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  stopifnot(class(pop)=="Pop")
  crossPlan = cbind(rep(1:pop@nInd,each=nProgeny),
                    rep(1:pop@nInd,each=nProgeny))
  rawPop = new("RawPop",
               nInd=nrow(crossPlan),
               nChr=pop@nChr,
               ploidy=pop@ploidy,
               nLoci=pop@nLoci,
               geno=cross2(pop@geno,crossPlan[,1],
                           pop@geno,crossPlan[,2],
                           simParam@genMaps,
                           simParam@recombRatio))
  return(newPop(rawPop=rawPop,id=id,
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
#' @param id optional id to assign to DH lines
#' @param useFemale should female recombination rates be used. 
#' This parameter has no effect if, recombRatio=1.
#' @param simParam an object of 'SimParam' class
#' 
#' @export
makeDH = function(pop,nDH,id=NULL,useFemale=TRUE,simParam){
  stopifnot(simParam@gender=="no")
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  stopifnot(class(pop)=="Pop")
  rawPop = new("RawPop",
               nInd=as.integer(pop@nInd*nDH),
               nChr=pop@nChr,
               ploidy=pop@ploidy,
               nLoci=pop@nLoci,
               geno=createDH2(pop@geno,nDH,
                              simParam@genMaps,
                              simParam@recombRatio,
                              useFemale))
  return(newPop(rawPop=rawPop,id=id,
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
#' @param id optional id to assign to progeny
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
pedigreeCross = function(pedigree,founders,id=NULL,simParam){
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
                                  simParam@genMaps))
  return(newPop(rawPop=rawPop,id=id,
                mother=as.character(sortedped@mother),
                father=as.character(sortedped@father),
                simParam=simParam))
}
