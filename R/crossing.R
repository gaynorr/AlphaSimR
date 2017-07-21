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
makeCross = function(pop,crossPlan,id=NULL,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  stopifnot(class(pop)=="Pop")
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],pop$id),
                      match(crossPlan[,2],pop$id))
    if(any(is.na(crossPlan))){
      stop("Failed to matched supplied IDs")
    }
  }
  geno = cross2(pop@geno,crossPlan[,1],
                pop@geno,crossPlan[,2],
                simParam@genMaps)
  if(simParam@gender=="no"){
    gender = rep("H",nrow(crossPlan))
  }else if(simParam@gender=="yes_rand"){
    gender = sample(c("M","F"),nrow(crossPlan),replace=TRUE)
  }else if(simParam@gender=="yes_sys"){
    gender = rep_len(c("M","F"),nrow(crossPlan))
  }else{
    stop(paste("no rules for gender type",simParam@gender))
  }
  rawPop = new("RawPop",
               nInd=nrow(crossPlan),
               nChr=pop@nChr,
               ploidy=pop@ploidy,
               nLoci=pop@nLoci,
               gender=gender,
               geno=geno)
  if(is.null(id)){
    lastId = get("LASTID",envir=.GlobalEnv)
    id = (1:rawPop@nInd) + lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  if(simParam@nTraits==0){
    gv = matrix(NA_real_,
                nrow=rawPop@nInd,
                ncol=0)
  }else{
    gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
    gv = do.call("cbind",gv)
  }
  output = new("Pop", rawPop,
               id=as.character(id),
               mother=pop@id[crossPlan[,1]],
               father=pop@id[crossPlan[,2]],
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits),
               ebv=matrix(NA_real_,
                          nrow=rawPop@nInd,
                          ncol=0))
  if(updateId){
    assign("LASTID",lastId,envir=.GlobalEnv)
  }
  return(output)
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
                     id=NULL,balance=TRUE,simParam=SIMPARAM){
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
  output = makeCross(pop=pop,crossPlan=crossPlan,id=id,
                     simParam=simParam)
  return(output)
}

#' @title Make designed crosses
#' 
#' @description
#' Makes crosses between two populations using a user supplied 
#' crossing plan.
#'
#' @param fPop an object of \code{\link{Pop-class}} for female parents.
#' @param mPop an object of \code{\link{Pop-class}} for male parents.
#' @param crossPlan a matrix with two column representing 
#' female and male parents. Either integers for the position in 
#' population or character strings for the IDs.
#' @param id optional ids to give to progeny
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#'
#' @export
makeCross2 = function(fPop,mPop,crossPlan,id=NULL,simParam=SIMPARAM){
  if(fPop@ploidy!=2){
    stop("Only works with diploids")
  }
  stopifnot(class(mPop)=="Pop",class(fPop)=="Pop")
  if(is.character(crossPlan)){ #Match by ID
    crossPlan = cbind(match(crossPlan[,1],fPop$id),
                      match(crossPlan[,2],mPop$id))
    if(any(is.na(crossPlan))){
      stop("Failed to matched supplied IDs")
    }
  }
  geno = cross2(fPop@geno,crossPlan[,1],
                mPop@geno,crossPlan[,2],
                simParam@genMaps)
  if(simParam@gender=="no"){
    gender = rep("H",nrow(crossPlan))
  }else if(simParam@gender=="yes_rand"){
    gender = sample(c("M","F"),nrow(crossPlan),replace=TRUE)
  }else if(simParam@gender=="yes_sys"){
    gender = rep_len(c("M","F"),nrow(crossPlan))
  }else{
    stop(paste("no rules for gender type",simParam@gender))
  }
  rawPop = new("RawPop",
               nInd=nrow(crossPlan),
               nChr=fPop@nChr,
               ploidy=fPop@ploidy,
               nLoci=fPop@nLoci,
               gender=gender,
               geno=geno)
  if(is.null(id)){
    lastId = get("LASTID",envir=.GlobalEnv)
    id = (1:rawPop@nInd) + lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  if(simParam@nTraits==0){
    gv = matrix(NA_real_,
                nrow=rawPop@nInd,
                ncol=0)
  }else{
    gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
    gv = do.call("cbind",gv)
  }
  output = new("Pop", rawPop,
               id=as.character(id),
               mother=fPop@id[crossPlan[,1]],
               father=mPop@id[crossPlan[,2]],
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits),
               ebv=matrix(NA_real_,
                          nrow=rawPop@nInd,
                          ncol=0))
  if(updateId){
    assign("LASTID",lastId,envir=.GlobalEnv)
  }
  return(output)
}

#' @title Make random crosses
#' 
#' @description 
#' A wrapper for \code{\link{makeCross2}} that randomly 
#' selects parental combinations for all possible combinantions between 
#' two populations.
#' 
#' @param fPop an object of \code{\link{Pop-class}} for female parents.
#' @param mPop an object of \code{\link{Pop-class}} for male parents.
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
randCross2 = function(fPop,mPop,nCrosses,nProgeny=1,
                     id=NULL,balance=TRUE,simParam=SIMPARAM){
  if(simParam@gender=="no"){
      crossPlan = sampAllComb(fPop@nInd,mPop@nInd,nCrosses)
  }else{
    female = which(fPop@gender=="F")
    if(length(female)==0){
      stop("population doesn't contain any females")
    }
    male = which(mPop@gender=="M")
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
  output = makeCross2(fPop=fPop,mPop=mPop,crossPlan=crossPlan,
                      id=id,simParam=simParam)
  return(output)
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
self = function(pop,nProgeny,id=NULL,simParam=SIMPARAM){
  stopifnot(simParam@gender=="no")
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  stopifnot(class(pop)=="Pop")
  crossPlan = cbind(rep(1:pop@nInd,each=nProgeny),
                    rep(1:pop@nInd,each=nProgeny))
  geno = cross2(pop@geno,crossPlan[,1],
                pop@geno,crossPlan[,2],
                simParam@genMaps)
  rawPop = new("RawPop",
               nInd=nrow(crossPlan),
               nChr=pop@nChr,
               ploidy=pop@ploidy,
               nLoci=pop@nLoci,
               gender=rep("H",nrow(crossPlan)),
               geno=geno)
  if(is.null(id)){
    lastId = get("LASTID",envir=.GlobalEnv)
    id = (1:rawPop@nInd) + lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  if(simParam@nTraits==0){
    gv = matrix(NA_real_,
                nrow=rawPop@nInd,
                ncol=0)
  }else{
    gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
    gv = do.call("cbind",gv)
  }
  output = new("Pop", rawPop,
               id=as.character(id),
               mother=rep(pop@mother,each=nProgeny),
               father=rep(pop@father,each=nProgeny),
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits),
               ebv=matrix(NA_real_,
                          nrow=rawPop@nInd,
                          ncol=0))
  if(updateId){
    assign("LASTID",lastId,envir=.GlobalEnv)
  }
  return(output)
}

#' @title Generates DH lines
#' 
#' @description Creates DH lines from each individual in a population. 
#' Only works when gender is "no".
#' 
#' @param pop an object of 'Pop' superclass
#' @param nDH total number of DH lines per individual
#' @param id optional id to assign to DH lines
#' @param simParam an object of 'SimParam' class
#' 
#' @export
makeDH = function(pop,nDH,id=NULL,simParam=SIMPARAM){
  stopifnot(simParam@gender=="no")
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  stopifnot(class(pop)=="Pop")
  #Should be replaced with something more efficient
  geno = createDH2(pop@geno,nDH,simParam@genMaps)
  rawPop = new("RawPop",
               nInd=as.integer(pop@nInd*nDH),
               nChr=pop@nChr,
               ploidy=pop@ploidy,
               nLoci=pop@nLoci,
               gender=rep("H",pop@nInd*nDH),
               geno=geno)
  if(is.null(id)){
    lastId = get("LASTID",envir=.GlobalEnv)
    id = (1:rawPop@nInd) + lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  if(simParam@nTraits==0){
    gv = matrix(NA_real_,
                nrow=rawPop@nInd,
                ncol=0)
  }else{
    gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
    gv = do.call("cbind",gv)
  }
  output = new("Pop", rawPop,
               id=as.character(id),
               mother=rep(pop@mother,each=nDH),
               father=rep(pop@father,each=nDH),
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits),
               ebv=matrix(NA_real_,
                          nrow=rawPop@nInd,
                          ncol=0))
  if(updateId){
    assign("LASTID",lastId,envir=.GlobalEnv)
  }
  return(output)
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
pedigreeCross = function(pedigree,founders,id=NULL,simParam=SIMPARAM){
  stopifnot(class(pedigree)=="Pedigree")
  stopifnot(class(founders)=="MapPop")
  
  sortedped = sortPed(pedigree)
  
  geno = crossPedigree(founders@geno,sortedped@father,
                sortedped@mother,
                simParam@genMaps)
  if(simParam@gender=="no"){
    gender = rep("H",sortedped@nInd)
  }else if(simParam@gender=="yes_rand"){
    gender = sample(c("M","F"),sortedped@nInd,replace=TRUE)
  }else if(simParam@gender=="yes_sys"){
    gender = rep_len(c("M","F"),sortedped@nInd)
  }else{
    stop(paste("no rules for gender type",simParam@gender))
  }
  rawPop = new("RawPop",
               nInd=sortedped@nInd,
               nChr=founders@nChr,
               ploidy=founders@ploidy,
               nLoci=founders@nLoci,
               gender=gender,
               geno=geno)
  if(is.null(id)){
    lastId = get("LASTID",envir=.GlobalEnv)
    id = (1:rawPop@nInd) + lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  if(simParam@nTraits==0){
    gv = matrix(NA_real_,
                nrow=rawPop@nInd,
                ncol=0)
  }else{
    gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
    gv = do.call("cbind",gv)
  }
  output = new("Pop", rawPop,
               id=as.character(id),
               mother=as.character(sortedped@mother),
               father=as.character(sortedped@father),
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits),
               ebv=matrix(NA_real_,
                          nrow=rawPop@nInd,
                          ncol=0))
  if(updateId){
    assign("LASTID",lastId,envir=.GlobalEnv)
  }
  return(output)
}
