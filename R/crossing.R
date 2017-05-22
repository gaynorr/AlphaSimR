#' @title Make designed crosses
#' 
#' @description
#' Makes crosses within a population using a user supplied 
#' crossing plan.
#'
#' @param pop an object of \code{\link{Pop-class}}
#' @param crossPlan a matrix with two column representing 
#' female and male parents
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
  gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
  gv = do.call("cbind",gv)
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
                          ncol=1))
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
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
randCross = function(pop,nCrosses,nProgeny=1,
                     id=NULL,simParam=SIMPARAM){
  if(simParam@gender=="no"){
    crossPlan = t(combn(pop@nInd,2))
  }else{
    female = which(pop@gender=="F")
    if(length(female)==0){
      stop("population doesn't contain any females")
    }
    male = which(pop@gender=="M")
    if(length(male)==0){
      stop("population doesn't contain any males")
    }
    crossPlan = expand.grid(female,male)
  }
  maxCrosses = nrow(crossPlan)
  if(maxCrosses>=nCrosses){
    take = sample.int(maxCrosses,nCrosses)
  }else{
    #Commented out warning to prevent confusion
    #warning("making duplicate crosses, because requested crosses exceeds unique combinations")
    take = rep_len(sample.int(maxCrosses,maxCrosses),nCrosses)
    take = take[sample.int(length(take),length(take))]
  }
  crossPlan = crossPlan[take,,drop=FALSE]
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
#' female and male parents
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
  gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
  gv = do.call("cbind",gv)
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
                          ncol=1))
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
#' @param simParam an object of \code{\link{SimParam-class}}
#' 
#' @return Returns an object of \code{\link{Pop-class}}
#' 
#' @export
randCross2 = function(fPop,mPop,nCrosses,nProgeny=1,
                     id=NULL,simParam=SIMPARAM){
  if(simParam@gender=="no"){
    crossPlan = expand.grid(1:fPop@nInd,1:mPop@nInd)
  }else{
    female = which(fPop@gender=="F")
    if(length(female)==0){
      stop("fPop doesn't contain any females")
    }
    male = which(mPop@gender=="M")
    if(length(male)==0){
      stop("mPop doesn't contain any males")
    }
    crossPlan = expand.grid(female,male)
  }
  maxCrosses = nrow(crossPlan)
  if(maxCrosses>=nCrosses){
    take = sample.int(maxCrosses,nCrosses)
  }else{
    #Commented out warning to prevent confusion
    #warning("making duplicate crosses, because requested crosses exceeds unique combinations")
    take = rep_len(sample.int(maxCrosses,maxCrosses),nCrosses)
    take = take[sample.int(length(take),length(take))]
  }
  crossPlan = crossPlan[take,,drop=FALSE]
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
  gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
  gv = do.call("cbind",gv)
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
                          ncol=1))
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
  gv = lapply(simParam@traits,getGv,pop=rawPop,w=0.5)
  gv = do.call("cbind",gv)
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
                          ncol=1))
  if(updateId){
    assign("LASTID",lastId,envir=.GlobalEnv)
  }
  return(output)
}
