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
  validObject(pop)
  stopifnot(class(pop)=="Pop")
  geno = cross2(pop@geno,crossPlan[,1],
                pop@geno,crossPlan[,2],
                simParam@genMaps)
  rawPop = new("RawPop",
               nInd = nrow(crossPlan),
               nChr = pop@nChr,
               ploidy = pop@ploidy,
               nLoci = pop@nLoci,
               gender = rep("H",nrow(crossPlan)),
               geno=geno)
  if(is.null(id)){
    lastId = get("LASTID",envir=.GlobalEnv)
    id = (1:rawPop@nInd) + lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  gv = lapply(simParam@traits,getGv,pop=rawPop,w=0)
  gv = do.call("cbind",gv)
  output = new("Pop", rawPop,
               id=as.character(id),
               mother=pop@id[crossPlan[,1]],
               father=pop@id[crossPlan[,2]],
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits))
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
#' The function does not consider gender.
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
  crossPlan = t(combn(pop@nInd,2))
  crossPlan = crossPlan[sample.int(nrow(crossPlan),nCrosses),]
  crossPlan = rep(crossPlan,nProgeny)
  output = makeCross(pop=pop,crossPlan=crossPlan,id=id,
                     simParam=simParam)
  return(output)
}

#' @title Self individuals
#' 
#' @description 
#' Creates selfed progeny from each individual in a 
#' population.
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
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  validObject(pop)
  stopifnot(class(pop)=="Pop")
  geno = cross2(pop@geno,crossPlan[,1],
                pop@geno,crossPlan[,2],
                simParam@genMaps)
  rawPop = new("RawPop",
               nInd = nrow(crossPlan),
               nChr = pop@nChr,
               ploidy = pop@ploidy,
               nLoci = pop@nLoci,
               gender = rep("H",nrow(crossPlan)),
               geno=geno)
  if(is.null(id)){
    lastId = get("LASTID",envir=.GlobalEnv)
    id = (1:rawPop@nInd) + lastId
    lastId = max(id)
    updateId = TRUE
  }else{
    updateId = FALSE
  }
  gv = lapply(simParam@traits,getGv,pop=rawPop,w=0)
  gv = do.call("cbind",gv)
  output = new("Pop", rawPop,
               id=as.character(id),
               mother=pop@id[crossPlan[,1]],
               father=pop@id[crossPlan[,2]],
               nTraits=simParam@nTraits,
               gv=gv,
               pheno=matrix(NA_real_,
                            nrow=rawPop@nInd,
                            ncol=simParam@nTraits))
  if(updateId){
    assign("LASTID",lastId,envir=.GlobalEnv)
  }
  return(output)
}

#' @title Generates DH lines
#' 
#' @description Creates DH lines from each individual in a population.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nDH total number of DH lines per individual
#' @param id optional id to assign to DH lines
#' @param simParam an object of 'SimParam' class
#' 
#' @export
makeDH = function(pop,nDH,id=NULL,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  #Should be replaced with something more efficient
  geno = cross2(pop@geno,rep(1:pop@nInd,nDH),
                pop@geno,rep(1:pop@nInd,nDH),
                simParam@genMaps)
  for(i in 1:pop@nChr){
    geno[[i]][[2]] = geno[[i]][[1]]
  }
  output = new("Pop")
  output@nInd=as.integer(pop@nInd*nDH)
  output@nChr=pop@nChr
  output@ploidy=pop@ploidy
  output@gender=rep("H",pop@nInd*nDH)
  output@geno=geno
  validObject(output)
  if(class(pop)=="PedPop"){
    if(is.null(id)) id=1:pop@nInd
    output = addPed(output,id=id,
                    par1=rep(pop@par1,nDH),
                    par2=rep(pop@par2,nDH),
                    simParam=simParam)
  }
  return(output)
}
