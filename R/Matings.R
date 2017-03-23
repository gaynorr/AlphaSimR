#' @title Make random crosses
#' 
#' @description Randomly selects parental combinations for crossing. Does not consider gender.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross
#' @param id optional id to assign to F1s
#' @param simParam an object of 'SimParam' class
#' 
#' @export
randCross = function(pop,nCrosses,nProgeny=1,
                     id=NULL,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  nInd = pop@nInd
  stopifnot(nCrosses<=(nInd*(nInd-1)/2))
  parComb = t(combn(1:nInd,2))
  parComb = parComb[sample.int(nInd*(nInd-1)/2,nCrosses),]
  femalePar = rep(parComb[,1],nProgeny)
  malePar = rep(parComb[,2],nProgeny)
  geno = cross2(pop@geno,femalePar,
                  pop@geno,malePar,
                  simParam@genMaps)
  output = new("Pop")
  output@nInd=as.integer(nCrosses*nProgeny)
  output@nChr=pop@nChr
  output@ploidy=pop@ploidy
  output@gender=rep("H",nCrosses*nProgeny)
  output@geno=geno
  validObject(output)
  if(class(pop)=="PedPop"){
    femalePar = pop@id[as.integer(femalePar)]
    malePar = pop@id[as.integer(malePar)]
  }
  if(is.null(id)){
    id=1:output@nInd
  }
  output = addPed(output,id=id,
                  par1=femalePar,par2=malePar,
                  simParam=simParam)
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

#' @title Self individuals
#' 
#' @description Creates selfed progeny from each individual in a population.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nProgeny total number of selfed progeny per individual
#' @param simParam an object of 'SimParam' class
#' 
#' @export
self = function(pop,nProgeny,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  geno = cross2(pop@geno,rep(1:pop@nInd,nProgeny),
                pop@geno,rep(1:pop@nInd,nProgeny),
                simParam@genMaps)
  output = new("Pop")
  output@nInd=as.integer(pop@nInd*nProgeny)
  output@nChr=pop@nChr
  output@ploidy=pop@ploidy
  output@gender=rep("H",pop@nInd*nProgeny)
  output@geno=geno
  validObject(output)
  return(output)
}

#' @title Test cross
#' 
#' @description Creates test crosses of all possible combinations
#' 
#' @param pop an object of 'Pop' superclass
#' @param testers an object of 'Pop' superclass
#' @param simParam an object of 'SimParam' class
#' 
#' @export
testCross = function(pop,testers,simParam=SIMPARAM){
  if(pop@ploidy!=2){
    stop("Only works with diploids")
  }
  fPar = rep(1:pop@nInd,each=testers@nInd)
  mPar = rep(1:testers@nInd,pop@nInd)
  geno = cross2(pop@geno,fPar,
                testers@geno,mPar,
                simParam@genMaps)
  if(class(pop)=="PedPop"){
    fPar = pop@id[fPar]
  }
  if(class(testers)=="PedPop"){
   mPar = testers@id[mPar] 
  }
  id = paste(fPar,mPar,sep="_")
  output = new("Pop")
  output@nInd=as.integer(pop@nInd*testers@nInd)
  output@nChr=pop@nChr
  output@ploidy=pop@ploidy
  output@gender=rep("H",pop@nInd*testers@nInd)
  output@geno=geno
  validObject(output)
  output = addPed(output,id=id,
                  par1=fPar,
                  par2=mPar,
                  simParam=simParam)
  return(output)
}

#' @title Calculate GCA
#' 
#' @description Calculate general combining ability of test crosses
#' 
#' @param pop an object of class 'PedPop' created by \code{\link{testCross}}
#' 
#' @export
calcGCA = function(pop,useGv=FALSE){
  if(useGv){
    y=pop@gv
  }else{
    y=pop@pheno
  }
  colnames(y) = paste0("Trait",1:pop@nTraits)
  output = list()
  output$inbreds=aggregate(y,list(inbred=factor(pop@par1,
                                                levels=unique(pop@par1))),
                           mean)
  output$inbreds$inbred = as.character(output$inbreds$inbred)
  output$testers=aggregate(y,list(tester=factor(pop@par2,
                                                levels=unique(pop@par2))),
                           mean)
  output$testers$tester = as.character(output$testers$tester)
  output$SCA=aggregate(y,list(inbred=factor(pop@par1,
                                            levels=unique(pop@par1)),
                              tester=factor(pop@par2,
                                            levels=unique(pop@par2))),
                       mean)
  output$SCA$inbred = as.character(output$SCA$inbred)
  output$SCA$tester = as.character(output$SCA$tester)
  return(output)
}

