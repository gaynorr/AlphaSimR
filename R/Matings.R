#' @title Make random crosses
#' 
#' @description Randomly selects parental combinations for crossing. Does not consider gender.
#' 
#' @param pop an object of 'Pop' superclass
#' @param nCrosses total number of crosses to make
#' @param nProgeny number of progeny per cross
#' @param simParam an object of 'SimParam' class
#' 
#' @export
randCross = function(pop,nCrosses,nProgeny=1,simParam=SIMPARAM){
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
  output = addPed(output,id=1:output@nInd,
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
#' @param simParam an object of 'SimParam' class
#' 
#' @export
makeDH = function(pop,nDH,simParam=SIMPARAM){
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
